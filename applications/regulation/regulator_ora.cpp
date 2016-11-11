#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/EntityDatabase.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/OverRepresentationAnalysis.h>
#include <genetrail2/core/BinomialTest.h>
#include <genetrail2/core/multiprecision.h>
#include <genetrail2/core/PValue.h>

#include <genetrail2/regulation/RegulatorCategoryFileReader.h>
#include <genetrail2/regulation/RegulatorEffectResult.h>

#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/prettywriter.h>

#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

using namespace GeneTrail;

namespace bpo = boost::program_options;

std::string input_, regulations_, output_, method_, correction_method_;

EntityDatabase db;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")(
	    "testset,t", bpo::value(&input_)->required(),
	    "A list of deregulated targets.")(
	    "regulations,r", bpo::value(&regulations_)->required(),
	    "A whitespace separated file containing regulator, target and scores.")(
	    "output,o", bpo::value(&output_)->required(),
	    "Output prefix for text files.")("method,m",
	                                     bpo::value(&method_)->required(),
	                                     "Method (ora, binomial-test).")(
	    "adjust,a", bpo::value(&correction_method_)->required(),
	    "P-value adjustment method.");

	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),
		           vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	return true;
}

OverRepresentationAnalysis getORA(const Category& reference, const GeneSet& test_set){
	if(method_ == "hyper"){
		return OverRepresentationAnalysis(reference, test_set.toCategory("test", &db), true);
	} else if(method_ == "fisher"){
		return OverRepresentationAnalysis(reference, test_set.toCategory("test", &db), false);	
	}
	return OverRepresentationAnalysis(reference, test_set.toCategory("test", &db));
}

std::vector<RegulatorEffectResult>
runORA(const Category& reference, const GeneSet& test_set,
       const std::map<std::string, Category>& categories)
{
	std::vector<RegulatorEffectResult> results;
	OverRepresentationAnalysis ora = getORA(reference, test_set);
	for(auto it = categories.begin(); it != categories.end(); ++it) {
		std::cout << "INFO: Processing: " << it->first << std::endl;
		if(ora.numberOfHits(it->second) < 1) {
			continue;
		}
		RegulatorEffectResult result;
		result.name = it->first;
		result.hits = ora.numberOfHits(it->second);
		result.expected_hits = ora.expectedNumberOfHits(it->second);
		result.p_value = ora.computeUpperTailedPValue(it->second);
		result.score = ora.computeScore(it->second);
		results.emplace_back(std::move(result));
	}
	return std::move(results);
}

std::vector<RegulatorEffectResult>
runBinomialTest(const Category& reference, const GeneSet& test_set,
                const std::map<std::string, Category>& categories)
{
	std::vector<RegulatorEffectResult> results;
	BinomialTest<big_float> b;
	Category test = test_set.toCategory("test", &db);
	for(auto it = categories.begin(); it != categories.end(); ++it) {
		std::cout << "INFO: Processing: " << it->first << std::endl;
		if(b.numberOfHits(test, it->second) < 1) {
			continue;
		}
		RegulatorEffectResult result;
		result.name = it->first;
		result.hits = b.numberOfHits(test, it->second);
		result.expected_hits = b.expectedNumberOfHits(reference, it->second);
		result.p_value =
		    b.computeUpperTailedPValue(reference, test, it->second).convert_to<double>();
		result.score = b.computeScore(reference, test, it->second).convert_to<double>();
		results.emplace_back(std::move(result));
	}

	return std::move(results);
}

struct get_p_value
{
	template <typename Result> double operator()(const Result& t) const
	{
		return t.corrected_p_value;
	}

	// Important: the method needs to return a mutable reference!
	template <typename Result> double& operator()(Result& t) const
	{
		return t.corrected_p_value;
	}
};

void writeJSONResults(std::vector<RegulatorEffectResult> results)
{
	rapidjson::StringBuffer sb;
	rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(sb);
	writer.StartArray();
	for(RegulatorEffectResult& res : results) {
		res.serializeJSON(writer, true);
	}
	writer.EndArray();

	std::ofstream out;
	out.open(output_);
	if(out.is_open()) {
		out << sb.GetString();
	} else {
		std::cerr << "Could not open file: " << output_ << "\n";
	}
	out.close();
}

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv)) {
		return -1;
	}

	RegulatorCategoryFileReader parser(regulations_, &db);
	parser.parse();

	Category reference = parser.getReference();
	std::map<std::string, Category> categories = parser.getCategories();

	GeneSetReader reader;
	GeneSet test_set = reader.readScoringFile(input_);

	std::vector<RegulatorEffectResult> results;
	if(method_ == "ora" || method_ == "hyper" || method_ == "fisher") {
		results = runORA(reference, test_set, categories);
	} else if(method_ == "binom") {
		results = runBinomialTest(reference, test_set, categories);
	}

	std::sort(results.begin(), results.end(), [](const auto& a, const auto& b) {
		return a.p_value < b.p_value;
	});
	
	int rank = 0;
	for(RegulatorEffectResult& r : results) {
		r.corrected_p_value = r.p_value;
		r.rank = ++rank;
	}

	if(correction_method_ != "no"){
	pvalue::adjustPValues(
	    results, get_p_value(),
	    pvalue::getCorrectionMethod(correction_method_).get());
	}
	
	writeJSONResults(results);

	return 0;
}
