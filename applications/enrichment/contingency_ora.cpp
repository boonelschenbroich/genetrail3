#include <iostream>
#include <fstream>
#include <tuple>

#include <boost/program_options.hpp>

#include <genetrail2/enrichment/ContingencyMatrixParser.h>
#include <genetrail2/enrichment/ContingencyEnrichmentResult.h>
#include <genetrail2/enrichment/ContingencyORA.h>

#include <genetrail2/core/multiprecision.h>
#include <genetrail2/core/PValue.h>

#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/prettywriter.h>

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string input_file = "", output_file, adjustment_method_ = "";

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
	    ("adjust,a", bpo::value(&adjustment_method_)->required()->default_value("benjamini-yekutieli"), "Method for multiple testing correction. (default: benjamini-yekutieli)")
		("input,i", bpo::value<std::string>(&input_file)->required(), "Input file (AB, notAB, AnotB, notAnotB).")
		("output,o", bpo::value<std::string>(&output_file)->required(), "Name of the output file.");

	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);
		bpo::notify(vm);
	}
	catch(bpo::error& e)
	{
		std::cerr << "ERROR: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	return true;
}

struct get_pvalue
{
	template <typename T> double operator()(const T& t) const
	{
		return t.p_value;
	}

	// Important: the method needs to return a mutable reference!
	template <typename T> double& operator()(T& t) const
	{
		return t.p_value;
	}
};

void writeJSONResults(std::vector<ContingencyEnrichmentResult>& results)
{
	std::sort(
	    results.begin(), results.end(),
	    [](const ContingencyEnrichmentResult& lhs, const ContingencyEnrichmentResult& rhs) {
		    return lhs.p_value < rhs.p_value;
		});

	rapidjson::StringBuffer sb;
	rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(sb);

	writer.StartArray();
	for(auto& res : results) {
		res.serializeJSON(writer);
	}
	writer.EndArray();

	std::ofstream out;
	out.open(output_file);
	if(out.is_open()) {
		out << sb.GetString();
	} else {
		std::cerr << "Could not open file: " << output_file << "\n";
	}
	out.close();
}

int main(int argc, char* argv[])
{ 
  	if(!parseArguments(argc, argv))
	{
		return -1;
	}

	std::cout << "INFO: Parsing contingency matrix ..." << std::endl;
	ContingencyMatrixParser<uint64_t> parser(input_file);
	auto ctables = parser.getContingencyTables();

	std::cout << "INFO: Setting up ORA ..." << std::endl;
	ContingencyORA<uint64_t, long double> ora (ctables);
	auto results = ora.run(true);

	std::cout << "INFO: Adjusting p-values ..." << std::endl;
	results = pvalue::adjustPValues(results, get_pvalue(), pvalue::getCorrectionMethod(adjustment_method_).get());

	std::cout << "INFO: Writing results ..." << std::endl;
	writeJSONResults(results);

    return 0;
}
