#include <genetrail2/core/Category.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/WeightedGeneSetEnrichmentAnalysis.h>
#include <genetrail2/core/WGSEAPermutationTest.h>
#include <genetrail2/core/PermutationTest.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/multiprecision.h>

#include "common.h"

using namespace GeneTrail;
namespace bpo = boost::program_options;
namespace bm = boost::math;

std::string json;
bool increasing = false, absolute = false;

Params p;

GeneSet test_set;
CategoryList cat_list;
std::vector<std::string> names;
std::vector<double> values;
size_t numberOfPermutations;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	addCommonCLIArgs(desc, p);
	desc.add_options()("identifier, d", bpo::value<std::string>(&p.identifier), "A file containing identifier line by line.")(
	                   "increasing,i", bpo::value(&increasing)->zero_tokens(), "Use increasingly sorted scores. (Decreasing is default)")(
	                   "absolute,abs", bpo::value(&absolute)->zero_tokens(), "Use decreasingly sorted absolute scores.")(
	                   "permutations,per", bpo::value<size_t>(&numberOfPermutations)->default_value(1000), "Number of permutations for p-value computation.");

	if(absolute && increasing) {
		std::cerr << "ERROR: Please specify only one option to sort the file."
		          << "\n";
	}

	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),
		           vm);
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

std::shared_ptr<EnrichmentResult>
computeEnrichment(const Category& c, const std::pair<int, std::string>& genes)
{
	WeightedGeneSetEnrichmentAnalysis<big_float> gsea(names, values);

	auto result = std::make_shared<EnrichmentResult>();
	result->name = c.name();
	result->reference = c.reference();

	double RSc = gsea.computeRunningSum(c).convert_to<double>();
	;
	result->enriched = RSc > 0.0;

	result->score = RSc;

	result->hits = genes.first;
	result->info = genes.second;
	return result;
}

void computePValues(AllResults& all_results)
{
	std::vector<TestResult<double>> results;
	for(const auto& it : all_results) {
		for(const auto& jt : it.second) {
			TestResult<double> t(it.first + "\t" + jt.first, jt.second->score, jt.second->hits);
			t.enriched = jt.second->enriched;
			results.emplace_back(t);
		}
	}

	WeightedGeneSetEnrichmentAnalysis<double> wgsea(names, values);
	WGSEAPermutationTest<double> pTest(results, wgsea, names.size(),
	                                   numberOfPermutations);
	std::vector<std::pair<std::string, double>> pvalues = pTest.computePValue();
	updatePValues(all_results, pvalues);
}

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv)) {
		return -1;
	}

	if(init(test_set, cat_list, p) != 0) {
		return -1;
	}

	std::vector<std::pair<std::string, double>> sorted_test_set;

	if(absolute) {
		sorted_test_set = test_set.getAbsoluteSortedScores();
	} else {
		sorted_test_set = test_set.getSortedScores(!increasing);
	}

	names.reserve(sorted_test_set.size());
	values.reserve(sorted_test_set.size());

	for(const auto& t : sorted_test_set) {
		names.emplace_back(t.first);
		values.emplace_back(t.second);
	}

	run(test_set, cat_list, p, true);
	return 0;
}

