#include <genetrail2/core/Category.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/HTest.h>
#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/OneSampleWilcoxonSignedRankTest.h>
#include <genetrail2/core/WilcoxonRankSumTest.h>
#include <genetrail2/core/OneSampleTTest.h>
#include <genetrail2/core/IndependentTTest.h>
#include <genetrail2/core/multiprecision.h>
#include <genetrail2/core/compat.h>

#include "common.h"

#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <iostream>
#include <fstream>
#include <cstdint>
#include <utility>
#include <map>
#include <cstdlib>
#include <vector>

using namespace GeneTrail;
namespace bpo = boost::program_options;
namespace bm = boost::math;

std::string method;

Params p;

GeneSet test_set;
CategoryList cat_list;
AllResults name_to_cat_results;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	addCommonCLIArgs(desc, p);
	desc.add_options()
		("method,met", bpo::value<std::string>(&method)->required(), "Method for gene set tesing.");

	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),vm);
		bpo::notify(vm);
	}
	catch(bpo::error& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	return true;
}

std::unique_ptr<EnrichmentResult> computeEnrichment(const Category& c, const std::pair<int,std::string>& genes)
{
	auto result = std::make_unique<EnrichmentResult>();
	result->name = c.name();
	result->reference = c.reference();

	std::vector<double> all_genes;
	all_genes.resize(test_set.size());

	std::vector<double> contained_genes;
	std::vector<double> not_contained_genes;

	for(int i = 0; i < test_set.size(); ++i) {
		all_genes[i] = test_set[i].second;
		if(c.contains(test_set[i].first)) {
			contained_genes.push_back(test_set[i].second);
		}else{
			not_contained_genes.push_back(test_set[i].second);
		}
	}

	result->hits = genes.first;
	result->info = genes.second;
	result->score = 0;
	result->pvalue = 1;

	if(contained_genes.size() <= 1) {
		std::cerr << "WARNING: Could not compute a p-value for category " + c.name() + ". The number of found genes is too small." << std::endl;
		return result;
	}

	if(method == "one-sample-t-test") {
		auto mean = statistic::mean<double, std::vector<double>::iterator>(all_genes.begin(), all_genes.end());
		OneSampleTTest<double> ttest(1e-4, mean);
		result->score = HTest::test(ttest, contained_genes.begin(), contained_genes.end());
		if(result->score < 0){
			result->pvalue = HTest::lowerTailedPValue(ttest, result->score);
		}else{
			result->pvalue = HTest::upperTailedPValue(ttest, result->score);
		}
		result->enriched = result->score > 0;

		return result;
	}

	if(not_contained_genes.size() <= 1){
		std::cerr << "WARNING: Could not compute a p-value for category " + c.name() + ". The number of found genes is too small." << std::endl;
		return result;
	}

	if(method == "two-sample-wilcoxon"){
		WilcoxonRankSumTest<big_float> wilcox(1e-4);
		auto score = HTest::test(wilcox, contained_genes.begin(), contained_genes.end(), not_contained_genes.begin(), not_contained_genes.end());
		result->score = score.convert_to<double>();;
		if(wilcox.enriched()){
			result->pvalue = HTest::upperTailedPValue(wilcox, score);
		}else{
			result->pvalue = HTest::lowerTailedPValue(wilcox, score);
		}
		result->enriched = wilcox.enriched();
	}else if(method == "two-sample-t-test"){
		IndependentTTest<double> ttest(1e-4);
		result->score = HTest::test(ttest, contained_genes.begin(), contained_genes.end(), not_contained_genes.begin(), not_contained_genes.end());
		if(result->score < 0){
			result->pvalue = HTest::lowerTailedPValue(ttest, result->score);
		}else{
			result->pvalue = HTest::upperTailedPValue(ttest, result->score);
		}
		result->enriched = result->score > 0;
	}

	return result;
}

void computePValues(AllResults& results)
{}

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv)) {
		return -1;
	}

	if(init(test_set,cat_list,p) != 0)
	{
		return -1;
	}

	run(test_set, cat_list, p);
}
