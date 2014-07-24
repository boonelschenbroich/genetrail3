#include <genetrail2/core/Category.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/HTest.h>
#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/OneSampleWilcoxonSignedRankTest.h>
#include <genetrail2/core/OneSampleTTest.h>

#include "common.h"

#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

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
using namespace boost::multiprecision;

std::string method;

Params p;

GeneSet<double> test_set;
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

std::shared_ptr<EnrichmentResult> computeEnrichment(const Category& c, const std::pair<int,std::string>& genes)
{
	auto result = std::make_shared<EnrichmentResult>();
	result->name = c.name();
	result->reference = c.reference();

	std::vector<double> all_genes;
	all_genes.resize(test_set.size());

	std::vector<double> contained_genes;

	for(int i = 0; i < test_set.size(); ++i) {
		all_genes.push_back(test_set[i].second);
		if(c.contains(test_set[i].first)) {
			contained_genes.push_back(test_set[i].second);
		}
	}

	if(method == "wilcoxon"){
		auto median = statistic::median<double, std::vector<double>::iterator>(all_genes.begin(), all_genes.end());
		OneSampleWilcoxonSignedRankTest<double, std::vector<double>::iterator> wilcox(1e-4, median);
		result->score = HTest::test(wilcox, contained_genes.begin(), contained_genes.end());
		if(!wilcox.enriched()){
			result->pvalue = HTest::lowerTailedPValue(wilcox, result->score);
		}else{
			result->pvalue = HTest::upperTailedPValue(wilcox, result->score);
		}
		result->enriched = wilcox.enriched();
	}else if(method == "t-test"){
		auto mean = statistic::mean<double, std::vector<double>::iterator>(all_genes.begin(), all_genes.end());
		OneSampleTTest<double, std::vector<double>::iterator> ttest(1e-4, mean);
		result->score = HTest::test(ttest, contained_genes.begin(), contained_genes.end());
		if(result->score < 0.0){
			result->pvalue = HTest::lowerTailedPValue(ttest, result->score);
		}else{
			result->pvalue = HTest::upperTailedPValue(ttest, result->score);
		}
		result->enriched = result->score > 0.0;
	}

	result->hits = genes.first;
	result->info = genes.second;
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
