#include <genetrail2/core/Category.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/GeneSetEnrichmentAnalysis.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/GeneSet.h>
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
#include <memory>

using namespace GeneTrail;
namespace bpo = boost::program_options;
namespace bm = boost::math;

std::string json;
bool increasing = false, absolute = false;

Params p;

GeneSet test_set;
GeneSet adapted_test_set;
CategoryList cat_list;
GeneSetEnrichmentAnalysis<big_float, int64_t> gsea;
std::vector<std::string> identifierOfTestSet;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	addCommonCLIArgs(desc, p);
	desc.add_options()
		("identifier, d", bpo::value<std::string>(&p.identifier), "A file containing identifier line by line.")
		("increasing,i", bpo::value(&increasing)->zero_tokens(), "Use increasingly sorted scores. (Decreasing is default)")
		("absolute,abs", bpo::value(&absolute)->zero_tokens(), "Use decreasingly sorted absolute scores.");

	if(absolute && increasing) {
		std::cerr << "ERROR: Please specify only one option to sort the file." << "\n";
	}

	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),vm);
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

GeneSet adapt_all_gene_sets(const Category& all_genes_of_database)
{
	adapted_test_set = adapt_gene_set(test_set, all_genes_of_database);
	identifierOfTestSet = getSortedIdentifier(adapted_test_set, p,  absolute, increasing);
	return adapted_test_set;
}

int intersectionSize(const Category& category, const std::vector<std::string>& testSet){
	int n = 0;
	for(auto gene : testSet ){
		if(category.contains(gene)){
			++n;
		}
	}
	return n;
}

std::unique_ptr<EnrichmentResult> computeEnrichment(const Category& c, const std::pair<int,std::string>& genes)
{
	auto result = std::make_unique<EnrichmentResult>(c);
	
	int RSc = gsea.computeRunningSum(c, identifierOfTestSet);
	result->enriched = RSc > 0;
	if(result->enriched){
		result->pvalue = gsea.computeRightPValue(identifierOfTestSet.size(), intersectionSize(c, identifierOfTestSet), RSc);
	}else{
		result->pvalue = gsea.computeLeftPValue(identifierOfTestSet.size(), intersectionSize(c, identifierOfTestSet), RSc);
	}

	result->hits = genes.first;
	result->info = genes.second;
	return result;
}

void computePValues(AllResults& results)
{
}

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv)) {
		return -1;
	}

	if(init(test_set,cat_list,p) != 0)
	{
		return -1;
	}

	adapted_test_set = test_set;
	identifierOfTestSet = getSortedIdentifier(test_set, p,  absolute, increasing);

	run(test_set, cat_list, p);
	return 0;
}
