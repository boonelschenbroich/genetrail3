#include <genetrail2/core/Category.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/GeneSetEnrichmentAnalysis.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GSEAResult.h>

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
#include <stdlib.h>
#include <memory>

using namespace GeneTrail;
namespace bpo = boost::program_options;
namespace bm = boost::math;
using namespace boost::multiprecision;

std::string json;
bool increasing = false, absolute = false;

Params p;

GeneSet<double> test_set;
CategoryList cat_list;
GeneSetEnrichmentAnalysis<cpp_dec_float_50, int64_t> gsea;
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
		std::cerr << "Error: Please specify only one option to sort the file." << "\n";
	}

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

std::shared_ptr<EnrichmentResult> computeEnrichment(const Category& c, std::pair<int,std::string> genes)
{
	std::cout << "INFO: Processing - " << c.name() << std::endl;
	std::shared_ptr<GSEAResult> result(new GSEAResult());
	result->name = c.name();
	result->reference = c.reference();

	auto enr = gsea.computePValue(c, identifierOfTestSet);
	result->pvalue = enr.first;
	result->points = enr.second;
	int abs_max = -1;
	bool enriched;
	for(auto p : enr.second){
		if(std::abs(p.second) > abs_max){
			abs_max = std::abs(p.second);
			enriched = p.second >= 0.0;
		}
	}

	result->enriched = enriched;
	result->hits = genes.first;
	result->info = genes.second;
	return result;
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

	identifierOfTestSet = getSortedIdentifier(test_set, p,  absolute, increasing);

	run(test_set, cat_list, p);
}
