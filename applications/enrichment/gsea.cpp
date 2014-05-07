#include <genetrail2/core/Category.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/GeneSetEnrichmentAnalysis.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/GeneSet.h>

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

using namespace GeneTrail;
namespace bpo = boost::program_options;
namespace bm = boost::math;
using namespace boost::multiprecision;

std::string json;
bool increasing = false, absolute = false;

Params p;

GeneSet<double> test_set;
CategoryList cat_list;
AllResults name_to_cat_results;
GeneSetEnrichmentAnalysis<cpp_dec_float_50, int64_t> gsea;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	addCommonCLIArgs(desc, p);
	desc.add_options()
		("identifier, d", bpo::value<std::string>(&p.identifier), "A file containing identifier line by line.")
		("increasing,i", bpo::value(&increasing)->zero_tokens(), "Use increasingly sorted scores. (Decreasing is default)")
		("absolute,abs", bpo::value(&absolute)->zero_tokens(), "Use decreasingly sorted absolute scores.")
		("json,j", bpo::value<std::string>(&json), "Output filename of json file.");

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

EnrichmentResult computeEnrichment(const Category& c, std::pair<int,std::string> genes)
{
	std::cout << "INFO: Processing " << c.name() << std::endl;
	EnrichmentResult result;
	result.name = c.name();
	result.reference = c.reference();

	std::vector<std::string> test;
	if(p.scores != ""){
		if(absolute){
			test = test_set.getIdentifier(test_set.getAbsoluteSortedScores());
		}else{
			if(increasing){
				test = test_set.getIdentifier(test_set.getIncreasinglySortedScores());
			}else{
				test = test_set.getIdentifier(test_set.getDecreasinglySortedScores());
			}
		}
	}else{
		test = test_set.getIdentifier(test_set.getScores());
	}

	auto enr = gsea.computePValue(c, test);
	result.pvalue = enr.first.convert_to<double>();

	int abs_max = -1;
	bool enriched;
	for(auto p : enr.second){
		if(std::abs(p.second) > abs_max){
			abs_max = std::abs(p.second);
			enriched = p.second >= 0.0;
		}
	}
	result.enriched = enriched;
	result.hits = genes.first;
	result.info = genes.second;
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

	run(test_set, cat_list, name_to_cat_results, p);
}
