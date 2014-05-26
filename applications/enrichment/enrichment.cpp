#include <genetrail2/core/Category.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GSEAResult.h>
#include <genetrail2/core/GeneLabelPermutationTest.h>
#include <genetrail2/core/Statistic.h>

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
#include <vector>

using namespace GeneTrail;
namespace bpo = boost::program_options;
namespace bm = boost::math;
using namespace boost::multiprecision;

std::string method;
int numberOfPermutations;
Params p;

GeneSet<double> test_set;
CategoryList cat_list;
GeneLabelPermutationTest<double, std::vector<double>::iterator> pTest;

typedef std::vector<double>::iterator _viter;
std::map<std::string, std::function<double(std::vector<double>::iterator, _viter)>>
methods({{"mean",		statistic::mean<double, _viter>},
         {"median",		statistic::median<double, _viter>},
         {"sum",		statistic::sum<double, _viter>},
         {"max-mean",	statistic::max_mean<double, _viter>}});

double apply(std::vector<double>& v, std::string method){
	return methods[method](v.begin(),v.end());
}

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	addCommonCLIArgs(desc, p);
	desc.add_options()
		("method,met", bpo::value<std::string>(&method)->default_value("no"),"Method for gene set tesing.")
		("permutations,per", bpo::value<int>(&numberOfPermutations)->default_value(1000),"Number of permutations for p-value computation.");

	if (methods.find(method) != methods.end())
	{
		std::cerr << "ERROR: Given method not defined!" << std::endl;
	}

	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),vm);
		bpo::notify(vm);
	}
	catch(bpo::error& e)
	{
		std::cerr << "ERROR: " << e.what() << std::endl;
		desc.print(std::cerr);
		return false;
	}

	return true;
}

std::shared_ptr<EnrichmentResult>
computeEnrichment(const Category& c, const std::pair<int, std::string>& genes)
{
	std::cout << "INFO: Processing - " << c.name() << std::endl;
	auto result = std::make_shared<EnrichmentResult>();
	result->name = c.name();
	result->reference = c.reference();

	std::vector<double> contained_genes;
	std::vector<double> not_contained_genes;

	for(int i = 0; i < test_set.size(); ++i) {
		if(c.contains(test_set[i].first)) {
			contained_genes.push_back(test_set[i].second);
		} else {
			not_contained_genes.push_back(test_set[i].second);
		}
	}

	double z = apply(contained_genes, method);
	bool enriched = z > apply(not_contained_genes, method);
	double p =
	    enriched
	        ? pTest.computeUpperTailedPValue(genes.first, z, methods[method])
	        : pTest.computeLowerTailedPValue(genes.first, z, methods[method]);

	result->pvalue = p;

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

	if(init(test_set, cat_list, p) != 0) {
		return -1;
	}

	std::vector<double> for_permutation;
	for(auto& i : test_set.getScores()) {
		for_permutation.push_back(i.second);
	}

	pTest = GeneLabelPermutationTest<double, std::vector<double>::iterator>(
	    for_permutation.begin(), for_permutation.end(), numberOfPermutations);

	run(test_set, cat_list, p);
}
