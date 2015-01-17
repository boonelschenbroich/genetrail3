#include <genetrail2/core/Category.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/PermutationTest.h>
#include <genetrail2/core/Statistic.h>
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
#include <vector>

using namespace GeneTrail;
namespace bpo = boost::program_options;
namespace bm = boost::math;

std::string method;
int numberOfPermutations;
Params p;

GeneSet test_set;
CategoryList cat_list;

typedef std::vector<double>::iterator _viter;
std::map<std::string, std::function<double(_viter, _viter)>>
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

std::unique_ptr<EnrichmentResult>
computeEnrichment(const Category& c, const std::pair<int, std::string>& genes)
{
	auto result = std::make_unique<EnrichmentResult>();
	result->name = c.name();
	result->reference = c.reference();

	std::vector<double> contained_genes;
	std::vector<double> all_genes;
	all_genes.resize(test_set.size());

	for(int i = 0; i < test_set.size(); ++i) {
		all_genes[i] = test_set[i].second;
		if(c.contains(test_set[i].first)) {
			contained_genes.push_back(test_set[i].second);
		}
	}

	result->score = apply(contained_genes, method);
	if(method == "mean"){
		result->enriched = result->score >= apply(all_genes, method);
	}else if(method == "median"){
		result->enriched = result->score >= apply(all_genes, method);
	}else {
		result->enriched = result->score >= 0.0;
	}
	result->hits = genes.first;
	result->info = genes.second;
	return result;
}

void computePValues(AllResults& all_results)
{
	std::vector<TestResult<double>> results;
	std::vector<double> scores;
	scores.reserve(test_set.size());
	for(int i = 0; i < test_set.size(); ++i) {
		scores.push_back(test_set[i].second);
	}

	for(const auto& it : all_results) {
		for(const auto& jt : it.second) {
			TestResult<double> t(it.first + "\t" + jt.first, jt.second->score, jt.second->hits);
			t.enriched = jt.second->enriched;
			results.emplace_back(t);
		}
	}
	PermutationTest<double> pTest(results, scores, numberOfPermutations);
	std::vector<std::pair<std::string, double>> pvalues = pTest.computePValue(methods[method]);
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

	run(test_set, cat_list, p, true);
}
