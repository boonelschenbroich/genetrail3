#include <genetrail2/core/EnrichmentAlgorithm.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/WilcoxonRankSumTest.h>
#include <genetrail2/core/OneSampleTTest.h>
#include <genetrail2/core/IndependentTTest.h>

#include "common.h"

#include <boost/program_options.hpp>
#include <iostream>

using namespace GeneTrail;
namespace bpo = boost::program_options;
namespace bm = boost::math;

std::string method;

Params p;

GeneSet test_set;
GeneSet adapted_test_set;
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

EnrichmentAlgorithmPtr getAlgorithm(const std::string& method, const Scores& scores, PValueMode mode)
{
	if(method == "one-sample-t-test") {
		return createEnrichmentAlgorithm<HTestEnrichment<OneSampleTTest<double>>>(mode, scores);
	} else if(method == "two-sample-t-test") {
		return createEnrichmentAlgorithm<HTestEnrichment<IndependentTTest<double>>>(mode, scores);
	} else if(method == "two-sample-wilcoxon") {
		return createEnrichmentAlgorithm<HTestEnrichment<WilcoxonRankSumTest<double>>>(mode, scores);
	}

	throw std::string("Invalid method " + method);
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

	Scores scores(test_set);

	auto algorithm = getAlgorithm(method, scores, p.pValueMode);

	run(scores, cat_list, algorithm, p, true);
}
