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

typedef Scores::ConstScoreIterator _viter;


bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	addCommonCLIArgs(desc, p);
	desc.add_options()
		("method,met", bpo::value<std::string>(&method)->default_value("no"),"Method for gene set tesing.")
		("permutations,per", bpo::value<int>(&numberOfPermutations)->default_value(1000),"Number of permutations for p-value computation.");

	//TODO: Validate arguments
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

	checkCLIArgs(p);

	return true;
}

EnrichmentAlgorithmPtr getAlgorithm(PValueMode mode, const std::string& method, const Scores& scores) {
	if(method == "mean") {
		return createEnrichmentAlgorithm<StatisticsEnrichment>(mode, statistic::mean<double, _viter>, scores);
	} else if(method == "median") {
		return createEnrichmentAlgorithm<StatisticsEnrichment>(mode, statistic::median<double, _viter>, scores);
	} else if(method == "sum") {
		return createEnrichmentAlgorithm<StatisticsEnrichment>(mode, statistic::sum<double, _viter>, scores);
	} else if(method == "max-mean") {
		return createEnrichmentAlgorithm<StatisticsEnrichment>(mode, statistic::max_mean<double, _viter>, scores);
	}

	//TODO: Throw the proper exception
	throw std::string("Unknown method: " + method);
}

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv)) {
		return -1;
	}

	if(init(test_set, cat_list, p) != 0) {
		return -1;
	}

	Scores scores(test_set);

	auto algorithm = getAlgorithm(p.pValueMode, method, scores);

	run(scores, cat_list, algorithm, p, true);
}
