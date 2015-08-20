#include <genetrail2/core/EnrichmentAlgorithm.h>
#include <genetrail2/core/GeneSet.h>

#include "common.h"

#include <boost/program_options.hpp>
#include <iostream>

namespace bpo = boost::program_options;

using namespace GeneTrail;

std::string method;

bool parseArguments(int argc, char* argv[], Params& p)
{
	bpo::variables_map vm;
	bpo::options_description desc;

	addCommonCLIArgs(desc, p);
	desc.add_options()
		("method,met", bpo::value<std::string>(&method)->default_value("no"),"Method for gene set tesing.");

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
	using _viter = Scores::ConstScoreIterator;

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
	Params p;
	if(!parseArguments(argc, argv, p)) {
		return -1;
	}

	GeneSet test_set;
	CategoryList cat_list;
	if(init(test_set, cat_list, p) != 0) {
		return -1;
	}

	Scores scores(test_set);

	auto algorithm = getAlgorithm(p.pValueMode, method, scores);

	run(scores, cat_list, algorithm, p, true);
}
