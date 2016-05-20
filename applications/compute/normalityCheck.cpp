#include <iostream>
#include <cstdlib>
#include <tuple>
#include <vector>
#include <functional>
#include <utility>

#include <boost/program_options.hpp>

#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/NormalityTest.h>

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string scores, distribution;

typedef std::function<bool(std::vector<double>)> ftype;

const double TOLERANCE = 0.05;

ftype gaussian = [](std::vector<double> v) {
	return NormalityTest::testForNormalDistribution(v.begin(), v.end(),
	                                                TOLERANCE);
};

ftype standardNormal = [](std::vector<double> v) {
	return NormalityTest::testForStandardNormalDistribution(v.begin(), v.end(),
	                                                        TOLERANCE);
};

std::map<std::string, ftype> distributions{{"gaussian", gaussian},
                                           {"standard-normal", standardNormal}};

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")(
	    "scores,s", bpo::value<std::string>(&scores)->required(),
	    "Name of the scoring file.")(
	    "distribution,d", bpo::value<std::string>(&distribution)->required(),
	    "Check if the scores follow the given distribution(gaussian, "
	    "standard-normal).");

	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),
		           vm);
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

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv)) {
		return -1;
	}

	GeneSetReader reader;
	GeneSet gene_set = reader.readScoringFile(scores);
	std::vector<std::pair<std::string, double>> genes = gene_set.getScores();
	std::vector<double> scores;
	scores.resize(gene_set.size());
	for(size_t i = 0; i < gene_set.size(); ++i) {
		scores[i] = genes[i].second;
	}

	return (distributions[distribution](scores)) ? 0 : 1;
}
