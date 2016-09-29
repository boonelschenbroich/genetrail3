#include <genetrail2/regulation/Regulator.h>
#include <genetrail2/regulation/RegulationScoreFileReader.h>
#include <genetrail2/regulation/RegulatorPermutationTest.h>
#include <genetrail2/core/Statistic.h>

#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>
#include <string>

using namespace GeneTrail;

namespace bpo = boost::program_options;

std::string scores_, out_;
size_t permutations_, randomSeed_;
bool useAbsoluteValues_;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")(
	    "scores,s", bpo::value(&scores_)->required(),
	    "A whitespace separated file containing regulator, target and scores.")(
	    "permutations,p", bpo::value(&permutations_)->required(),
	    "Number of permutations used to compute a p-value.")(
	    "seed,e", bpo::value(&randomSeed_)->required(),
	    "Random seed used to perform permutations.")(
	    "abs,a", bpo::value(&useAbsoluteValues_)->default_value(false)->zero_tokens(),
	    "Use absolute values.")(
	    "output,o", bpo::value(&out_)->required(),
	    "Output prefix for text files.");

	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),
		           vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	return true;
}

void write(const std::vector<Regulator<double>>& regulators)
{
	std::ofstream out;
	out.open(out_);
	if(out.is_open()) {
		for(size_t i = 0; i < regulators.size(); ++i) {
			out << regulators[i].name() << "\t"
				<< regulators[i].original_score() << "\t"
				<< regulators[i].numberOfTargets() << "\t"
			    << regulators[i].estimateUpperTailedPValue() << "\n";
		}
	} else {
		std::cerr << "Could not open file: " << out_ << "\n";
	}
	out.close();
}

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv)) {
		return -1;
	}

	RegulationScoreFileReader<double> parser(scores_, permutations_, useAbsoluteValues_);
	std::vector<Regulator<double>> regulators =
	    parser.parse(statistic::mean<double, std::vector<double>::iterator>);
	std::vector<double> scores = parser.getScores();

	RegulatorPermutationTest<double> test(regulators, scores, permutations_,
	                                      randomSeed_);
	regulators =
	    test.execute(statistic::mean<double, std::vector<double>::iterator>);

	write(regulators);

	return 0;
}
