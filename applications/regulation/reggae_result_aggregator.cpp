#include <boost/program_options.hpp>

#include <genetrail2/regulation/RegulatorEffectResultAggregator.h>
#include <genetrail2/core/macros.h>

using namespace GeneTrail;

namespace bpo = boost::program_options;

std::string in_, out_, adjust_, aggregate_;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
	("input,i", bpo::value(&in_)->required(), "Path of input file.")
	("output,o", bpo::value(&out_)->required(), "Path of output file.")
	("adjust,a", bpo::value(&adjust_)->required(), "Method to adjust p-values.")
	("aggregate,g", bpo::value(&aggregate_)->required(), "Method to aggregate p-values.")
	;

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

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv)) {
		return -1;
	}
	
	RegulatorEffectResultAggregator aggregator;
	aggregator.read_results(in_);
	std::string rank_aggregator = "sum";
	aggregator.aggregateRanks(rank_aggregator);
	aggregator.aggregatePValues(aggregate_);
	aggregator.adjustPValues(adjust_);
	aggregator.write(out_);
		
	return 0;
}
