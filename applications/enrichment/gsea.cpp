#include <genetrail2/core/EnrichmentAlgorithm.h>
#include <genetrail2/core/GeneSet.h>

#include "common.h"

#include <boost/program_options.hpp>
#include <iostream>

using namespace GeneTrail;
namespace bpo = boost::program_options;

bool increasing = false, absolute = false;

bool parseArguments(int argc, char* argv[], Params& p)
{
	bpo::variables_map vm;
	bpo::options_description desc;

	addCommonCLIArgs(desc, p);
	desc.add_options()
		("identifier, d", bpo::value(&p.identifier), "A file containing identifier line by line.")
		("increasing,i",  bpo::value(&increasing)->zero_tokens(), "Use increasingly sorted scores. (Decreasing is default)")
		("absolute,abs",  bpo::value(&absolute)->zero_tokens(), "Use decreasingly sorted absolute scores.");

	if(absolute && increasing) {
		std::cerr << "ERROR: Please specify only one option to sort the file."
		          << "\n";
	}

	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),
		           vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "ERROR: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	return true;
}

void prepareScores(Scores& scores) {
	if(absolute) {
		std::transform(scores.scores().begin(), scores.scores().end(),
		               scores.scores().begin(),
		               static_cast<double (*)(double)>(std::abs));
	}

	scores.sortByScore(increasing ? Order::Increasing : Order::Decreasing);
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

	if(p.identifier == "") {
		prepareScores(scores);
	}

	// TODO: Improve this interface.
	auto gsea = createEnrichmentAlgorithm<KolmogorovSmirnov>(
	    p.pValueMode, scores.indices().begin(), scores.indices().end());

	run(scores, cat_list, gsea, p, true);
	return 0;
}
