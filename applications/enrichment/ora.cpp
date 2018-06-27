#include <genetrail2/core/EntityDatabase.h>
#include <genetrail2/core/Exception.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/OverRepresentationAnalysis.h>

#include <genetrail2/enrichment/common.h>
#include <genetrail2/enrichment/EnrichmentAlgorithm.h>
#include <genetrail2/enrichment/SetLevelStatistics.h>
#include <genetrail2/enrichment/Parameters.h>

#include <boost/program_options.hpp>
#include <iostream>

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string reference, hypothesis;

bool parseArguments(int argc, char* argv[], Params& p)
{
	bpo::variables_map vm;
	bpo::options_description desc;

	addCommonCLIArgs(desc, p);
	desc.add_options()
		("identifier, d", bpo::value(&p.identifier_), "A file containing identifier line by line.")
		("reference, r", bpo::value<std::string>(&reference)->required(), "A file containing identifier line by line.")
		("hypothesis, h", bpo::value<std::string>(&hypothesis)->default_value("two-sided"), "Null hypothesis that should be used.");

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

	return checkCLIArgs(p);
}

int main(int argc, char* argv[])
{
	Params p;
	if(!parseArguments(argc, argv, p)) {
		return -1;
	}

	GeneSet test_set;
	GeneSet reference_set;
	CategoryList cat_list;

	if(init(test_set, cat_list, p) != 0)
	{
		return -1;
	}

	GeneSetReader reader;
	try
	{
		reference_set = reader.readGeneList(reference);
	}
	catch(IOError& exn)
	{
		std::cerr << "ERROR: Failed to read reference set. Reason: " << exn.what() << std::endl;
		return -1;
	}

	auto db = std::make_shared<EntityDatabase>();
	NullHypothesis hypothesis_;
	if (hypothesis == "upper-tailed") {
		hypothesis_ = NullHypothesis::UPPER_TAILED;
	} else if (hypothesis == "lower-tailed") {
		hypothesis_ = NullHypothesis::LOWER_TAILED;
	} else {
		hypothesis_ = NullHypothesis::TWO_SIDED;
	}

	auto db = std::make_shared<EntityDatabase>();
	auto enrichmentAlgorithm = createEnrichmentAlgorithm<Ora>(p.pValueMode, reference_set.toCategory(db, "reference"), test_set.toCategory(db, "test"), hypothesis_);

	Scores scores(test_set, db);
	run(scores, cat_list, enrichmentAlgorithm, p, true);

	return 0;
}
