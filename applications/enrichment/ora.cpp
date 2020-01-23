#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>
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
#include <fstream>
#include <chrono> 

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string reference, hypothesis, preComputedPValues;
bool usePreComputedPValues = false;

bool parseArguments(int argc, char* argv[], Params& p)
{
	bpo::variables_map vm;
	bpo::options_description desc;

	addCommonCLIArgs(desc, p);
	desc.add_options()
		("identifier, d", bpo::value(&p.identifier_), "A file containing identifier line by line.")
		("reference, r", bpo::value<std::string>(&reference)->required(), "A file containing identifier line by line.")
		("hypothesis, h", bpo::value<std::string>(&hypothesis)->default_value("two-sided"), "Null hypothesis that should be used.")
		("use-precomputed-p-values, u", bpo::value<bool>(&usePreComputedPValues)->zero_tokens(), "Use precomputed p-values.")
		("precomputed-p-values, p", bpo::value<std::string>(&preComputedPValues), "Precomputed p-values.");

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

    if (!usePreComputedPValues) {
		auto enrichmentAlgorithm = createEnrichmentAlgorithm<Ora>(p.pValueMode, reference_set.toCategory(db, "reference"), test_set.toCategory(db, "test"), hypothesis_);
		Scores scores(test_set, db);
		run(scores, cat_list, enrichmentAlgorithm, p, true);
	} else {
		//auto start = std::chrono::high_resolution_clock::now();
		DenseMatrixReader reader;
		std::ifstream file(preComputedPValues);
		if(!file) {
			std::cerr << "Could not open " << preComputedPValues << " for reading." << std::endl;
			return -1;
		}
		//This is needed since our matrix does not contain any row/col names
		unsigned int opts = DenseMatrixReader::NO_OPTIONS;
		DenseMatrix p_values = reader.read(file, opts);
		file.close();
		//auto finish = std::chrono::high_resolution_clock::now();
		//std::chrono::duration<double> elapsed = finish - start;
		//std::cout << "Elapsed time (Loading Matrix): " << elapsed.count() << " s\n";

		p.verbose = false;
                //start = std::chrono::high_resolution_clock::now();
		auto enrichmentAlgorithm = createEnrichmentAlgorithm<PreprocessedORA>(p.pValueMode, reference_set.toCategory(db, "reference"), test_set.toCategory(db, "test"), hypothesis_, p_values, p.reducedOutput);
		Scores scores(test_set, db);
		run(scores, cat_list, enrichmentAlgorithm, p, true);
		//finish = std::chrono::high_resolution_clock::now();
		//elapsed = finish - start;
		//std::cout << "Elapsed time (Enrichment analysis): " << elapsed.count() << " s\n";
	}

	return 0;
}
