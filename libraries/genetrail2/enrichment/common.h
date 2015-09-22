#ifndef GT2_APPLICATIONS_ENRICHMENT_COMMON_H
#define GT2_APPLICATIONS_ENRICHMENT_COMMON_H

#include "CommandLineInterface.h"

#include <genetrail2/core/EnrichmentAlgorithm.h>
#include <genetrail2/core/MatrixHTest.h>
#include <genetrail2/core/macros.h>
#include <genetrail2/core/PValue.h>

#include <boost/program_options.hpp>

#include <list>
#include <utility>
#include <string>
#include <memory>
#include <random>

namespace bpo = boost::program_options;

using namespace GeneTrail;

namespace GeneTrail
{
	class GeneSet;
	class EnrichmentResult;
	class FilePath;
	class DirectoryPath;
	class Params;

	using EnrichmentAlgorithmPtr = std::unique_ptr<EnrichmentAlgorithm>;
}

typedef std::map<std::string, std::shared_ptr<EnrichmentResult>> Results;
typedef std::map<std::string, Results> AllResults;
typedef std::vector<std::pair<std::string, double>> PValueList;

using CategoryList = std::list<std::pair<std::string, std::string>>;


/**
 * This function initializes the needed attributes.
 *
 * @param test_set Test set to be filled
 * @param cat_list CategoryList to be filled
 * @param p Parameter object
 * @return -1 if an error occurred and 0 if not
 */
GT2_EXPORT int init(GeneSet& test_set, CategoryList& cat_list, const Params& p);

/**
 * This function runs the entire pipeline.
 *
 * @param test_set Test set for which the computation should be started
 * @param cat_list List of categories for the computation
 * @param p
 */
GT2_EXPORT void run(Scores& test_set, CategoryList& cat_list, EnrichmentAlgorithmPtr& algorithm, const Params& p, bool computePValues);

#endif // GT2_APPLICATIONS_ENRICHMENT_COMMON_H
