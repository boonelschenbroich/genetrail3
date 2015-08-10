#ifndef GT2_APPLICATIONS_ENRICHMENT_COMMON_H
#define GT2_APPLICATIONS_ENRICHMENT_COMMON_H

#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

#include <genetrail2/core/Exception.h>
#include <genetrail2/core/EnrichmentAlgorithm.h>
#include <genetrail2/core/EnrichmentResult.h>
#include <genetrail2/core/Category.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/PValue.h>

#include <list>
#include <utility>
#include <fstream>
#include <string>
#include <iostream>
#include <memory>
#include <cassert>

using namespace GeneTrail;
namespace bpo = boost::program_options;

namespace GeneTrail
{
	class GeneSet;
	class EnrichmentAlgorithm;

	using EnrichmentAlgorithmPtr = std::unique_ptr<EnrichmentAlgorithm>;
}

typedef std::map<std::string, std::shared_ptr<EnrichmentResult>> Results;
typedef std::map<std::string, Results> AllResults;
typedef std::vector<std::pair<std::string, double>> PValueList;

struct Params
{
	std::string algorithm;
	double significance;
	std::string categories;
	std::string scores;
	std::string identifier;
	size_t minimum;
	size_t maximum;
	std::string out;
	std::string adjustment;
	bool runSeparately;
	PValueMode pValueMode;
	size_t numPermutations;
	std::string dataMatrixPath;
	std::string groups;
	std::string scoringMethod;
};

using CategoryList = std::list<std::pair<std::string, std::string>>;

namespace GeneTrail {
	void validate(boost::any&, const std::vector<std::string>&, PValueMode*,
	              int);
}

/**
 * This function adds common arguments to the BOOST commandline parser.
 *
 * @param desc BOOST program-options description
 * @param p Parameter object
 */
void addCommonCLIArgs(bpo::options_description& desc, Params& p);

/**
 * This function initializes the needed attributes.
 *
 * @param test_set Test set to be filled
 * @param cat_list CategoryList to be filled
 * @param p Parameter object
 * @return -1 if an error occurred and 0 if not
 */
int init(GeneSet& test_set, CategoryList& cat_list, const Params& p);

/**
 * This function returns the sorted identifier of a GeneSet.
 *
 * @param test_set GeneSet
 * @param alsolute
 * @param increasing
 * @return Sorted list of identifier.
 */
std::vector<std::string> getSortedIdentifier(GeneSet& test_set, const Params& p, bool absolute, bool increasing);

/**
 * This function runs the entire pipeline.
 *
 * @param test_set Test set for which the computation should be started
 * @param cat_list List of categories for the computation
 * @param p
 */
void run(GeneSet& test_set, CategoryList& cat_list, EnrichmentAlgorithmPtr& algorithm, const Params& p, bool computePValues);

#endif // GT2_APPLICATIONS_ENRICHMENT_COMMON_H
