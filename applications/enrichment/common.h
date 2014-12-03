#ifndef GT2_APPLICATIONS_ENRICHMENT_COMMON_H
#define GT2_APPLICATIONS_ENRICHMENT_COMMON_H

#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

#include <genetrail2/core/Exception.h>
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
}

typedef std::map<std::string, std::shared_ptr<EnrichmentResult>> Results;
typedef std::map<std::string, Results> AllResults;
typedef std::vector<std::pair<std::string, double>> PValueList;

struct Params
{
	double significance;
	std::string categories;
	std::string scores;
	std::string identifier;
	int minimum;
	int maximum;
	std::string out;
	std::string adjustment;
	bool runSeparately;
};

/**
 * This function adds common arguments to the BOOST commandline parser.
 *
 * @param desc BOOST program-options description
 * @param p Parameter object
 */
void addCommonCLIArgs(bpo::options_description& desc, Params& p);

/**
 * This function parses the given input file and saves the information as GeneSet object.
 *
 * @param test_set GeneSet object to be filled by this function
 * @param p Parameter object
 */
void readTestSet(GeneSet& test_set, const Params& p);

typedef std::list<std::pair<std::string, std::string>> CategoryList;

/**
 * This function parses the given category file and returns the contained information in a list.
 *
 * @param catfile_list File containing the categories for which pvalues should be computed.
 * @param single_cati TODO
 * @return List of (category name, path to category file) pairs.
 */
CategoryList getCategoryList(const std::string& catfile_list, const std::string& single_cat);


/**
 * This function converts a AllResults obeject into a PValueList object.
 *
 * @param results AllResults& object to be converted
 * @return List of (category name,  pvalue) pairs.
 */
PValueList resultVector(const AllResults& results);

/**
 * This function converts a Results obeject into a PValueList object.
 *
 * @param results Results& object to be converted
 * @return List of (category name,  pvalue) pairs.
 */
PValueList resultVector(const Results& results);

/**
 *
 * @param c Category that should be processed
 * @param test_set TestSet for which pvalue should be computed based on the category
 * @param p Parameter object
 * @return
 */
std::pair<bool, std::pair<int, std::string>> processCategory(Category& c, GeneSet& test_set, const Params& p);

/**
 * Writes the results for each database in a separate file.
 *
 * @param output_dir Path and prefix of the output files.
 * @param all_results AllResults object containing all results.
 */
void writeFiles(const std::string& output_dir, const AllResults& all_results);

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
 * This function computes all results for a given Category.
 *
 * @param c Category
 * @param genes
 * @return Shared pointer to an Enrichment result object.
 */
std::shared_ptr<EnrichmentResult> computeEnrichment(const Category& c, const std::pair<int, std::string>& genes);

/**
 * This function updates all pvalues in a Results object.
 *
 * @params results Results object
 * @param pvalues This new pvalues
 */
void updatePValues(Results& results, const PValueList& pvalues);

/**
 * This function updates all pvalues in a AllResults object.
 *
 * @params results Results object
 * @param pvalues This new pvalues
 */
void updatePValues(AllResults& results, const PValueList& pvalues);

/**
 * This function adjusts all pvalues in a AllResults object combined.
 *
 * @params results Results object
 * @param pvalues This new pvalues
 */
void adjustCombined(AllResults& results, const Params& p);

/**
 * This function adjusts all pvalues in a AllResults object.
 * All databases will be separately adjusted.
 *
 * @params results Results object
 * @param pvalues This new pvalues
 */
void adjustSeparately(AllResults& results, const Params& p);

/**
 *
 */
void computePValues(AllResults& results);

/**
 * This function runs the entire pipeline.
 *
 * @param test_set Test set for which the computation should be started
 * @param cat_list List of categories for the computation
 * @param p
 */
void run(GeneSet& test_set, CategoryList& cat_list, const Params& p);

/**
 * This function runs the entire pipeline.
 *
 * @param test_set Test set for which the computation should be started
 * @param cat_list List of categories for the computation
 * @param p
 * @param Flag indicating if the p-value should be computed separately.
 */
void run(GeneSet& test_set, CategoryList& cat_list, const Params& p, bool computePValue);

#endif // GT2_APPLICATIONS_ENRICHMENT_COMMON_H

