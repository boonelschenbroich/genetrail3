#ifndef ENRICH_COMMON_H
#define ENRICH_COMMON_H

#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

#include <genetrail2/core/Exception.h>
#include <genetrail2/core/EnrichmentResult.h>
#include <genetrail2/core/Category.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/PValue.h>

#include <list>
#include <utility>
#include <fstream>
#include <string>
#include <iostream>

using namespace GeneTrail;

namespace bpo = boost::program_options;

typedef std::map<std::string, EnrichmentResult> Results;
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
};

void addCommonCLIArgs(bpo::options_description& desc, Params& p);

int readTestSet(GeneSet<double>& test_set, const Params& p);

typedef std::list<std::pair<std::string, std::string>> CategoryList;
CategoryList getCategoryList(const std::string& catfile_list, const std::string& single_cat);

PValueList resultVector(const AllResults& results);

std::pair<bool, std::pair<int, std::string>> processCategory(Category& c, GeneSet<double>& test_set, const Params& p);

std::map<std::string, std::vector<EnrichmentResult>> splitDatabases(AllResults& all_results, const PValueList& pvalues);

void writeFile(const std::string& output_dir, const std::map<std::string,std::vector<EnrichmentResult>>& databases);

int init(GeneSet<double>& test_set, CategoryList& cat_list, const Params& p);

EnrichmentResult computeEnrichment(const Category& c, std::pair<int, std::string> genes);

void run(GeneSet<double>& test_set, CategoryList& cat_list, AllResults& name_to_cat_results, const Params& p);

#endif // ENRICH_COMMON_H

