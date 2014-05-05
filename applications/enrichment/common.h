#ifndef ENRICH_COMMON_H
#define ENRICH_COMMON_H

#include <boost/program_options.hpp>

#include <list>

namespace bpo = boost::program_options;

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
void getAdjustedPValues();

typedef std::list<std::pair<std::string, std::string>> CategoryList;
CategoryList getCategoryList(const std::string& catfile_list, const std::string& single_cat);

#endif // ENRICH_COMMON_H

