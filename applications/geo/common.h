#ifndef GT2_APPLICATIONS_COMPUTE_SCORES_COMMON_H
#define GT2_APPLICATIONS_COMPUTE_SCORES_COMMON_H

#include <boost/program_options.hpp>

#include <iostream>

#include <genetrail2/core/GEO.h>
#include <genetrail2/core/GEOGPLParser.h>

namespace bpo = boost::program_options;

struct Params
{
	std::string geo = "";
	std::string geo_dir = "";
	std::string methodToHandleDuplicates = "";
	std::string output_file = "";
	bool gds = false;
};

bool parseArguments(int argc, char* argv[], Params& p);

void annotateAndWriteGEO(GeneTrail::GEOMap& geo, const Params& p);

#endif //GT2_APPLICATIONS_COMPUTE_SCORES_COMMON_H
