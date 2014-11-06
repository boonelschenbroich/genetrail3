#include <genetrail2/core/Category.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/OverRepresentationAnalysis.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/ORAResult.h>

#include "common.h"

#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <iostream>
#include <fstream>
#include <cstdint>
#include <utility>
#include <tuple>
#include <map>
#include <functional>

using namespace GeneTrail;
namespace bpo = boost::program_options;
namespace bm = boost::math;
using namespace boost::multiprecision;

std::string reference;
Params p;

GeneSet test_set;
GeneSet reference_set;
CategoryList cat_list;
OverRepresentationAnalysis ora;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	addCommonCLIArgs(desc, p);
	desc.add_options()
		("identifier, d", bpo::value<std::string>(&p.identifier), "A file containing identifier line by line.")
		("reference, r", bpo::value<std::string>(&reference)->required(), "A file containing identifier line by line.");

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

	return true;
}

std::shared_ptr<EnrichmentResult> computeEnrichment(const Category& c, const std::pair<int, std::string>& genes)
{
	auto result = std::make_shared<ORAResult>();
	*result = ora.computePValue(c);
	return result;
}

void computePValues(AllResults& results)
{
}

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv)) {
		return -1;
	}

	if(init(test_set,cat_list,p) != 0)
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

	ora = OverRepresentationAnalysis(reference_set.toCategory("reference"), test_set.toCategory("test"));

	run(test_set, cat_list, p);

	return 0;
}

