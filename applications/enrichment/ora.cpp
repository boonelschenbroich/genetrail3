#include <genetrail2/core/Category.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/OverRepresentationAnalysis.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/ScoringFile.h>
#include <genetrail2/core/TestSet.h>

#include "common.h"

#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <iostream>
#include <fstream>
#include <cstdint>
#include <utility>
#include <tuple>
#include <map>

using namespace GeneTrail;
namespace bpo = boost::program_options;
namespace bm = boost::math;
using namespace boost::multiprecision;

std::string reference;
Params p;

struct OraResult
{
	OraResult() : expected_hits(0.0), hits(0)
	{
	}

	std::string name;
	double pvalue;
	std::string reference;
	double expected_hits;
	unsigned int hits;
	std::string info;
};

typedef std::map<std::string, OraResult> Results;
typedef std::map<std::string, Results> AllResults;
typedef std::vector<std::pair<std::string, double>> PValueList;

PValueList resultVector(const AllResults& results)
{
	PValueList result;

	for(const auto& it : results) {
		for(const auto& jt : it.second) {
			result.push_back(std::make_pair(jt.second.name, jt.second.pvalue));
		}
	}

	return result;
}

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	addCommonCLIArgs(desc, p);
	desc.add_options()
		("reference, r", bpo::value<std::string>(&reference)->required(), "A file containing identifier line by line.");

	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),
		           vm);
		bpo::notify(vm);
	}
	catch(bpo::error& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	if(p.scores != "" && p.identifier != "") {
		std::cerr << "ERROR: Please specify only one input file." << std::endl;
		return false;
	} else if(p.scores == "" && p.identifier == "") {
		std::cerr << "ERROR: Please specify a input file." << std::endl;
		return false;
	}

	return true;
}

bool processCategory(const Category& c, const TestSet& test_set)
{
	int hits = 0;
	for(auto s : test_set) {
		if(c.contains(s.first)) {
			++hits;
		}
	}

	return p.minimum <= hits && hits <= p.maximum;
}

OraResult computeEnrichment(const Category& c, const OverRepresentationAnalysis& ora)
{
	std::cout << "INFO: Processing " << c.name() << std::endl;

	OraResult result;

	result.name = c.name();
	result.reference = c.reference();

	auto enr = ora.computePValue(c);

	//TODO: Currently there is no fourth argument
	result.pvalue        = std::get<0>(enr);
	result.expected_hits = std::get<1>(enr);
	result.info          = std::get<2>(enr);
	result.hits          = std::get<3>(enr);

	return result;
}

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv)) {
		return -1;
	}

	TestSet test_set;
	try
	{
		std::ifstream bla(p.scores);
		test_set = TestSet::readSet(bla);
	}
	catch(IOError& exn)
	{
		std::cerr << "ERROR: failed to read test set. Reason: " << exn.what()
		          << std::endl;
		return -1;
	}

	GeneSetReader reader;
	auto reference_set = reader.readGeneSet(reference, "reference");

	CategoryList cat_list;
	try
	{
		//TODO: Add single category feature
		cat_list = getCategoryList(p.categories, std::string());
	}
	catch(IOError& exn)
	{
		std::cerr << "ERROR: failed to read categories. Reason: " << exn.what()
		          << std::endl;
		return -1;
	}

	AllResults name_to_cat_results;

	OverRepresentationAnalysis ora(reference_set, test_set);
	for(const auto& cat : cat_list) {
		try
		{
			GMTFile input(cat.second);

			Results name_to_result;
			while(input) {
				Category c = input.read();

				if(!processCategory(c, test_set)) {
					continue;
				}

				name_to_result.insert(std::make_pair(c.name(), computeEnrichment(c, ora)));
			}

			name_to_cat_results.insert(std::make_pair(cat.first, name_to_result));
		}
		catch(IOError& exn)
		{
			std::cerr << "WARNING: could not process category file " << cat.first
			          << " skipping! " << std::endl;
		}

	}

	auto results = resultVector(name_to_cat_results);
	results = pvalue<double>::adjustPValues(results, p.adjustment);

	return 0;
}

