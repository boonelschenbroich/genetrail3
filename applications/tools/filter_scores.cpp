#include <genetrail2/core/EntityDatabase.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetFilters.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/GeneSetWriter.h>
#include <genetrail2/core/Scores.h>

#include <boost/program_options.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace GeneTrail;

namespace bpo = boost::program_options;

using GSFPointer = std::unique_ptr<GeneSetFilter::GeneSetFilter>;

class InvalidFilterSyntax : public std::exception
{
	public:
	explicit InvalidFilterSyntax(const std::string& msg) noexcept : message_(msg) {}

	const char* what() const noexcept override { return message_.c_str(); }

	private:
	std::string message_;
};

GSFPointer parseFilter(const std::string& filter)
{
	std::vector<std::string> splitted;
	boost::split(splitted, filter, boost::is_any_of(":"));

	if(splitted.size() != 2) {
		throw InvalidFilterSyntax("The filter '" + filter +
		                          "' could not be parsed.");
	}

	double parameter;
	try {
		parameter = boost::lexical_cast<double>(splitted[1]);
	} catch(boost::bad_lexical_cast& e) {
		throw InvalidFilterSyntax("Could not convert the filter parameter '" +
		                          splitted[1] + "' to a number.");
	}

	if(splitted[0] == "larger") {
		return std::make_unique<GeneSetFilter::LargerFilter>(parameter);
	}

	if(splitted[0] == "smaller") {
		return std::make_unique<GeneSetFilter::SmallerFilter>(parameter);
	}

	if(splitted[0] == "absLarger") {
		return std::make_unique<GeneSetFilter::AbsLargerFilter>(parameter);
	}

	if(splitted[0] == "absSmaller") {
		return std::make_unique<GeneSetFilter::AbsSmallerFilter>(parameter);
	}

	if(splitted[0] == "upperQuantile") {
		return std::make_unique<GeneSetFilter::UpperQuantileFilter>(parameter);
	}

	if(splitted[0] == "lowerQuantile") {
		return std::make_unique<GeneSetFilter::LowerQuantileFilter>(parameter);
	}

	if(splitted[0] == "quantile") {
		return std::make_unique<GeneSetFilter::QuantileFilter>(parameter);
	}

	throw InvalidFilterSyntax("Unknown filter '" + splitted[0] + "'.");
}

std::vector<GSFPointer> parseFilters(const std::vector<std::string>& filters)
{
	std::vector<GSFPointer> result;
	result.reserve(filters.size());

	for(const auto& f : filters) {
		result.push_back(parseFilter(f));
	}

	return result;
}

int main(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	std::string infile, outfile;
	std::vector<std::string> filters;

	desc.add_options()
		("help,", "Display this message")
		("input,i", bpo::value<std::string>(&infile)->required(),
			"A file containing a whitespace separated list of identifiers and "
			"scores.")
		("output,o", bpo::value<std::string>(&outfile)->required(),
			"The filtered scores.")
		("filter,f", bpo::value<std::vector<std::string>>(&filters)->required(),
			"The filters that should be applied.");

	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),
		           vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return -1;
	}

	try {
		auto geneSetFilters = parseFilters(filters);

		GeneSetReader reader;
		auto gene_set = reader.readScoringFile(infile);

		for(const auto& filter : geneSetFilters) {
			gene_set.filter(filter.get());
		}

		auto db = std::make_shared<EntityDatabase>();
		GeneSetWriter writer;
		writer.writeScoringFile(Scores(gene_set, db), outfile);
	} catch(InvalidFilterSyntax& e) {
		std::cerr << e.what() << std::endl;
		return -2;
	}

	return 0;
}

