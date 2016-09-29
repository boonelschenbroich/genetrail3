#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/EntityDatabase.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/OverRepresentationAnalysis.h>
#include <genetrail2/regulation/RegulatorCategoryFileReader.h>

#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <map>

using namespace GeneTrail;

namespace bpo = boost::program_options;

std::string input_, regulations_, output_;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")(
	    "testset,t", bpo::value(&input_)->required(),
	    "A list of deregulated targets.")(
	    "regulations,r", bpo::value(&regulations_)->required(),
	    "A whitespace separated file containing regulator, target and scores.")(
	    "output,o", bpo::value(&output_)->required(),
	    "Output prefix for text files.");

	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),
		           vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	return true;
}

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv)) {
		return -1;
	}

	EntityDatabase db;

	RegulatorCategoryFileReader parser(regulations_, &db);
	parser.parse();

	Category reference = parser.getReference();
	std::map<std::string, Category> categories = parser.getCategories();

	GeneSetReader reader;
	GeneSet test_set = reader.readGeneList(input_);

	OverRepresentationAnalysis ora(reference, test_set.toCategory("test", &db));

	std::ofstream out;
	out.open(output_);
	if(out.is_open()) {
		for(auto it = categories.begin(); it != categories.end(); ++it) {
			std::cout << "INFO: Processing: " << it->first << std::endl;
			
			out << it->first << "\t" << ora.numberOfHits(it->second)
			    << "\t" << ora.expectedNumberOfHits(it->second) << "\t"
			    << ora.computeUpperTailedPValue(it->second) << "\n";
		}
	}
	out.close();

	return 0;
}
