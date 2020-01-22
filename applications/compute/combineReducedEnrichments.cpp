#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

#include <genetrail2/core/Exception.h>
#include <genetrail2/core/CombineReducedEnrichments.h>


using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string samples = "", output = "";

bool parseArguments(int argc, char* argv[]){
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("samples,s", bpo::value<std::string>(&samples)->required(), "Path to a file listing the samples along with their output directories of the enrichment analyses that should be combined")
		("out_files,o", bpo::value<std::string>(&output)->required(), "Path to a file containing the measured category databases along with a path to an file for the created matrix (for each category)");

	try{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);
		bpo::notify(vm);
	} catch(bpo::error& e){
		std::cerr << "ERROR: " << e.what() << std::endl;
		desc.print(std::cerr);
		return false;
	}
	return true;
}

int main(int argc, char* argv[]){
	if(!parseArguments(argc, argv)) return -1;

	try{
		CombineReducedEnrichments c;
		c.writeFiles(samples, output);
	} catch(const IOError& e) {
		std::cerr << "ERROR: " << e.what() << std::endl;
		return -1;
	}
	return 0;
}
