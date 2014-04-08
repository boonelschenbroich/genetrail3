#include "GeneSetReader.h"

using namespace GeneTrail;

Category GeneSetReader::readGeneSet(const std::string& path, const std::string name){
	Category c(name);
	std::ifstream input(path);
	if(!input) {
		std::cerr << "ERROR: Could not open " << path << " for reading." << std::endl;
	}
	std::string gene;
	for(std::string gene; getline(input, gene);) {
		boost::trim(gene);
		c.insert(gene);
	}
	return c;
}

std::vector<std::string> GeneSetReader::readGeneList(const std::string& path){
	std::vector<std::string> list;
	std::ifstream input(path);
	if(!input) {
		std::cerr << "ERROR: Could not open " << path << " for reading." << std::endl;
	}
	std::string gene;
	for(std::string gene; getline(input, gene);) {
		boost::trim(gene);
		list.push_back(gene);
	}
	return list;
}
