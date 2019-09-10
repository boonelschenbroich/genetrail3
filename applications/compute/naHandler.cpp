#include <iostream>
#include <fstream>
#include <functional> 
#include <map>
#include <string> 

#include <boost/program_options.hpp>

#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrixWriter.h>
#include <genetrail2/core/NAHandler.h>

#include "../matrixTools.h"

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string input = "", output = "", strategy = "";
MatrixReaderOptions options;

bool parseArguments(int argc, char* argv[]){
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("input,i", bpo::value<std::string>(&input)->required(), "Name of the input file with NAs.")
		("strategy,s", bpo::value<std::string>(&strategy)->required(), "The strategy for handling NAs.")
		("output,o", bpo::value<std::string>(&output)->required(), "Name of the output file.");

	try{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);
		bpo::notify(vm);
	} catch(bpo::error& e){
		std::cerr << "ERROR: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}
	return true;
}

bool handleNAs(DenseMatrix& m){
	NAHandler h;
	if(strategy.compare("remove") == 0){
		h.handle(m, NAStrategyRemove());
	} else if(strategy.compare("mean") == 0){
		h.handle(m, NAStrategyMean());
	} else if(strategy.compare("median") == 0){
		h.handle(m, NAStrategyMedian());
	} else if(strategy.compare("constantzero") == 0){
		h.handle(m, NAStrategyZero());
	} else {
		std::cerr << "ERROR: Strategy is unknown." << std::endl;
		return false;
	}
	return true;
}

void writeTransformedMatrix(DenseMatrix& m){
	DenseMatrixWriter writer;
	std::filebuf fb2;
	
	fb2.open (output,std::ios::out);
	std::ostream os(&fb2);

	writer.writeText(os, m);
	fb2.close();
}


int main(int argc, char* argv[]){ 
  	if(!parseArguments(argc, argv)) return -1;
	DenseMatrix m(0,0);
	try {m = readDenseMatrix(input, options);}
	catch(const IOError& e) {
		std::cerr << "ERROR: Could not read from matrix file " << input << std::endl;
		return -1;
	}
	if(!handleNAs(m)) return -1;
	try {writeTransformedMatrix(m);}
	catch(const IOError& e) {
	  	std::cerr << "ERROR: Could not write transformed matrix in file " << output << std::endl;
		return -1;
	}
	return 0;
}
