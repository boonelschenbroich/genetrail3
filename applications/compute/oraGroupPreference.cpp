#include <iostream>
#include <fstream>
#include <map>

#include <boost/program_options.hpp>

#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrixWriter.h>
#include <genetrail2/core/MetadataReader.h>
#include <genetrail2/core/ORAGroupPreference.h>

#include "../matrixTools.h"

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string matrix = "", metadata = "", column = "", output = "";
double threshold = 0.05;
MatrixReaderOptions options;

bool parseArguments(int argc, char* argv[]){
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("matrix,m", bpo::value<std::string>(&matrix)->required(), "Name of a (category x sample) matrix file containing p-values for the category in the respective sample.")
		("metadata,e", bpo::value<std::string>(&metadata)->required(), "Name of a tab-separated metadata file in which the first column is a sample and the following columns are metadata information about that sample. One column has to store the name of the group to which the sample belongs. This file needs to have a header that has one element less than the following rows.")
		("threshold,t", bpo::value<double>(&threshold)->required(), "The p-value threshold to name a p-value 'significant'.")
		("column,c", bpo::value<std::string>(&column)->required(), "Name of the column in the metadata file that stores group information.")
		("output,o", bpo::value<std::string>(&output)->required(), "An output file for the (category x group) matrix storing a p-value on how significant a category is only present in the respective group.");

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

void writeMatrix(DenseMatrix& m){
	DenseMatrixWriter writer;
	std::filebuf fb;
	
	fb.open(output,std::ios::out);
	std::ostream os(&fb);

	writer.writeText(os, m);
	fb.close();
}

int main(int argc, char* argv[]){
	if(!parseArguments(argc, argv)) return -1;

	try{
		MetadataReader reader;
		std::ifstream strm(metadata);
		auto meta = reader.readMetadataFile(strm, column);
		options.split_only_tab = true;
		auto m = readDenseMatrix(matrix, options);
		
		ORAGroupPreference ogp(threshold);
		DenseMatrix result(0,0);
		ogp.calculatePreference(m, meta, result);
		
		writeMatrix(result);
	} catch(const IOError& e) {
		std::cerr << "ERROR: " << e.what() << std::endl;
		return -1;
	}
	return 0;
}
