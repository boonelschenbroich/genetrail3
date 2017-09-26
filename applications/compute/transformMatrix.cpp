#include <iostream>
#include <fstream>
#include <functional> 
#include <map>
#include <string> 

#include <boost/program_options.hpp>

#include <genetrail2/core/MatrixTransformation.h>
#include <genetrail2/core/Matrix.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrixWriter.h>

#include "../matrixTools.h"

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string matrix = "", output = "", transformation = "";

bool use_ranks = false;
DenseMatrix valueMatrix(0,0);
MatrixReaderOptions matrixOptions;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("matrix,m", bpo::value<std::string>(&matrix)->required(), "Name of the matrix file.")
		("output,o", bpo::value<std::string>(&output)->required(), "Name of the output file.")
		("ranks,s", bpo::value<bool>(&use_ranks)->default_value(false)->zero_tokens(), "Transform matrix values to ranks for each column.")
		("transformation,t", bpo::value<std::string>(&transformation)->default_value("default"), "Method to transform the values.")
		("no-row-names,r", bpo::value<bool>(&matrixOptions.no_rownames)->default_value(false)->zero_tokens(), "Does the file contain row names.")
		("no-col-names,c", bpo::value<bool>(&matrixOptions.no_colnames)->default_value(false)->zero_tokens(), "Does the file contain column names.")
		("add-col-name,a", bpo::value<bool>(&matrixOptions.additional_colname)->default_value(false)->zero_tokens(), "File containing two lines specifying which rownames belong to which group.");

	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);
		bpo::notify(vm);
	}
	catch(bpo::error& e)
	{
		std::cerr << "ERROR: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	return true;
}

bool applyTransformation()
{
  
  	if(use_ranks) 
	  valuesToRanks(valueMatrix);
	
	//no transformation
	if(transformation.compare("default")== 0)
	  return true;
	
	else if(transformation.compare("upweightEnds")== 0) {
	  if(!use_ranks) {
	    std::cerr << "ERROR: upweightEnds should only be used on ranks" << std::endl;
	    return false;
	  }
	  upweightEnds(valueMatrix);
	}
	
	else if(transformation.compare("upweightTail")== 0) {
	  if(!use_ranks) {
	    std::cerr << "ERROR: upweightTail should only be used on ranks" << std::endl;
	    return false;
	  }
	  upweightTail(valueMatrix);
	}
	
	else if(transformation.compare("abs")== 0){
	  abs(valueMatrix);
	}else if(transformation.compare("sqrt")== 0){
	  sqrt(valueMatrix);
	}else if(transformation.compare("log")== 0){
	  log(valueMatrix);
	}else if(transformation.compare("log2")== 0){
	  log2(valueMatrix);
	}else if(transformation.compare("log10")== 0){
	  log10(valueMatrix);
	}else if(transformation.compare("pow2")== 0){
	  pow2(valueMatrix);
	}else {
	  std::cerr << "ERROR: Transformation is unknown." << std::endl;
	  return false;
	}
	
	return true;
}

void readValueMatrix()
{
	valueMatrix = readDenseMatrix(matrix,matrixOptions);
}


void writeTransformedMatrix()
{
	DenseMatrixWriter writer;
	std::filebuf fb2;
	
	fb2.open (output,std::ios::out);
	std::ostream os(&fb2);

	writer.writeText(os,valueMatrix);
	fb2.close();

}


int main(int argc, char* argv[])
{ 
  	if(!parseArguments(argc, argv))
	{
		return -1;
	}
	
	if(matrixOptions.additional_colname && matrixOptions.no_rownames) {
		std::cerr << "Conflicting arguments. Additional colnames can only be "
		             "specified if row names are present!" << std::endl;
		return -2;
	}

	
	
	try {readValueMatrix();}
	catch(const IOError& e) {
		std::cerr << "ERROR: Could not read from matrix file " << matrix << std::endl;
		return -3;
	}
	
	if(!applyTransformation())
	{
		return -4;
	}
	
	try {writeTransformedMatrix();}
	catch(const IOError& e) {
	  	std::cerr << "ERROR: Could not write transformed matrix in file " << output << std::endl;
		return -5;
	}
	
	

  return 0;
  
}
