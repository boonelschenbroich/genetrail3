#include <iostream>
#include <fstream>
#include <tuple>

#include <boost/program_options.hpp>

#include <genetrail2/core/MatrixNormalization.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrixWriter.h>
#include <genetrail2/core/TextFile.h>
#include <genetrail2/core/DenseColumnSubset.h>

#include "../matrixTools.h"


using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string matrix = "", output = "", normalization = "", samples = "";
DenseMatrix valueMatrix(0,0);
MatrixReaderOptions matrixOptions;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("matrix,m", bpo::value<std::string>(&matrix)->required(), "Name of the matrix file.")
		("output,o", bpo::value<std::string>(&output)->required(), "Name of the output file.")
		("normalization,n", bpo::value<std::string>(&normalization)->required(), "Method to normalize the values (gauss, poisson, zscore).")
		("no-row-names,r", bpo::value<bool>(&matrixOptions.no_rownames)->default_value(false)->zero_tokens(), "Does the file contain row names.")
		("no-col-names,c", bpo::value<bool>(&matrixOptions.no_colnames)->default_value(false)->zero_tokens(), "Does the file contain column names.")
		("add-col-name,a", bpo::value<bool>(&matrixOptions.additional_colname)->default_value(false)->zero_tokens(), "Additional column names")
		("samples,s", bpo::value<std::string>(&samples)->required(), "File containing one line specifying which rownames should be used in the analysis.");


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

void readValueMatrix()
{
	valueMatrix = readDenseMatrix(matrix, matrixOptions);
}


bool writeNormalizedMatrix(DenseColumnSubset subset)
{
	MatrixNormalization norm;
	DenseMatrixWriter writer;
	std::filebuf fb2;
	
	fb2.open (output,std::ios::out);
	std::ostream os(&fb2);
	
	if(normalization.compare("none") == 0) {
	  writer.writeText(os, subset);
	  fb2.close();
	  return true;
	}
	
	else if(normalization.compare("gauss")== 0) {
	  GaussEstimator gauss;
	  writer.writeText(os,norm.normalizeMatrix(subset,gauss));
	  fb2.close();
	  return true;
	}
	else if(normalization.compare("poisson") == 0) {
	  PoissonEstimator poisson;
	  writer.writeText(os,norm.normalizeMatrix(subset,poisson));
	  fb2.close();
	  return true;
	}
	else if(normalization.compare("zscore") == 0) {
	  
	  if(valueMatrix.cols() < 3) {
	    std::cerr << "WARNING: Z-Score for < 3 samples will always return 0 as standard deviation is 0." << std::endl;
	  }

	  ExclusiveZScore zscore;
	  writer.writeText(os,norm.normalizeMatrix(subset,zscore));
	  fb2.close();
	  return true;
	}
	else {
	  fb2.close();
	  return false;
	}
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
	
	//Decide which samples will be used for the normalization
	std::vector<std::string> sample, reference;
	TextFile t(samples, ",", std::set<std::string>());
	
	try {
		sample = t.read();
	} catch(const IOError& e) {
		std::cerr << "ERROR: Could not read from samples file " << samples << std::endl;
		return -4;
	}

	
	try {
	  
	   std::tuple<DenseColumnSubset, DenseColumnSubset> subset = splitMatrix(valueMatrix, reference, sample);
	   DenseColumnSubset inputSamples = std::get<1>(subset);
	   
	   if (!writeNormalizedMatrix(inputSamples)) {
	      std::cerr << "ERROR: Could not parse valid normalization-method" << "\n";
	      return -5;
	   }
	} catch(const EmptyGroup& e) { 
		std::cerr << "ERROR: " << e.what() << "\n";
		return -6;
	} catch(const IOError& e) {
	  	std::cerr << "ERROR: Could not write normalized matrix in file " << output << std::endl;
		return -7;
	}
	
	

  return 0;
  
}
