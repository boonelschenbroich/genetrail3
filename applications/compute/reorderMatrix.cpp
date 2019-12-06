#include <iostream>
#include <fstream>
#include <functional> 
#include <map>
#include <string> 

#include <boost/program_options.hpp>

#include <genetrail2/core/TextFile.h>
#include <genetrail2/core/Matrix.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrixWriter.h>

#include "../matrixTools.h"

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string matrix_name = "", output = "", sample_names = "";

MatrixReaderOptions matrixOptions;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("matrix,m", bpo::value<std::string>(&matrix_name)->required(), "Name of the matrix file.")
		("output,o", bpo::value<std::string>(&output)->required(), "Name of the output matrix file.")
		("sample-names,s", bpo::value<std::string>(&sample_names)->required(), "Samples ...")
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

	DenseMatrix matrix(0,0);
	try {
		matrix = readDenseMatrix(matrix_name, matrixOptions);
	} catch(const IOError& e) {
		std::cerr << "ERROR: Could not read from matrix file " << matrix_name << std::endl;
		return -3;
	}

	TextFile t(sample_names, ",", std::set<std::string>());
	std::vector<std::string> samples;
	try {
		samples = t.read();
	} catch(const IOError& e) {
		std::cerr << "ERROR: Could not read from group file " << sample_names << std::endl;
		return -5;
	}

	DenseColumnSubset dcsubset = DenseColumnSubset(&matrix, getIndices(matrix, samples, "reordered"));
	DenseMatrixWriter writer;
	std::filebuf fb2;
	
	fb2.open (output, std::ios::out);
	std::ostream os(&fb2);

	writer.writeText(os, dcsubset);
	fb2.close();

  	return 0;
  
}
