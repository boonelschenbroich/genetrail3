#include <iostream>
#include <cmath>
#include <cstdlib>
#include <tuple>

#include <boost/program_options.hpp>

#include <genetrail2/core/MatrixHTest.h>
#include <genetrail2/core/GeneSetWriter.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrixWriter.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/GroupedScores.h>

#include "../matrixTools.h"

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string matrix = "", output = "", method = "";
std::vector<std::string> sample_groups, reference_groups;

MatrixReaderOptions options;

bool parseArguments(int argc, char* argv[]){
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("matrix,m", bpo::value<std::string>(&matrix)->required(), "Name of a text file containing mean expression values as a matrix.")
		("output,o", bpo::value<std::string>(&output)->required(), "Name of the output file.")
		("sampleGroups,s", bpo::value<std::vector<std::string>>(&sample_groups)->multitoken()->required(), "List of sample groups. The mean expression value for the groups must be accessible in the matrix")
		("referenceGroups,r", bpo::value<std::vector<std::string>>(&reference_groups)->multitoken()->required(), "List of reference groups. The mean expression value for the groups must be accessible in the matrix")
		("method,m", bpo::value<std::string>(&method)->required(), "Method used for scoring.");

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

	std::cout << sample_groups.size() << ", " << reference_groups.size() << std::endl;
	
	try{
		options.split_only_tab = true;
		auto m = readDenseMatrix(matrix, options);
		
		auto result = DenseMatrix(0,0);
		GroupedScores calculator;
		calculator.calculateGroupedScores(m, sample_groups, reference_groups, method, result);
		writeMatrix(result);
	} catch (std::invalid_argument e){
		std::cout << e.what() << std::endl;
		return -1;
	}

	return 0;
}
