#include <iostream>
#include <cmath>
#include <cstdlib>
#include <tuple>

#include <boost/program_options.hpp>

#include <genetrail2/core/MatrixHTest.h>
#include <genetrail2/core/GeneSetWriter.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/TextFile.h>

#include "../matrixTools.h"

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string expr1 = "", expr2 = "", output = "", method = "", groups = "";
bool binary = false;

MatrixReaderOptions matrixOptions;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("expression-matrix-1,1", bpo::value<std::string>(&expr1)->required(), "Name of a text file containing expression values as a matrix.")
		("expression-matrix-2,2", bpo::value<std::string>(&expr2), "Name of a text file containing expression values as a matrix.(optional)")
		("output,o", bpo::value<std::string>(&output)->required(), "Name of the output file.")
		("groups,g", bpo::value<std::string>(&groups)->required(), "File containing two lines specifying which rownames belong to which group.")
		("no-row-names,r", bpo::value<bool>(&matrixOptions.no_rownames)->default_value(false)->zero_tokens(), "Does the file contain row names.")
		("no-col-names,c", bpo::value<bool>(&matrixOptions.no_colnames)->default_value(false)->zero_tokens(), "Does the file contain column names.")
		("add-col-name,a", bpo::value<bool>(&matrixOptions.additional_colname)->default_value(false)->zero_tokens(), "File containing two lines specifying which rownames belong to which group.")
		("method,m", bpo::value<std::string>(&method)->required(), "Method used for scoring.");

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

	TextFile t(groups, ",", std::set<std::string>());
	DenseMatrix matrix(0,0);

	try {
		matrix = buildDenseMatrix(expr1, expr2, matrixOptions);
	} catch(const IOError& e) {
		std::cerr << "ERROR: Could not open input data matrix for reading." << std::endl;
		return -4;
	}

	std::vector<std::string> reference, sample;
	try {
		reference = t.read();
		sample = t.read();
	} catch(const IOError& e) {
		std::cerr << "ERROR: Could not read from group file " << groups << std::endl;
		return -5;
	}

	try {
		auto subset = splitMatrix(matrix, reference, sample);

		MatrixHTest htest;
		auto gene_set = htest.test(method, std::get<0>(subset), std::get<1>(subset));

		for(const auto& score : gene_set.scores()) {
			if(isnan(score)) {
				std::cerr
				    << "WARNING: NaNs generated during score computation.\n";
				break;
			}
		}

		GeneSetWriter writer;
		writer.writeScoringFile(gene_set, output);
	} catch(const EmptyGroup& e) {
		std::cerr << "ERROR: " << e.what() << "\n";
		return -3;
	} catch(const std::invalid_argument& e) {
		std::cerr << "ERROR: Unknown method '" << e.what() << "'\n";
		return -6;
	}

	return 0;
}
