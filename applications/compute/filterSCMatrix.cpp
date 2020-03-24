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
#include <genetrail2/core/SCMatrixFilter.h>

#include "../matrixTools.h"

using namespace GeneTrail;
namespace bpo = boost::program_options;

FilterParams params;
std::string matrix = "";
MatrixReaderOptions options;

bool parseArguments(int argc, char* argv[]){
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("matrix,m", bpo::value<std::string>(&matrix)->required(), "Name of a text file containing unfiltered expression values as a matrix.")
		("min-total-count", bpo::value<double>(&params.min_total_count)->default_value(0), "The minimal number of counts allowed for a cell to pass.")
		("max-total-count", bpo::value<double>(&params.max_total_count)->default_value(std::numeric_limits<double>::max()), "The maximal number of counts allowed for a cell to pass.")
		("min-features", bpo::value<double>(&params.min_features)->default_value(0), "The minimal number of non-zero genes allowed for a cell to pass.")
		("max-features", bpo::value<double>(&params.max_features)->default_value(std::numeric_limits<double>::max()), "The maximal number of non-zero genes allowed for a cell to pass.")
		("max-mitochondrial", bpo::value<double>(&params.max_mito)->default_value(1), "The maximal percentage of mitochondrial counts allowed for a cell to pass.")
		("mitochondrial-genes", bpo::value<std::string>(&params.mito_genes)->required(), "A file containing mitochondrial genes as symbols.")
		("statistics-file,e", bpo::value<std::string>(&params.out_statistics)->required(), "Name of the resulting statistics file.")
		("output,o", bpo::value<std::string>(&params.out_matrix)->required(), "Name of the filtered output file.");

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

int main(int argc, char* argv[]){
	if(!parseArguments(argc, argv)) return -1;

	try{
		std::cout << "Reading matrix..." << std::endl;
		options.split_only_tab = true;
		auto m = readDenseMatrix(matrix, options);
		
		SCMatrixFilter filter;
		filter.filterMatrix(m, params);
	} catch (std::invalid_argument e){
		std::cout << e.what() << std::endl;
		return -1;
	}

	return 0;
}
