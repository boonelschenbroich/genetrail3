#include <boost/program_options.hpp>

#include <genetrail2/core/Exception.h>
#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseColumnSubset.h>
#include <genetrail2/core/TextFile.h>

#include <genetrail2/regulation/RegulationFileParser.h>
#include <genetrail2/regulation/RegulationFile.h>
#include <genetrail2/regulation/RegulatorEffectResult.h>
#include <genetrail2/regulation/CorrelationSetAnalysis.h>
#include <genetrail2/regulation/RegulatorAssociationScore.h>

#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility> 
#include <map> 

#include "common.h"
#include "../matrixTools.h"

using namespace GeneTrail;

namespace bpo = boost::program_options;

std::string scores_ = "", matrix_ = "", regulations_ = "", out_ = "", method_ = "", adjustment_method_ = "";
size_t seed_, permutations_;
bool upper_tailed_, abs_;

MatrixReaderOptions matrixOptions;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
					  ("scores,s", bpo::value(&scores_)->required(), "A whitespace separated file containing genes of interest.")
					  ("matrix,x", bpo::value(&matrix_)->default_value(""), "A whitespace separated file containing expression values for all genes.")
					  ("no-row-names,w", bpo::value<bool>(&matrixOptions.no_rownames)->default_value(false)->zero_tokens(), "Does the file contain row names.")
					  ("no-col-names,c", bpo::value<bool>(&matrixOptions.no_colnames)->default_value(false)->zero_tokens(), "Does the file contain column names.")
					  ("add-col-name,n", bpo::value<bool>(&matrixOptions.additional_colname)->default_value(false)->zero_tokens(), "File containing two lines specifying which rownames belong to which group.")
					  ("regulations,r", bpo::value(&regulations_)->default_value(""), "A whitespace separated file containing regulator and target.")
					  ("seed,e", bpo::value(&seed_)->default_value(0), "Random seed used for pertubation.")
					  ("permutations,p", bpo::value(&permutations_)->default_value(0), "Number of bootstrapping runs.")
					  ("method,m", bpo::value(&method_)->required(), "The method that should be applied (pearson-correlation, spearman-correlation, kendall-correlation).")
					  ("adjust,a", bpo::value(&adjustment_method_)->required(), "Method to adjust p-values.")
					  ("output,o", bpo::value(&out_)->required(), "Output prefix for text files.")
					  ("upper-tailed,u", bpo::value(&upper_tailed_)->default_value(false)->zero_tokens(), "Calculate an upper-tailed p-value.")
					  ("abs,b", bpo::value(&upper_tailed_)->default_value(false)->zero_tokens(), "Use absolute value.");
	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "ERROR: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	return true;
}

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv)) {
		return -1;
	}

	if(matrixOptions.additional_colname && matrixOptions.no_rownames) {
		std::cerr << "Conflicting arguments. Additional colnames can only be "
		             "specified if row names are present!" << std::endl;
		return -2;
	}

	std::cout << "INFO: Parsing data matrix" << std::endl;
	DenseMatrix matrix(0, 0);
	try {
		matrix = readDenseMatrix(matrix_, matrixOptions);
	} catch(const IOError& e) {
		std::cerr << "ERROR: Could not open input data matrix for reading."
		          << std::endl;
		return -4;
	}
	
	std::cout << "INFO: Parsing test set" << std::endl;
	GeneSetReader reader;
	GeneSet test_set = reader.readScoringFile(scores_);
	std::vector<std::string> sorted_target_names = test_set.getSortedIdentifier(true);
	std::vector<size_t> sorted_targets = translate_test_set(&matrix, sorted_target_names);
		
	std::cout << "INFO: Parsing regulations" << std::endl;
	std::unordered_set<size_t> tset(sorted_targets.begin(), sorted_targets.end());
	MatrixNameDatabase name_database(&matrix);
	RegulationFileParser<MatrixNameDatabase, double> parser(name_database, tset, regulations_, 0.0);
	RegulationFile<double>& regulationFile = parser.getRegulationFile();
	
	std::cout << "INFO: Performing correlation set analysis" << std::endl;
	CorrelationSetAnalysis<double> csa(&matrix, name_database, sorted_targets, regulationFile);
	std::vector<RegulatorEffectResult> results;
	if(method_ == "pearson_correlation") {
		results = csa.run(PearsonCorrelation(), seed_, permutations_, upper_tailed_, abs_);
	} else if(method_ == "spearman_correlation") {
		results = csa.run(SpearmanCorrelation(), seed_, permutations_, upper_tailed_, abs_);
	} else if(method_ == "kendall_correlation") {
		results = csa.run(KendallCorrelation(), seed_, permutations_, upper_tailed_, abs_);
	}
		
	std::cout << "INFO: Adjusting p-values" << std::endl;
	adjustPValues(results, adjustment_method_);

	std::cout << "INFO: Writing results" << std::endl;
	write(results, out_, true);
	
	return 0;
}
