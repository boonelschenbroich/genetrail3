#include <boost/program_options.hpp>

#include <genetrail2/core/Exception.h>
#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrix.h>

#include <genetrail2/regulation/RegulatorGeneAssociationEnrichmentAnalysis.h>
#include <genetrail2/regulation/RegulatorGeneAssociationEnrichmentAlgorithms.h>
#include <genetrail2/regulation/RegulationFileParser.h>
#include <genetrail2/regulation/RegulatorImpactScore.h>

#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_set>

#include "../matrixTools.h"

using namespace GeneTrail;

namespace bpo = boost::program_options;

std::string scores_, matrix_, regulations_, out_, method_, adjustment_method_, impact_score_, confidence_interval_;
bool useAbsoluteValues_, decreasingly_, perturbation_, json_;
size_t seed_, bootstrapping_runs_;
double alpha_;

MatrixReaderOptions matrixOptions;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
	("scores,s", bpo::value(&scores_)->required(), "A whitespace separated file containing deregulated targets.")
	("matrix,x", bpo::value(&matrix_)->required(), "A whitespace separated file containing expression values for all genes.")
	("no-row-names,w", bpo::value<bool>(&matrixOptions.no_rownames)->default_value(false)->zero_tokens(), "Does the file contain row names.")
	("no-col-names,c", bpo::value<bool>(&matrixOptions.no_colnames)->default_value(false)->zero_tokens(), "Does the file contain column names.")
	("add-col-name,n", bpo::value<bool>(&matrixOptions.additional_colname)->default_value(false)->zero_tokens(), "File containing two lines specifying which rownames belong to which group.")
	("decreasingly,d", bpo::value(&decreasingly_)->default_value(false)->zero_tokens(), "Sort testset decreasingly (default: increasingly).")
	("regulations,r", bpo::value(&regulations_)->required(), "A whitespace separated file containing regulator and target.")
	("abs,a", bpo::value(&useAbsoluteValues_)->default_value(false)->zero_tokens(), "Use absolute correlations.")
	("method,m", bpo::value(&method_)->required(), "The method that should be applied (ks-test, wrs-test).")
	("output,o", bpo::value(&out_)->required(), "Output prefix for text files.")
	("seed,e", bpo::value(&seed_)->required(), "Random seed used for pertubation.")
	("bootstrap,b", bpo::value(&bootstrapping_runs_)->required(), "Number of bootstrapping runs.")
	("alpha,l", bpo::value(&alpha_)->required()->default_value(0.1), "Alpha level of confidence interval.")
	("adjust,u", bpo::value(&adjustment_method_)->required()->default_value("benjamini-yekutieli"), "Alpha level of confidence interval.")
	("json,j", bpo::value(&json_)->default_value(false)->zero_tokens(), "Output file in .json format (default: .tsv).")
	("impact,p", bpo::value(&impact_score_)->required(), "Method that should be used to compute the impact of regulators.")
	("confidence-intervals, v", bpo::value(&confidence_interval_)->required(), "Method that should be used to compute confidence intervals.")
	;

	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),
		           vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	return true;
}

std::vector<size_t>
translate_test_set(const DenseMatrix* matrix,
                   const std::vector<std::string>& test_set_names)
{
	std::vector<size_t> test_set;
	test_set.reserve(test_set_names.size());
	for(const std::string name : test_set_names) {
		test_set.emplace_back(matrix->rowIndex(name));
	}
	return test_set;
}

template <typename RegulatorImpactScore>
void run(RegulatorEnrichmentAnalysis<double>& analysis, RegulatorImpactScore impactScore)
{
	if(method_ == "ks-test") {
		std::cout << "INFO: Performing KS-test" << std::endl;
		analysis.run(KSTest(), impactScore);
	} else if(method_ == "wrs-test") {
		std::cout << "INFO: Performing WRS-test" << std::endl;
		analysis.run(WRSTest(), impactScore);
	} else if(method_ == "em-test") {
		std::cout << "INFO: Performing EM-test" << std::endl;
		analysis.run(EMTest(), impactScore);
	} else {
		std::cerr << "ERROR: Method '" << method_ << "' not known!" << std::endl;  
	}
}

void run(RegulatorEnrichmentAnalysis<double>& analysis)
{
	if(impact_score_ == "pearson_correlation") {
		run(analysis, PearsonCorrelation());
	} else if(impact_score_ == "spearman_correlation") {
		run(analysis, SpearmanCorrelation());
	} else if(impact_score_ == "kendall_correlation") {
		run(analysis, KendallCorrelation());
	} else {
		std::cerr << "ERROR: Impact score '" << impact_score_ << "' not known!" << std::endl;  
	}
}

void write(RegulatorEnrichmentAnalysis<double>& analysis){
	if(json_) {
		analysis.writeJSONResults(out_, alpha_, confidence_interval_);
	} else {
		analysis.writeResults(out_, alpha_, confidence_interval_);
	}
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
	std::vector<std::string> sorted_target_names =
	    test_set.getSortedIdentifier(decreasingly_);
	std::vector<size_t> sorted_targets =
	    translate_test_set(&matrix, sorted_target_names);

	std::cout << "INFO: Parsing regulations" << std::endl;
	std::unordered_set<size_t> tset(sorted_targets.begin(),
	                                sorted_targets.end());
	RegulationFileParser<double> parser(&matrix, tset, regulations_, 0.0);

	RegulatorEnrichmentAnalysis<double> analysis(
	    sorted_targets, parser, &matrix, seed_, useAbsoluteValues_,
	    bootstrapping_runs_);
	run(analysis);

	std::cout << "INFO: Adjusting p-values" << std::endl;
	analysis.adjustPValues(adjustment_method_);

	std::cout << "INFO: Writing results" << std::endl;
	write(analysis);
		
	return 0;
}
