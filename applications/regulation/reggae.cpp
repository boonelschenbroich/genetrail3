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
#include <genetrail2/regulation/RegulationBootstrapper.h>
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

std::string scores_, associations_, matrix_, regulations_, out_, method_, adjustment_method_,
    impact_score_, confidence_interval_;
bool useAbsoluteValues_, decreasingly_, perturbation_, json_;
size_t seed_, bootstrapping_runs_;
double alpha_;

MatrixReaderOptions matrixOptions;

bool checkIfFileExists(bpo::options_description& desc, std::string fname){
	std::ifstream infile(fname);
	if(!infile.good()){
		std::cerr << "ERROR: File '" + fname + "' does not exist. \n";
		desc.print(std::cerr);
		return false;
	}
	return true;
}

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
					  ("scores,s", bpo::value(&scores_)->required(), "A whitespace separated file containing deregulated targets.")
					  ("associations,t", bpo::value(&associations_)->default_value(""), "A whitespace separated file regulator target and association scores.") 
					  ("matrix,x", bpo::value(&matrix_)->default_value(""), "A whitespace separated file containing expression values for all genes.")
					  ("no-row-names,w", bpo::value<bool>(&matrixOptions.no_rownames)->default_value(false)->zero_tokens(), "Does the file contain row names.")
					  ("no-col-names,c", bpo::value<bool>(&matrixOptions.no_colnames)->default_value(false)->zero_tokens(), "Does the file contain column names.")
					  ("add-col-name,n", bpo::value<bool>(&matrixOptions.additional_colname)->default_value(false)->zero_tokens(), "File containing two lines specifying which rownames belong to which group.")
					  ("decreasingly,d", bpo::value(&decreasingly_)->default_value(false)->zero_tokens(), "Sort testset decreasingly (default: increasingly).")
					  ("regulations,r", bpo::value(&regulations_)->default_value(""), "A whitespace separated file containing regulator and target.")
					  ("abs,a", bpo::value(&useAbsoluteValues_)->default_value(false)->zero_tokens(), "Use absolute correlations.")
					  ("method,m", bpo::value(&method_)->required(), "The method that should be applied (ks-test, wrs-test).")
					  ("output,o", bpo::value(&out_)->required(), "Output prefix for text files.")
					  ("seed,e", bpo::value(&seed_)->default_value(0), "Random seed used for pertubation.")
					  ("bootstrap,b", bpo::value(&bootstrapping_runs_)->default_value(0), "Number of bootstrapping runs.")
					  ("alpha,l", bpo::value(&alpha_)->required()->default_value(0.1), "Alpha level of confidence interval.")
					  ("adjust,u", bpo::value(&adjustment_method_)->required()->default_value("benjamini-yekutieli"), "Alpha level of confidence interval.")
					  ("json,j", bpo::value(&json_)->default_value(false)->zero_tokens(), "Output file in .json format (default: .tsv).")
					  ("impact,p", bpo::value(&impact_score_)->required(), "Method that should be used to compute the impact of regulators.")
					  ("confidence-intervals, v", bpo::value(&confidence_interval_)->required(), "Method that should be used to compute confidence intervals.");
	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "ERROR: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}
	
	if(associations_ != "" && checkIfFileExists(desc, associations_)){
	} else if (matrix_ != "" && checkIfFileExists(desc, matrix_) && regulations_ != "" && checkIfFileExists(desc, regulations_)) {
	} else {
		std::cerr << "ERROR: Either a valid 'associations' file or valid 'matrix' and 'regulations' have to be specified. \n";
		desc.print(std::cerr);
		return false;
	}
	
	if(matrix_ != ""){
		if(seed_ == 0){
			std::cerr << "ERROR: Seed needs to be specified. \n";
			desc.print(std::cerr);
			return false;
		}
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

struct MatrixNameDatabase
{
	MatrixNameDatabase(DenseMatrix* matrix) : matrix_(matrix) {}

	std::string operator()(size_t index) { return matrix_->rowNames()[index]; }

	size_t operator()(std::string name) { return matrix_->rowIndex(name); }

	size_t size() { return matrix_->rows(); }

	DenseMatrix* matrix_;
};

struct MapNameDatabase
{
	MapNameDatabase() {}
	
	MapNameDatabase(std::string file) { initialize(file); }

	std::string operator()(size_t index) { return index2name[index]; }

	size_t operator()(std::string name)
	{
		std::map<std::string, size_t>::iterator iter = name2index.find(name);
		if(iter == name2index.end()) {
			index2name.emplace_back(name);
			name2index[name] = index2name.size() - 1;
		}
		return name2index[name];
	}

	size_t size() { return index2name.size(); }
	
	void initialize(std::string file){
		std::cout << "INFO: Initializing name database" << std::endl;
		std::ifstream input(file);
		if(!input) {
			throw GeneTrail::IOError("File (" + file +
			                         ") is not open for reading");
		}

		std::vector<std::string> sline(2);
		for(std::string line; getline(input, line);) {
			boost::split(sline, line, boost::is_any_of(" \t"));
			(*this)(sline[0]);
			(*this)(sline[1]);
		}
	}

	std::map<std::string, size_t> name2index;
	std::vector<std::string> index2name;
};

struct DummyBootstrapper
{
	using Regulation = std::tuple<size_t, size_t, double>;
	DummyBootstrapper() {}

	void create_bootstrap_sample() {}

	template <typename RegulatorImpactScore>
	void perform_bootstrapping_run(std::vector<Regulation>&, bool,
	                               RegulatorImpactScore)
	{
	}
};

template <typename REGGAEAnalysis, typename RegulatorImpactScore>
void run(REGGAEAnalysis& analysis, RegulatorImpactScore impactScore)
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
		std::cerr << "ERROR: Method '" << method_ << "' not known!"
		          << std::endl;
	}
}

template <typename REGGAEAnalysis> void run(REGGAEAnalysis& analysis)
{
	if(impact_score_ == "pearson_correlation") {
		run(analysis, PearsonCorrelation());
	} else if(impact_score_ == "spearman_correlation") {
		run(analysis, SpearmanCorrelation());
	} else if(impact_score_ == "kendall_correlation") {
		run(analysis, KendallCorrelation());
	} else {
		std::cerr << "ERROR: Impact score '" << impact_score_ << "' not known!"
		          << std::endl;
	}
}

template <typename REGGAEAnalysis> void write(REGGAEAnalysis& analysis)
{
	if(json_) {
		analysis.writeJSONResults(out_, alpha_, confidence_interval_);
	} else {
		analysis.writeResults(out_, alpha_, confidence_interval_);
	}
}

int runGeneExpressionAnalysis()
{
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
	MatrixNameDatabase name_database(&matrix);
	RegulationFileParser<MatrixNameDatabase, double> parser(name_database, tset,
	                                                        regulations_, 0.0);
	RegulationFile<double>& regulationFile = parser.getRegulationFile();
	RegulationBootstrapper<double> bootstrapper(&matrix, seed_);
	RegulatorGeneAssociationEnrichmentAnalysis<RegulationBootstrapper<double>,
	                                           MatrixNameDatabase, double>
	    analysis(sorted_targets, regulationFile, bootstrapper, name_database,
	             useAbsoluteValues_, bootstrapping_runs_);
	run(analysis);

	std::cout << "INFO: Adjusting p-values" << std::endl;
	analysis.adjustPValues(adjustment_method_);

	std::cout << "INFO: Writing results" << std::endl;
	write(analysis);

	return 0;
}

int runAssociationAnalysis()
{
	std::cout << "INFO: Parsing test set" << std::endl;
	GeneSetReader reader;
	GeneSet test_set = reader.readScoringFile(scores_);
	std::vector<std::string> sorted_target_names =
	    test_set.getSortedIdentifier(decreasingly_);
	MapNameDatabase name_database(associations_);
	std::vector<size_t> sorted_targets;
	sorted_targets.reserve(sorted_target_names.size());
	for(const std::string& s : sorted_target_names) {
		sorted_targets.emplace_back(name_database(s));
	}

	std::cout << "INFO: Parsing regulations" << std::endl;
	std::unordered_set<size_t> tset(sorted_targets.begin(),
	                                sorted_targets.end());
	RegulationFileParser<MapNameDatabase, double> parser(name_database, tset, associations_, 0.0);
	RegulationFile<double>& regulationFile = parser.getRegulationFile();
	DummyBootstrapper bootstrapper;
	RegulatorGeneAssociationEnrichmentAnalysis<DummyBootstrapper,
	                                           MapNameDatabase, double>
	    analysis(sorted_targets, regulationFile, bootstrapper, name_database,
	             useAbsoluteValues_, bootstrapping_runs_);
	run(analysis);

	std::cout << "INFO: Adjusting p-values" << std::endl;
	analysis.adjustPValues(adjustment_method_);

	std::cout << "INFO: Writing results" << std::endl;
	write(analysis);
	return 0;
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

	if(matrix_ != "") {
		return runGeneExpressionAnalysis();
	} else {
		return runAssociationAnalysis();
	}
}
