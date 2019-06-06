#include <boost/program_options.hpp>
#include <genetrail2/core/Exception.h>
#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/NameDatabases.h>
#include <genetrail2/core/TextFile.h>
#include <genetrail2/core/MatrixHTest.h>
#include <genetrail2/core/Scores.h>

#include <genetrail2/regulation/RegulatorGeneAssociationEnrichmentAnalysis.h>
#include <genetrail2/regulation/RegulatorGeneAssociationEnrichmentAlgorithms.h>
#include <genetrail2/regulation/RegulationFileParser.h>
#include <genetrail2/regulation/RegulationBootstrapperMicro.h>
#include <genetrail2/regulation/RegulationBootstrapper.h>
#include <genetrail2/regulation/RegulatorAssociationScore.h>
#include <genetrail2/regulation/RTINetwork.h>

#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_set>
#include <iostream>
#include "common.h"
#include "../matrixTools.h"

using namespace GeneTrail;

namespace bpo = boost::program_options;

std::string scores_, associations_, matrix_, regulations_, out_, method_, adjustment_method_,groups_,
    impact_score_, confidence_interval_, network_, output_score_file_,matrix_micro_;
bool normalize_scores_, useAbsoluteValues_, decreasingly_, perturbation_, json_, fill_blanks_, sort_correlations_decreasingly_, normalize_impact_scores_,fold_change_;
size_t seed_, bootstrapping_runs_, max_regulators_per_target_ = 0;
double alpha_;

MatrixReaderOptions matrixOptions;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
					  ("scores,s", bpo::value(&scores_)->required(), "A whitespace separated file containing deregulated targets. (test set)")
					  ("associations,t", bpo::value(&associations_)->default_value(""), "A whitespace separated file regulator, target and association scores. (only needed if '--regulations' and '--matrix' options are not provided)") 
					  ("matrix,x", bpo::value(&matrix_)->default_value(""), "A whitespace separated file containing expression values for all genes. (Only needed if '--matrix' is used)")
					  ("matrix-micro,f",bpo::value(&matrix_micro_)->default_value(""),"A whitespace separated file containing the expression values for all mircoRNAs.")   
					  ("groups,g", bpo::value<std::string>(&groups_), "File containing two lines specifying which colnames belong to which group. First line corresponds to test samples, second line to reference samples. (Only needed if matrices are given.)")
					  ("fold-change,i", bpo::value(&fold_change_)->default_value(false)->zero_tokens(), "Should the fold change of miRNAs be calculated in order to exclude ones with small fold change?")
					  ("no-row-names,w", bpo::value<bool>(&matrixOptions.no_rownames)->default_value(false)->zero_tokens(), "Does the matrix file contain row names? (Only needed if '--matrix' is used)")
					  ("no-col-names,c", bpo::value<bool>(&matrixOptions.no_colnames)->default_value(false)->zero_tokens(), "Does the matrix file contain column names? (Only needed if '--matrix' is used)")
					  ("add-col-name,n", bpo::value<bool>(&matrixOptions.additional_colname)->default_value(false)->zero_tokens(), "Does the matrix file contain two rows for column names? (Only needed if '--matrix' is used)")
					  ("decreasingly,d", bpo::value(&decreasingly_)->default_value(false)->zero_tokens(), "Should the testset be sorted decreasingly? (default: increasingly)")
					  ("regulations,r", bpo::value(&regulations_)->default_value(""), "A whitespace separated file containing regulator and target. (RTI database)")
					  ("abs,a", bpo::value(&useAbsoluteValues_)->default_value(false)->zero_tokens(), "Should absolute values be used for association scores?")
					  ("sort-rtis-decreasingly,e", bpo::value(&sort_correlations_decreasingly_)->default_value(false)->zero_tokens(), "Should the association scores for each target be sorted decreasingly? (default: increasingly)")
					  ("method,m", bpo::value(&method_)->required(), "The method that should be applied (ks-test, wrs-test).")
					  ("output,o", bpo::value(&out_)->required(), "Output prefix for text files.")
					  ("seed,e", bpo::value(&seed_)->default_value(0), "Random seed used for pertubation.")
					  ("bootstrap,b", bpo::value(&bootstrapping_runs_)->default_value(0), "Number of bootstrapping runs.")
					  ("alpha,l", bpo::value(&alpha_)->required()->default_value(0.1), "Alpha level of confidence interval.")
					  ("adjust,u", bpo::value(&adjustment_method_)->required()->default_value("benjamini-yekutieli"), "Method for multiple testing correction. (default: benjamini-yekutieli)")
					  ("json,j", bpo::value(&json_)->default_value(false)->zero_tokens(), "Output file in .json format (default: .tsv).")
					  ("impact,p", bpo::value(&impact_score_)->default_value("pearson_correlation"), "Method that should be used to compute the impact of regulators. (pearson_correlation, spearman_correlation, kendall_correlation, network_formular)")
					  ("normalize-association-scores,k", bpo::value(&normalize_scores_)->default_value(false)->zero_tokens(), "Should the association scores be normalized? (default: false) (only used if impact score is calculated '--impact/--matrix')")
					  ("confidence-intervals,v", bpo::value(&confidence_interval_), "Method that should be used to compute confidence intervals. (percentile, bca)")
					  ("max-regulator-per-target,y", bpo::value(&max_regulators_per_target_), "The maximum number of regulators that are allowed to influence a regulator. (optional)")
					  ("fill-blanks,f", bpo::value(&fill_blanks_)->default_value(false)->zero_tokens(), "Should blanks be introduced if a target has a fewer regulators than all other targets. (optional)")
					  ("rti-network-dir,z", bpo::value(&network_)->default_value(""), "Output prefix for RTI network. (optional)")
					  ("create-score-file,i", bpo::value(&output_score_file_)->default_value(""), "Filename for scorefile that should be generated based on the REGGAE results. This is needed to build the RTI network iteratively. (optional)");
	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "ERROR: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}
	
	if(impact_score_ == "network_formular" && (associations_== "" || !checkIfFileExists(desc, associations_))){
	  return false;
	}
	if(fold_change_ && (groups_ == "" || !checkIfFileExists(desc,groups_))){
	  return false;
	}
	if(associations_ != "" && checkIfFileExists(desc, associations_)){
	} else if (regulations_ != "" && checkIfFileExists(desc, regulations_) && matrix_micro_ != "" && checkIfFileExists(desc, matrix_micro_) && matrix_ != "" && checkIfFileExists(desc, matrix_)&& groups_ != "" && checkIfFileExists(desc, groups_)) {
	} else {
		std::cerr << "ERROR: Either a valid 'associations' file or valid 'matrix', 'microRNA matrix' and 'regulations' have to be specified. \n";
		desc.print(std::cerr);
		return false;
	}
	
	if(matrix_ != "" && bootstrapping_runs_ != 0){
		if(seed_ == 0){
			std::cerr << "ERROR: Seed needs to be specified. \n";
			desc.print(std::cerr);
			return false;
		}
	}
	return true;
}

struct DummyBootstrapper
{
	using Regulation = std::tuple<size_t, size_t, double>;
	DummyBootstrapper() {}

	void create_bootstrap_sample() {}

	template <typename RegulatorImpactScore>
	void perform_bootstrapping_run(RegulationFile<double>, std::vector<Regulation>& regulations,
								   bool, bool use_absolute_values,
	                               bool sort_decreasingly, RegulatorImpactScore)
	{
		// Sort values decreasingly
		if(use_absolute_values) {
			if (sort_decreasingly) {
				std::sort(regulations.begin(), regulations.end(),
						[](const Regulation& a, const Regulation& b) {
					return std::abs(std::get<2>(a)) > std::abs(std::get<2>(b));
				});
			} else {
				std::sort(regulations.begin(), regulations.end(),
						[](const Regulation& a, const Regulation& b) {
					return std::abs(std::get<2>(a)) < std::abs(std::get<2>(b));
				});
			}
		} else {
			if (sort_decreasingly) {
				std::sort(regulations.begin(), regulations.end(),
						[](const Regulation& a, const Regulation& b) {
					return std::get<2>(a) > std::get<2>(b);
				});
			} else {
				std::sort(regulations.begin(), regulations.end(),
						[](const Regulation& a, const Regulation& b) {
					return std::get<2>(a) < std::get<2>(b);
				});
			}
		}	
	}
};

template <typename REGGAEAnalysis, typename RegulatorImpactScore>
std::vector<RegulatorEffectResult> run(REGGAEAnalysis& analysis,
                                       RegulatorImpactScore impactScore)
{
	std::vector<RegulatorEffectResult> results;
	if(method_ == "ks-test") {
		std::cout << "INFO: Performing KS-test" << std::endl;
		results = analysis.run(KSTest(), impactScore);
	} else if(method_ == "wrs-test") {
		std::cout << "INFO: Performing WRS-test" << std::endl;
		results = analysis.run(WRSTest(), impactScore);
	} else if(method_ == "em-test") {
		std::cout << "INFO: Performing EM-test" << std::endl;
		results = analysis.run(EMTest(), impactScore);
	} else {
		std::cerr << "ERROR: Method '" << method_ << "' not known!"
		          << std::endl;
	}
	return results;
}

template <typename REGGAEAnalysis>
std::vector<RegulatorEffectResult> run(REGGAEAnalysis& analysis)
{
	std::vector<RegulatorEffectResult> results;
	if(impact_score_ == "pearson_correlation") {
		results = run(analysis, PearsonCorrelation());
	} else if(impact_score_ == "spearman_correlation") {
		results = run(analysis, SpearmanCorrelation());
	} else if(impact_score_ == "kendall_correlation") {
		results = run(analysis, KendallCorrelation());
	} else if(impact_score_ == "network_formular") {
		 results = run(analysis, NetworkFormular());
	} else {
		std::cerr << "ERROR: Impact score '" << impact_score_ << "' not known!"
		          << std::endl;
	}
	return results;
}

void calculate_bootstrap_parameters(std::vector<RegulatorEffectResult>& results,
                                    double alpha, const std::string& ci_method)
{
	for(auto& res : results) {
		if(res.name != ""){
			res.calculate_bootstrap_parameters(alpha, ci_method);
		}
	}
}

void write_score_file(std::vector<RegulatorEffectResult>& results, const std::string& file) {
	std::ofstream out;
	out.open(file);
	for(auto& r : results) {
		if(r.corrected_p_value < 0.05) {
        	out << r.name << "\t" << r.score << std::endl;
		}
	}
	out.close();
}


int runAssociationAnalysis()
{	
	std::cout << "INFO: Parsing test set" << std::endl;
	GeneSetReader reader;
	GeneSet test_set = reader.readScoringFile(scores_);
	std::vector<std::string> sorted_target_names = test_set.getSortedIdentifier(decreasingly_);
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
	if(max_regulators_per_target_ == 0){
		max_regulators_per_target_  = regulationFile.maxNumberOfRegulators();
	}

	DummyBootstrapper bootstrapper;
	RegulatorGeneAssociationEnrichmentAnalysis<DummyBootstrapper,
	                                           MapNameDatabase, double>
	    analysis(sorted_targets, regulationFile, bootstrapper, name_database, normalize_scores_,
	             useAbsoluteValues_, sort_correlations_decreasingly_, fill_blanks_, bootstrapping_runs_, max_regulators_per_target_);
	std::vector<RegulatorEffectResult> results = run(analysis);

	std::cout << "INFO: Adjusting p-values" << std::endl;
	adjustPValues(results, adjustment_method_);
	
	std::cout << "INFO: Calculating confidence intervals" << std::endl;
	calculate_bootstrap_parameters(results, alpha_, confidence_interval_);

	std::cout << "INFO: Writing results" << std::endl;
	write(results, out_, json_);
	if(network_ != "") {
		RTINetwork<MapNameDatabase, double> network(name_database, results, regulationFile);
		network.build_network();
		network.save(network_);
	}

	if(output_score_file_ != "") {
		write_score_file(results, output_score_file_);
	}
	return 0;
}


int matrixAnalysis()
{
    //read in gene matrix
    std::cout << "INFO: Parsing data matrix" << std::endl;
    DenseMatrix matrix(0, 0);
    try {
	  matrix = readDenseMatrix(matrix_, matrixOptions);
	} catch(const IOError& e) {
	  std::cerr << "ERROR: Could not open input data matrix for reading."<< std::endl;
	  return -4;
	}
	
    //read in microRNA matrix
    std::cout << "INFO: Parsing microRNA matrix" << std::endl;
    DenseMatrix microMatrix(0, 0);
    try {
      microMatrix = readDenseMatrix(matrix_micro_, matrixOptions);
    } catch(const IOError& e) {
      std::cerr << "ERROR: Could not open input micoRNA matrix for reading."<< std::endl;
      return -5;
    }
      TextFile t(groups_, ",", std::set<std::string>());
      std::vector<std::string> sample, ref;
	try {
		sample = t.read();
		ref = t.read();
	} catch(const IOError& e) {
		std::cerr << "ERROR: Could not read from group file " << groups_ << std::endl;
		return -5;
	}
	//check fold change of miRNAs and eliminate miRNAS with small/high foldChange
	if(fold_change_ ){
	//create miRNA set with all computed scores
	  auto subset = splitMatrix(microMatrix, sample, ref);

	  //calculate scores for miRNA
	  MatrixHTest htest;
	  Scores mirna_set = htest.test("log-mean-fold-quotient", std::get<0>(subset), std::get<1>(subset));
	  std::vector<unsigned int> vec_indices;
	  EntityDatabase db = mirna_set.getEntityDatabase();
	  for(size_t i  = 0; i < mirna_set.size(); ++i){
	    if(mirna_set[i].score() >= -0.2 && mirna_set[i].score() <= 0.2){//*mirna_set[i].score() == 0 )
	      vec_indices.push_back(microMatrix.rowIndex(mirna_set[i].name(db)));
	    }
	  }
	  
	  microMatrix.removeRows(vec_indices);
    }
    //sort columns of matrices in same order
    std::vector<std::string> orientationControl;
    std::vector<std::string> orientationDisease;


    //check if colomn exists in both matrices, if yes save the indices of both occurrences into orientation
    if(matrix.cols() > microMatrix.cols()){
         for(size_t c = 0; c < microMatrix.cols(); ++c){
	    if(matrix.hasCol(microMatrix.colName(c))){
	      auto it = find(ref.begin(),ref.end(),microMatrix.colName(c));
	      if(it != ref.end()){
		  orientationControl.push_back(microMatrix.colName(c)); 
	      }else{
		  orientationDisease.push_back(microMatrix.colName(c)); 
	      }
	    }
	}
    }else{
      for(size_t c = 0; c < matrix.cols(); ++c){
	    if(microMatrix.hasCol(matrix.colName(c))){
	      auto it = find(ref.begin(),ref.end(),matrix.colName(c));
	      if(it != ref.end()){
		  orientationControl.push_back(matrix.colName(c)); 
	      }else{
		  orientationDisease.push_back(matrix.colName(c)); 
	      }
	    }
	}
    }
    size_t change =orientationDisease.size()+1;
    for(size_t s = 0; s < orientationControl.size(); ++s){
      orientationDisease.push_back(orientationControl[s]);
    }
    microMatrix.sortRows(orientationDisease);
    matrix.sortRows(orientationDisease);

    
   
    //read in scoring file and create the test set with the genes in correct order (decreasingly/increasingly)
    std::cout << "INFO: Parsing test set" << std::endl;
    GeneSetReader reader;
    GeneSet test_set = reader.readScoringFile(scores_); 
    std::vector<std::string> sorted_target_names = test_set.getSortedIdentifier(decreasingly_);
    
    //sorted_target contains the rownumbers of the matrix in the order of the sorted target names
    std::vector<size_t> sorted_targets = translate_test_set(&matrix, sorted_target_names);
    std::unordered_set<size_t> tset(sorted_targets.begin(), sorted_targets.end());
    MatrixNameDatabase name_database_matrix(&matrix);
    MatrixNameDatabase name_database_micro(&microMatrix);  

    
    std::cout << "INFO: Parsing regulations" << std::endl;
    RegulationFileParser<MatrixNameDatabase, double> parser(name_database_matrix, name_database_micro, tset, regulations_, 0.0);
    RegulationFile<double>& regulationFile = parser.getRegulationFile();
    
    if(max_regulators_per_target_ == 0){
	  max_regulators_per_target_  = regulationFile.maxNumberOfRegulators();
    }
    
    RegulationBootstrapperMicro<double> bootstrapper(&matrix,&microMatrix, seed_,change);
    RegulatorGeneAssociationEnrichmentAnalysis<RegulationBootstrapperMicro<double>,MatrixNameDatabase, double> analysis(sorted_targets, regulationFile, bootstrapper, name_database_micro, normalize_scores_,
	             useAbsoluteValues_, sort_correlations_decreasingly_, fill_blanks_, bootstrapping_runs_, max_regulators_per_target_);
    std::vector<RegulatorEffectResult> results = run(analysis);
    
    std::cout << "INFO: Adjusting p-values" << std::endl;
    adjustPValues(results, adjustment_method_);
	
    std::cout << "INFO: Calculating confidence intervals" << std::endl;
    calculate_bootstrap_parameters(results, alpha_, confidence_interval_);

    std::cout << "INFO: Writing results" << std::endl;
    write(results, out_, json_);

    if(output_score_file_ != "") {
		write_score_file(results, output_score_file_);
    }

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
	if(associations_ != "") {
		return runAssociationAnalysis();

	}else{
		return matrixAnalysis();
	}

}
