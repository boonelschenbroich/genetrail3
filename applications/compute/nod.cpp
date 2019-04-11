 
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
#include <genetrail2/core/WilcoxonSignedRankSumTest.h>
#include <genetrail2/core/ExactQuantileTest.h>


#include <genetrail2/regulation/RegulatorGeneAssociationEnrichmentAnalysis.h>
#include <genetrail2/regulation/RegulatorGeneAssociationEnrichmentAlgorithms.h>
#include <genetrail2/regulation/RegulationFileParser.h>
#include <genetrail2/regulation/RegulationBootstrapperMicro.h>
#include <genetrail2/regulation/RegulationBootstrapper.h>
#include <genetrail2/regulation/RegulatorAssociationScore.h>
#include <genetrail2/regulation/RTINetwork.h>
#include <genetrail2/regulation/RegulatorEffectResult.h>
#include <genetrail2/regulation/RegulationFile.h>

#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/prettywriter.h>


#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <algorithm>
#include <vector>
#include <map>
#include <chrono>
#include <unordered_set>
#include "../matrixTools.h"
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>


using namespace GeneTrail;
using namespace std::chrono;

namespace bpo = boost::program_options;


bool checkIfFileExists(bpo::options_description& desc, std::string fname)
{
	std::ifstream infile(fname);
	if(!infile.good()) {
		std::cerr << "ERROR: File '" + fname + "' does not exist. \n";
		desc.print(std::cerr);
		return false;
	}
	return true;
}

std::string scores_,scoresMirna_, matrix_, regulations_, out_, adjustment_method_,matrix_micro_,groups_;
bool normalize_scores_, useAbsoluteValues_, decreasingly_, perturbation_, json_, fill_blanks_, sort_correlations_decreasingly_, normalize_impact_scores_;
size_t max_regulators_per_target_,seed_ = 0;
double threshold_ = -0.6;

MatrixReaderOptions matrixOptions;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
					  ("scores,s", bpo::value(&scores_)->required(), "A whitespace separated file containing deregulated targets. (test set)")
					  ("scoresMirna", bpo::value(&scoresMirna_)->required(), "A whitespace separated file containing deregulated miRNAs. (test set for miRNAs)")
					  ("threshold,t", bpo::value(&threshold_)->default_value(-0.6), "A threshold for the correlation. default is -0.6")
					  ("matrix,x", bpo::value(&matrix_)->default_value(""), "A whitespace separated file containing expression values for all genes. (Only needed if '--matrix' is used)")
					  ("matrix-micro,h",bpo::value(&matrix_micro_)->default_value(""),"A whitespace separated file containing the expression values for all mircoRNAs.")
					  ("no-row-names,w", bpo::value<bool>(&matrixOptions.no_rownames)->default_value(false)->zero_tokens(), "Does the matrix file contain row names? (Only needed if '--matrix' is used)")
					  ("no-col-names,c", bpo::value<bool>(&matrixOptions.no_colnames)->default_value(false)->zero_tokens(), "Does the matrix file contain column names? (Only needed if '--matrix' is used)")
					  ("add-col-name,n", bpo::value<bool>(&matrixOptions.additional_colname)->default_value(false)->zero_tokens(), "Does the matrix file contain two rows for column names? (Only needed if '--matrix' is used)")
					  ("decreasingly,d", bpo::value(&decreasingly_)->default_value(false)->zero_tokens(), "Should the testset be sorted decreasingly? (default: increasingly)")
					  ("regulations,r", bpo::value(&regulations_)->default_value(""), "A whitespace separated file containing regulator and target. (RTI database)")
					  ("groups,g", bpo::value<std::string>(&groups_), "File containing two lines specifying which colnames belong to which group. First line corresponds to test samples, second line to reference samples. (Only needed if matrices are given.)")
					  ("abs,a", bpo::value(&useAbsoluteValues_)->default_value(false)->zero_tokens(), "Should absolute values be used for association scores?")
					  ("seed,e", bpo::value(&seed_)->default_value(0), "Random seed used for pertubation.")
					  ("sort-rtis-decreasingly,e", bpo::value(&sort_correlations_decreasingly_)->default_value(false)->zero_tokens(), "Should the association scores for each target be sorted decreasingly? (default: increasingly)")
					  ("output,o", bpo::value(&out_)->required(), "Output prefix for text files.")
					  ("adjust,u", bpo::value(&adjustment_method_)->required()->default_value("benjamini-yekutieli"), "Method for multiple testing correction. (default: benjamini-yekutieli)")
					  ("json,j", bpo::value(&json_)->default_value(false)->zero_tokens(), "Output file in .json format (default: .tsv).")
					  ("normalize-association-scores,k", bpo::value(&normalize_scores_)->default_value(false)->zero_tokens(), "Should the association scores be normalized? (default: false) (only used if impact score is calculated '--impact/--matrix')")
					  ("max-regulator-per-target,y", bpo::value(&max_regulators_per_target_), "The maximum number of regulators that are allowed to influence a regulator. (optional)")
					  ("fill-blanks,f", bpo::value(&fill_blanks_)->default_value(false)->zero_tokens(), "Should blanks be introduced if a target has a fewer regulators than all other targets. (optional)")
					  ;
	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "ERROR: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}
	
	if(scores_ == "" || checkIfFileExists(desc, regulations_)|| regulations_== "" || checkIfFileExists(desc, regulations_) || matrix_ =="" ||matrix_micro_ == "" || checkIfFileExists(desc, matrix_) || checkIfFileExists(desc, matrix_micro_)|| groups_== "" || checkIfFileExists(desc, groups_)){
	  return true;
	}
	return false;
	
	



}

//write results to json file
void writeResults(std::vector<RegulatorEffectResult>& results_,
                      const std::string& out_)
{
	std::vector<RegulatorEffectResult> results;
	for(size_t i = 0; i < results_.size(); ++i) {
		if(results_[i].name == "") {
			continue;
		}
		results.emplace_back(results_[i]);
	}
	std::sort(
	    results.begin(), results.end(),
	    [](const RegulatorEffectResult& lhs, const RegulatorEffectResult& rhs) {
		    return lhs.score == rhs.score ? abs(lhs.score) > abs(rhs.score)
		                                      : lhs.score > rhs.score;
		});

	rapidjson::StringBuffer sb;
	rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(sb);

	size_t rank = 1;
	writer.StartArray();
	for(RegulatorEffectResult& res : results) {
		res.rank = rank;
		res.serializeJSON(writer);
		++rank;
	}
	writer.EndArray();

	std::ofstream out;
	out.open(out_);
	if(out.is_open()) {
		out << sb.GetString();
	} else {
		std::cerr << "Could not open file: " << out_ << "\n";
	}
	out.close();
}





std::vector<size_t> translate_test_set(const DenseMatrix* matrix,
                   const std::vector<std::string>& test_set_names)
{
	std::vector<size_t> test_set;
	test_set.reserve(test_set_names.size());
	for(const std::string name : test_set_names) {
		test_set.emplace_back(matrix->rowIndex(name));

	}
	return test_set;
}

//count which miRNAs how often regulate a gene independently
std::map<size_t,size_t> findNODValues(std::vector<std::vector<size_t>>& regulatorsPerGene){
    std::map<size_t,size_t> nodValues;
    for(size_t gene = 0; gene < regulatorsPerGene.size(); ++gene){
      if(regulatorsPerGene[gene].size() == 1){
	if(nodValues.find(regulatorsPerGene[gene].front()) == nodValues.end()){
	  nodValues[regulatorsPerGene[gene].front()] = 1;
	}else{
	  size_t old = nodValues[regulatorsPerGene[gene].front()];
	  nodValues[regulatorsPerGene[gene].front()] = old+1;
	}
      }
    }
    return nodValues;
}
                
        
void adjustPValues(std::vector<RegulatorEffectResult>& results_,
                   const std::string& method)
{
	std::vector<std::pair<size_t, double>> p_values;
	p_values.reserve(results_.size());
	for(size_t i = 0; i < results_.size(); ++i) {
		if(results_[i].name == "") {
			continue;
		}
		p_values.emplace_back(std::make_pair(i, results_[i].p_value));
	}

	// Adjust p-values
	p_values = pvalue::adjustPValues(p_values, pvalue::get_second(),
	                                 pvalue::getCorrectionMethod(method).get());

	// Update p-values
	for(const auto& pair : p_values) {
		results_[pair.first].corrected_p_value = pair.second;
	}
}  

//deletes all rows/miRNAs in matrix that do not exist in scoresMirna file
void adjustMatrix(std::string& file, DenseMatrix* microMatrix){
	std::ifstream input(file);
	if(!input) {
		throw GeneTrail::IOError("File (" + file +
		                         ") is not open for reading");
	}
	std::vector<std::string> sline(2);
	std::vector<unsigned int> indicesKeep;
	for(std::string line; getline(input, line);) {
	  boost::split(sline, line, boost::is_any_of(" \t"), boost::token_compress_on);
	  if(microMatrix->hasRow(sline[0])){
	    indicesKeep.push_back(microMatrix->rowIndex(sline[0]));
	  }
	}
	std::vector<unsigned int> indicesRemove;
	for(unsigned int i = 0; i < microMatrix->rows(); ++i) {
	  if(std::find(indicesKeep.begin(),indicesKeep.end(),i) == indicesKeep.end()){
	    indicesRemove.push_back(i);
	  }
	}
	microMatrix->removeRows(indicesRemove);
}

//calls correlation and nod count
int NODAnalysis()
{
    auto start = high_resolution_clock::now();

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
    //remove rows in microMatrix that are not in scoresMirna
    adjustMatrix(scoresMirna_,&microMatrix);

    std::vector<std::string> orientationControl;
    std::vector<std::string> orientationDisease;

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
    std::cout << "INFO: Parsing scoring file" << std::endl;
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
    //calculate pearson correlation for valid microRNA target pairs
    size_t biggest_target_idx = 0;
    for(size_t targetname : sorted_targets) {
	if(!regulationFile.checkTarget(targetname)) {
		continue;
	}
	biggest_target_idx = std::max(targetname, biggest_target_idx);
	auto& regulations = regulationFile.target2regulations(targetname);
	bootstrapper.perform_bootstrapping_run(regulationFile, regulations, normalize_scores_, useAbsoluteValues_,sort_correlations_decreasingly_, PearsonCorrelation());
    }
    std::cout << "INFO: Calculating NOD values" << std::endl;

    std::vector<std::vector<size_t>> regulatorsPerGene(biggest_target_idx+1);
    for(size_t t : sorted_targets) {
	if(!regulationFile.checkTarget(t)) {
 		continue;
 	}
 	std::vector<std::tuple<size_t, size_t, double>>& regulators = regulationFile.target2regulations(t);
 	for(auto& r : regulators){
 	  if(!regulationFile.checkRegulator(std::get<0>(r)) || !regulationFile.checkTarget(std::get<1>(r))) {
 		continue;
 	  }
 	  if(std::get<2>(r) <= threshold_){
 	    regulatorsPerGene[std::get<1>(r)].emplace_back(std::get<0>(r));
 	  }
 	}
    }
    std::map<size_t,size_t> nodValues = findNODValues(regulatorsPerGene);
    
    //find largest index of regulator 
    size_t biggest_regulator_idx = 0;
    for(size_t t : sorted_targets) {
	if(!regulationFile.checkTarget(t)) {
 		continue;
 	}
 	std::vector<std::tuple<size_t, size_t, double>>& regulators = regulationFile.target2regulations(t);
 	for(auto& allRegulators : regulators){
 	  size_t regulator_i = std::get<0>(allRegulators);
 	  biggest_regulator_idx = std::max(regulator_i, biggest_regulator_idx);
 	}
    }
    std::vector<double> onlyNOD;
    for(auto& val : nodValues){
	onlyNOD.push_back(val.second);
    }
    std::sort(onlyNOD.begin(), onlyNOD.end());
 
    std::vector<RegulatorEffectResult> results;
    bool empty = results.size() == 0;
    if(empty || results.size() < biggest_regulator_idx + 1) {
 	results.resize(biggest_regulator_idx + 1);
    }
    boost::math::normal dist(0,1);
    for(auto& val : nodValues){
      if(empty){
 	    results[val.first].name = name_database_micro(val.first);
 	    results[val.first].hits = 0;
 	    results[val.first].score = val.second;
	    WilcoxonSignedRankSumTest<double> Wtest(1e-4,val.second);
	    double zscore = Wtest.test(onlyNOD.begin(),onlyNOD.end());
	    std::cout << "INFO: Calculating p-value for "<< results[val.first].name << std::endl;
	    double pvalue = boost::math::cdf(dist, zscore);
	    results[val.first].p_value = pvalue;
 	}
    }
   
    std::cout << "INFO: Adjusting p-values" << std::endl;
    adjustPValues(results, adjustment_method_);

    std::cout << "INFO: Writing results" << std::endl;
    writeResults(results, out_);
    auto stop = high_resolution_clock::now();
    auto dur = duration_cast<microseconds>(stop-start);
    std::cout << dur.count() << std::endl;
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

	return NODAnalysis();

}
