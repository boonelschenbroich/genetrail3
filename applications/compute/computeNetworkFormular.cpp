#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <tuple>

#include <boost/program_options.hpp>

#include <genetrail2/regulation/RegulationFileParser.h>
#include <genetrail2/regulation/RegulationFileWriter2.h>

#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/NameDatabases.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/TextFile.h>
#include <genetrail2/core/MatrixHTest.h>
#include <genetrail2/core/Scores.h>
#include <genetrail2/core/EntityDatabase.h>
#include <chrono>


using namespace std::chrono;

#include "../matrixTools.h"

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string expr1 = "", output = "", method = "", groups = "", regulations = "",scores = "",scores_mirna = "";
bool binary = false, fold_change_ = false;

MatrixReaderOptions matrixOptions;


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


bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("expression-matrix,e", bpo::value<std::string>(&expr1), "Name of a text file containing expression values for miRNAs as a matrix.")
		("output,o", bpo::value<std::string>(&output)->required(), "Name of the output file.")
		("groups,g", bpo::value<std::string>(&groups), "File containing two lines specifying which colnames belong to which group.")
		("no-row-names,r", bpo::value<bool>(&matrixOptions.no_rownames)->default_value(false)->zero_tokens(), "Does the file contain row names.")
		("no-col-names,c", bpo::value<bool>(&matrixOptions.no_colnames)->default_value(false)->zero_tokens(), "Does the file contain column names.")
		("add-col-name,a", bpo::value<bool>(&matrixOptions.additional_colname)->default_value(false)->zero_tokens(), "File containing two lines specifying which rownames belong to which group.")
		("method,m", bpo::value<std::string>(&method), "Method used for scoring.")
		("fold-change,f", bpo::value(&fold_change_)->default_value(false)->zero_tokens(), "Should the fold change of miRNAs be calculated in order to exclude ones with small fold change?")
		("regulations,u", bpo::value<std::string>(&regulations)->required(), "A whitespace separated file containing miRNA-target pairs.")
		("scores,s", bpo::value(&scores)->required(), "A whitespace separated file containing deregulated targets.")
		("scores-mirna,i", bpo::value(&scores_mirna), "A whitespace separated file containing deregulated miRNAs.");



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
	if(fold_change_ && groups != "" && checkIfFileExists(desc,groups) && expr1 != ""  && checkIfFileExists(desc,expr1) && method != ""){
	  return true;
	}else if(expr1 != ""  && checkIfFileExists(desc,expr1) && method != "" && groups != "" && checkIfFileExists(desc,groups)){
	  return true;
	}else if (scores_mirna != "" && checkIfFileExists(desc,scores_mirna)){
	  return true;	  
	}
	return false;
}


double calculateSum(std::vector<std::tuple<size_t,size_t,double>> regulations, MapNameDatabase& name_database, std::map<std::string, double>& map){
	double sum = 0;
	for(std::tuple<size_t,size_t,double> tuple1 : regulations){
	      std::string name_target_x = name_database(std::get<1>(tuple1));

	      if(map.find(name_target_x) == map.end()){
			continue;
		}
		sum += map.at(name_target_x);
	      }
	return sum;
}
	
int runScoreAnalysis(std::string& scores_mirna,std::string& scores,std::string& regulations,std::string& output){
		
		MapNameDatabase name_database(regulations);
	  
		//read in scoring file of genes for the following calculation
		GeneSetReader reader;
		GeneSet test_set = reader.readScoringFile(scores); 
		std::vector<std::string> target_names = test_set.getIdentifier();
		std::vector<size_t> targets;
		targets.reserve(target_names.size());
		for(const std::string& s : target_names) {
			targets.emplace_back(name_database(s));
		}
		std::unordered_set<size_t> tset(targets.begin(), targets.end());

		//create map gene name-score to look up the scores fast
		std::vector<std::pair<std::string,double>> pairs = test_set.getIncreasinglySortedScores();
		std::map<std::string, double> map;
		for(auto val : pairs){
		  map.insert(val);
		}		
		
		//read in scoring file for miRNAs
		GeneSetReader reader_mirna;
		GeneSet test_set_mirna = reader_mirna.readScoringFile(scores_mirna); 
		std::vector<std::string> target_names_mirna = test_set_mirna.getIdentifier();

		//create map mirna name-score to look up the scores fast
		std::vector<std::pair<std::string,double>> pairs_mirna = test_set_mirna.getIncreasinglySortedScores();
		std::map<std::string, double> map_mirna;
		for(auto val_mirna : pairs_mirna){
		  map_mirna.insert(val_mirna);
		}		
		
		//read in the regulation file to know the regulation target pairs
		RegulationFileParser<MapNameDatabase, double> parser(name_database, tset, regulations, 0.0);
		RegulationFile<double>& regulationFile = parser.getRegulationFile();

		std::vector<std::tuple<std::string, std::string, double>> results;
		std::map<std::string, double> mirna_sum;

		for(std::string& name_mirna : target_names_mirna) {

			if(!regulationFile.checkRegulator(name_database(name_mirna))){
			  continue;
			}
			auto& regulations = regulationFile.regulator2regulations(name_database(name_mirna));
			if(regulations.size() >  regulationFile.maxNumberOfTargets()){
			  continue;
			}
			double score_mirna = map_mirna.at(name_mirna);
			double sum = 0;
			for(std::tuple<size_t,size_t,double> tuple : regulations){
			      size_t target = std::get<1>(tuple);

			      std::string name_tar = name_database(target);
			      if(map.find(name_tar) == map.end()){
				continue;
			      }
			      double score_target = map.at(name_tar);
			      if(mirna_sum.find(name_mirna) == mirna_sum.end()){
				sum = calculateSum(regulations,name_database,map);
				mirna_sum[name_mirna] = sum;
			      }else{
				sum = mirna_sum.at(name_mirna);
			      }
			      double x = (score_target/sum)*score_mirna;
			      results.push_back(std::make_tuple(name_mirna,name_tar,x));
			}
		}
		RegulationFileWriter2 writer;
		writer.write(results,output,"\t");
		return 0;
}


int runMatrixAnalysis(std::string& expr1, std::string& method, std::string& groups, std::string& regulations, std::string& scores_, std::string& output){
	  
		auto start = high_resolution_clock::now();

		TextFile t(groups, ",", std::set<std::string>());
		DenseMatrix matrix(0,0);
		
		try {
			matrix = readDenseMatrix(expr1, matrixOptions);
		} catch(const IOError& e) {
			std::cerr << "ERROR: Could not open input data matrix for reading." << std::endl;
			return -4;
		}
		//divide coloumns into reference and sample
		std::vector<std::string> reference, sample;
		try {
			reference = t.read();
			sample = t.read();
		} catch(const IOError& e) {
			std::cerr << "ERROR: Could not read from group file " << groups << std::endl;
			return -5;
		}
		EntityDatabase db;
		//calculate fold changes for miRNA to exclude some miRNAs
		if(fold_change_){
	  
		    auto subset = splitMatrix(matrix, reference, sample);
		    MatrixHTest htest_foldChange;
		    Scores mirna_set_foldChange = htest_foldChange.test("log-mean-fold-quotient", std::get<0>(subset), std::get<1>(subset));
		    std::vector<unsigned int> vec_indices;
		    db = mirna_set_foldChange.getEntityDatabase();
		    for(size_t i  = 0; i < mirna_set_foldChange.size(); ++i){
		      if(/*mirna_set_foldChange[i].score() == 0 ){*/mirna_set_foldChange[i].score() >= -0.2 && mirna_set_foldChange[i].score() <= 0.2){

			vec_indices.push_back(matrix.rowIndex(mirna_set_foldChange[i].name(db)));
		      }
		    }
		    matrix.removeRows(vec_indices);
		}	
		//create miRNA set with all computed scores
		auto subset = splitMatrix(matrix, reference, sample);
		//calculate scores for miRNA
		MatrixHTest htest;
		Scores mirna_set = htest.test(method, std::get<0>(subset), std::get<1>(subset));

		MapNameDatabase name_database(regulations);
	  
		//read in scoring file of genes for the following calculation
		GeneSetReader reader;
		GeneSet test_set = reader.readScoringFile(scores_); 
		std::vector<std::string> target_names = test_set.getIdentifier();

		std::vector<size_t> targets;
		targets.reserve(target_names.size());
		for(const std::string& s : target_names) {
			targets.emplace_back(name_database(s));
		}
		std::unordered_set<size_t> tset(targets.begin(), targets.end());
		//create map gene name-score to look up the scores fast
		std::vector<std::pair<std::string,double>> scores = test_set.getIncreasinglySortedScores();
		std::map<std::string, double> map;
		for(auto val : scores){
		  map.insert(val);
		}		
		//read in the regulation file in order to know the regulation target pairs
		RegulationFileParser<MapNameDatabase, double> parser(name_database, tset, regulations, 0.0);
		RegulationFile<double>& regulationFile = parser.getRegulationFile();
		
		//for every mirna-score calculate the formular to all known targets
		db = mirna_set.getEntityDatabase();
		
		std::vector<std::tuple<std::string, std::string, double>> results;
		std::map<std::string, double> mirna_sum;

		for(size_t i  = 0; i < mirna_set.size(); ++i) {

			Score sc = mirna_set[i];
			if(isnan(sc.score())) {
				std::cerr
				    << "WARNING: NaNs generated during score computation.\n";
				break;
			}

			std::string name_mirna = sc.name(db);
			if(!regulationFile.checkRegulator(name_database(name_mirna))){
			  continue;
			}
			auto& regulations = regulationFile.regulator2regulations(name_database(name_mirna));
			if(regulations.size() >  regulationFile.maxNumberOfTargets()){
			  continue;
			}
			double sum = 0;

			for(std::tuple<size_t,size_t,double> tuple : regulations){
			      size_t target = std::get<1>(tuple);

			      std::string name_tar = name_database(target);
			      if(map.find(name_tar) == map.end()){
				continue;
			      }
			      double score_target = map.at(name_tar);
			      if(mirna_sum.find(name_mirna) == mirna_sum.end()){
				sum = calculateSum(regulations,name_database,map);
				mirna_sum[name_mirna] = sum;
			      }else{
				sum = mirna_sum.at(name_mirna);
			      }
			      double x = (score_target/sum)*sc.score();
			      results.push_back(std::make_tuple(name_mirna,name_tar,x));
			}
		}
		RegulationFileWriter2 writer;
		writer.write(results,output,"\t");
		auto stop = high_resolution_clock::now();
		auto dur = duration_cast<microseconds>(stop-start);
		std::cout << dur.count() << std::endl;
		return 0;
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
	
	if(scores_mirna == ""){
	  return runMatrixAnalysis(expr1,method,groups,regulations,scores,output);
	}else{
	  return runScoreAnalysis(scores_mirna,scores,regulations,output);
	}
}