#include <iostream>
#include <fstream>
#include <functional> 
#include <map>
#include <string> 
#include <vector>
#include <algorithm>
#include <utility>
#include <fstream>
#include <stdlib.h> 

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <genetrail2/core/Matrix.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrixWriter.h>
#include <genetrail2/core/Entropy.h>
#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/MatrixIterator.h>

#include "../matrixTools.h"

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string matrix = "", output = "", transformation = "";

bool use_ranks = false;
DenseMatrix valueMatrix(0,0);
MatrixReaderOptions matrixOptions;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("matrix,m", bpo::value<std::string>(&matrix)->required(), "Name of the matrix file.")
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

void readValueMatrix()
{
	valueMatrix = readDenseMatrix(matrix,matrixOptions);
}

void greedy() {
	std::vector<double> target_scores(valueMatrix.cols());
	std::vector<double> reg_scores(valueMatrix.cols());
	while(true) {
		std::string s;
		std::cout << "Bitte Gene eingeben (Target|Regulator,Regulator;Regulator,...):";
		std::cin >> s;
		if (s == "quit()") {
			break;
		}
		std::vector<std::string> ss;
		boost::split(ss, s, boost::is_any_of("|"));
		std::string target = ss[0];

		std::vector<std::string> regulator_groups;
		boost::split(regulator_groups, ss[1], boost::is_any_of(";"));
		
		// Get values
		RowMajorMatrixIterator<Matrix> ref_it(&valueMatrix, valueMatrix.rowIndex(target));
		std::copy(ref_it->begin(), ref_it->end(), target_scores.begin());

		std::ofstream out;
  		out.open (s + ".R");
		out << "library(\"ggplot2\")" << std::endl;
		std::vector<std::string> curves;

		// Get self information
		std::vector<std::vector<double>> scores;
		scores.emplace_back(target_scores);
		out << target << " <- c(";
		curves.emplace_back(target);
		std::vector<double> own_entropy = std::get<0>(Entropy::conditional_cummulative_entropy_estimator(target_scores, scores));
		for(size_t i=0; i<own_entropy.size(); ++i) {
				own_entropy[i] = (own_entropy.back() - own_entropy[i]) / own_entropy.back();
				out << own_entropy[i];
				if(i < (own_entropy.size()-1)) {
					out << ", ";
				}
		}
		out << ")" << std::endl;

		for(auto regg : regulator_groups) {
				std::vector<std::string> regulators;
				boost::split(regulators, regg, boost::is_any_of(","));

				//Get values for regulators
				std::vector<std::vector<double>> regulator_scores;
				//regulator_scores.emplace_back(target_scores);
				for (auto& reg : regulators) {
						RowMajorMatrixIterator<Matrix> reg_it(&valueMatrix, valueMatrix.rowIndex(reg));
						std::copy(reg_it->begin(), reg_it->end(), reg_scores.begin());
						regulator_scores.emplace_back(reg_scores);
				}

				// Calculate entropy
				std::string tmp = regg;
				boost::replace_all(tmp, ",", ".");
				
				
				/////////////////////////////////////////////////////////////////
				out << tmp << " <- c(";
				curves.emplace_back(tmp);
				std::vector<double> entropy = Entropy::conditional_cummulative_entropy_estimator(target_scores, regulator_scores);
				for(size_t i=0; i<entropy.size(); ++i) {
						entropy[i] = (entropy.back() - entropy[i]) / entropy.back();
						out << entropy[i];
						if(i < (entropy.size()-1)) {
								out << ", ";
						}
				}
				out << ")" << std::endl;

				/////////////////////////////////////////////////////////////////

				std::cout << std::endl;
				std::cout << target << "|" << regg << " - Mean: " << statistic::mean<double>(entropy.begin(), entropy.end()) << std::endl;
				std::cout << target << "|" << regg << " - Median: " << statistic::median<double>(entropy.begin(), entropy.end()) << std::endl << std::endl;
		}

		out << "df <- data.frame(";
		for (auto c : curves) {
				out << c << "=" << c << ",";
		}
		out << " Bins=c(length(" + curves[0] + "):1))" << std::endl;

		out << "pdf(\"" << s << ".pdf\")" << std::endl;
		out << "ggplot(df, aes(Bins))";
		out << + " + geom_line(data=data.frame(X=c(0,length(" << curves[0] << ")), Y=c(0,1)), aes(x=X, y=Y), colour='black')";
		for (auto c : curves) {
				out << " + geom_line(aes(y = " << c << ", colour = \"" << c << "\"))";
		}
		out << std::endl;
		
		out << "ggplot(df, aes(Bins))";
		out << + " + geom_line(data=data.frame(X=c(0,length(" << curves[0] << ")), Y=c(0,1)), aes(x=X, y=Y), colour='black')";
		for (auto c : curves) {
				out << " + geom_line(aes(y = " << c << ", colour = \"" << c << "\"))";
		}
		out << " + xlim(0, length(" << curves[0] << ")/3)";
		out << std::endl;

		out << "ggplot(df, aes(Bins))";
		out << + " + geom_line(data=data.frame(X=c(0,length(" << curves[0] << ")), Y=c(0,1)), aes(x=X, y=Y), colour='black')";
		for (auto c : curves) {
				out << " + geom_line(aes(y = " << c << ", colour = \"" << c << "\"))";
		}
		out << " + xlim(0, 100)";
		out << std::endl;

		out << "ggplot(df, aes(Bins))";
		out << + " + geom_line(data=data.frame(X=c(0,length(" << curves[0] << ")), Y=c(0,1)), aes(x=X, y=Y), colour='black')";
		for (auto c : curves) {
				out << " + geom_line(aes(y = " << c << ", colour = \"" << c << "\"))";
		}
		out << " + xlim(0, 50)";
		out << std::endl;

		out << "ggplot(df, aes(Bins))";
		out << + " + geom_line(data=data.frame(X=c(0,length(" << curves[0] << ")), Y=c(0,1)), aes(x=X, y=Y), colour='black')";
		for (auto c : curves) {
				out << " + geom_line(aes(y = " << c << ", colour = \"" << c << "\"))";
		}
		out << " + xlim(0, 30)";
		out << std::endl;

		out << "dev.off()" << std::endl;

		out.close();
		std::string command = "Rscript '" + s + ".R'";
		system(command.c_str());
	}
}

void kruskall() {
	std::vector<double> target_scores(valueMatrix.cols());
	std::vector<double> reg_scores(valueMatrix.cols());
	while(true) {
		std::string s;
		std::cout << "Bitte Gene eingeben (Target|Regulator,Regulator;Regulator,...):";
		std::cin >> s;
		if (s == "quit()") {
			break;
		}
		std::vector<std::string> ss;
		boost::split(ss, s, boost::is_any_of("|"));
		std::string target = ss[0];

		std::vector<std::string> regulator_groups;
		boost::split(regulator_groups, ss[1], boost::is_any_of(";"));
		
		// Get values
		RowMajorMatrixIterator<Matrix> ref_it(&valueMatrix, valueMatrix.rowIndex(target));
		std::copy(ref_it->begin(), ref_it->end(), target_scores.begin());

		std::ofstream out;
  		out.open ("Kruskall_" + s + ".R");
		out << "library(\"ggplot2\")" << std::endl;
		std::vector<std::string> curves;

		// Get self information
		std::vector<std::vector<double>> scores;
		scores.emplace_back(target_scores);
		out << target << " <- c(";
		curves.emplace_back(target);
		std::vector<double> own_entropy = Entropy::conditional_cummulative_entropy_kruskall_estimator(target_scores, scores);
		for(size_t i=0; i<own_entropy.size(); ++i) {
				own_entropy[i] = (own_entropy.back() - own_entropy[i]) / own_entropy.back();
				out << own_entropy[i];
				if(i < (own_entropy.size()-1)) {
					out << ", ";
				}
		}
		out << ")" << std::endl;

		for(auto regg : regulator_groups) {
				std::vector<std::string> regulators;
				boost::split(regulators, regg, boost::is_any_of(","));

				//Get values for regulators
				std::vector<std::vector<double>> regulator_scores;
				//regulator_scores.emplace_back(target_scores);
				for (auto& reg : regulators) {
						RowMajorMatrixIterator<Matrix> reg_it(&valueMatrix, valueMatrix.rowIndex(reg));
						std::copy(reg_it->begin(), reg_it->end(), reg_scores.begin());
						regulator_scores.emplace_back(reg_scores);
				}

				// Calculate entropy
				std::string tmp = regg;
				boost::replace_all(tmp, ",", ".");
				out << tmp << " <- c(";
				curves.emplace_back(tmp);
				std::vector<double> entropy = Entropy::conditional_cummulative_entropy_kruskall_estimator(target_scores, regulator_scores);
				for(size_t i=0; i<entropy.size(); ++i) {
						entropy[i] = (entropy.back() - entropy[i]) / entropy.back();
						out << entropy[i];
						if(i < (entropy.size()-1)) {
								out << ", ";
						}
				}
				out << ")" << std::endl;

				std::cout << std::endl;
				std::cout << target << "|" << regg << " - Mean: " << statistic::mean<double>(entropy.begin(), entropy.end()) << std::endl;
				std::cout << target << "|" << regg << " - Median: " << statistic::median<double>(entropy.begin(), entropy.end()) << std::endl << std::endl;
		}

		out << "df <- data.frame(";
		for (auto c : curves) {
				out << c << "=" << c << ",";
		}
		out << " Bins=c(length(" + curves[0] + "):1))" << std::endl;

		out << "pdf(\"" << "Kruskall_" << s << ".pdf\")" << std::endl;
		out << "ggplot(df, aes(Bins))";
		out << + " + geom_line(data=data.frame(X=c(0,length(" << curves[0] << ")), Y=c(0,1)), aes(x=X, y=Y), colour='black')";
		for (auto c : curves) {
				out << " + geom_line(aes(y = " << c << ", colour = \"" << c << "\"))";
		}
		out << std::endl;
		out << "ggplot(df, aes(Bins))";
		out << + " + geom_line(data=data.frame(X=c(0,length(" << curves[0] << ")), Y=c(0,1)), aes(x=X, y=Y), colour='black')";
		for (auto c : curves) {
				out << " + geom_line(aes(y = " << c << ", colour = \"" << c << "\"))";
		}
		out << " + xlim(0, length(" << curves[0] << ")/10)";
		out << std::endl;
		out << "dev.off()" << std::endl;

		out.close();
		std::string command = "Rscript 'Kruskall_" + s + ".R'";
		system(command.c_str());
	}
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
	
	try {
		readValueMatrix();
	}	catch(const IOError& e) {
		std::cerr << "ERROR: Could not read from matrix file " << matrix << std::endl;
		return -3;
	}

	while(true) {
		std::cout << "Please select an algorithm (greedy,kruskall): ";
		std::string s;
		std::cin >> s;
		if (s == "greedy") {
			greedy();
		} else if (s == "kruskall") {
			kruskall();
		}
	}

  return 0;
  
}
