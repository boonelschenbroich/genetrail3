#include <genetrail2/core/Category.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/OverRepresentationAnalysis.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/ScoringFile.h>

#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <iostream>
#include <fstream>
#include <cstdint>
#include <utility>
#include <tuple>
#include <map>

using namespace GeneTrail;
namespace bpo = boost::program_options;
namespace bm = boost::math;
using namespace boost::multiprecision;

std::string categories, scores = "", adjustment, out, identifier = "", json, reference;
double significance;
int minimum, maximum;

std::map<std::string,std::vector<std::pair<std::string, double>>> all_results;
std::map<std::string,std::map<std::string, std::string>> all_name2reference;
std::map<std::string,std::map<std::string, double>> all_name2ExpectedHits;
std::map<std::string,std::map<std::string, int>> all_name2hits;
std::map<std::string,std::map<std::string, std::string>> all_name2info;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("significance,sig", bpo::value<double>(&significance)->default_value(0.01),"The critical value for rejecting the H0 hypothesis.")
		("categories,gmt", bpo::value<std::string>(&categories)->required(), "A .gmt file containing the categories to be tested.")
		("scores,sco", bpo::value<std::string>(&scores), "A whitespace seperated file containing identifier and scores.")
		("identifier, id", bpo::value<std::string>(&identifier), "A file containing identifier line by line.")
		("reference, ref", bpo::value<std::string>(&reference)->required(), "A file containing identifier line by line.")
		("adjustment,adj", bpo::value<std::string>(&adjustment)->default_value("no"),"P-value adjustment method for multiple testing.")
		("minimum,min", bpo::value<int>(&minimum)->default_value(0),"Minimum number of genes allowed in categories.")
		("maximum,max", bpo::value<int>(&maximum)->default_value(1000),"Maximum number of genes allowed in categories.")
		("output,out", bpo::value<std::string>(&out), "Output prefix for text files.");

	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),
		           vm);
		bpo::notify(vm);
	}
	catch(bpo::error& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}
	return true;
}

Category parseTestSet(){
	GeneSetReader reader;
	if( scores != "" ){
		ScoringFile<double> sf = reader.readScoringFile<double>(scores);
		return sf.convert("test");
	}else if(identifier != ""){
		return reader.readGeneSet(identifier,"test");
	}
	Category c("test");
	return c;
}

int main(int argc, char* argv[])
{
	if(! parseArguments(argc, argv)){
		return -1;
	}

	if(scores != "" && identifier != ""){
		std::cerr << "ERROR: Please specify only one input file." << "\n";
	}else if (scores == "" && identifier == ""){
		std::cerr << "ERROR: Please specify a input file." << "\n";
	}

	auto test_set = parseTestSet();
	GeneSetReader reader;
	auto reference_set = reader.readGeneSet(reference,"reference");

	std::ifstream cat(categories);
	if(!cat) {
		std::cerr << "ERROR: Could not open " << categories << " for reading." << std::endl;
	}

	for(std::string line; getline(cat, line);) {

		std::vector<std::pair<std::string, double>> results;
		std::map<std::string, std::string> name2reference;
		std::map<std::string, double> name2ExpectedHits;
		std::map<std::string, int> name2hits;
		std::map<std::string, std::string> name2info;

		std::vector<std::string> strs;
		boost::split(strs,line,boost::is_any_of("\t"));

		std::string categoryName = strs[0];
		std::string category = strs[1];

		// Compute the enrichment
		GMTFile input(category);
		while(input) {
			Category c = input.read();
			int hits  = 0;
			for(auto s : test_set){
				if(c.contains(s)){
					++hits;
				}
			}
			std::cout << "INFO: Processing " << c.name() << ". " << hits << std::endl;
			if(minimum <= hits && hits <= maximum){
				name2reference[c.name()] = c.reference();
				name2hits[c.name()] = hits;

				OverRepresentationAnalysis ora;
				auto enr = ora.computePValue(c, reference_set, test_set);

				name2info[c.name()] = std::get<2>(enr);
				name2ExpectedHits[c.name()] = std::get<1>(enr);
				results.emplace_back(c.name(), std::get<0>(enr));
			}
		}

		// Sort p-values
		std::sort(results.begin(), results.end(),
		          [](const std::pair<std::string, double>& a,
		             const std::pair<std::string, double>& b) {
			return a.second < b.second;
		});

		std::vector<std::pair<std::string, double>> adj = pvalue<double>::adjustPValues(results, adjustment);

		all_results[categoryName] = adj;
		all_name2reference[categoryName] = name2reference;
		all_name2ExpectedHits[categoryName] = name2ExpectedHits;
		all_name2hits[categoryName] = name2hits;
		all_name2info[categoryName] = name2info;
	}

	for(auto p : all_results){
		std::ofstream myfile2;
		myfile2.open (out + "." + p.first + ".txt");
		for(size_t i = 0; i < all_results[p.first].size(); ++i) {
			if(all_results[p.first][i].second <= significance) {
				myfile2 << all_results[p.first][i].first << "\t" << all_name2reference[p.first][all_results[p.first][i].first] << "\t" << all_name2ExpectedHits[p.first][all_results[p.first][i].first] << "\t" << all_name2hits[p.first][all_results[p.first][i].first] << "\t" << all_results[p.first][i].second << "\t" << all_name2info[p.first][all_results[p.first][i].first] << "\n";
			}
		}
		myfile2.close();
	}

	return 0;
}
