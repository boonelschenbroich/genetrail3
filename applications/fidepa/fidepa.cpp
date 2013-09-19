#include "../../libraries/core/src/scoring_file_parser.h"
#include "../../libraries/core/src/graph_parser.h"
#include "../../libraries/core/src/graph_processor.h"
#include "../../libraries/core/src/pathfinder.h"
#include "../../libraries/core/src/commandline_parser.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <string>
#include <cstring>
#include <tuple>
#include <fstream>
#include <set>
#include <vector>
#include <iostream>
#include <ostream>
#include <sstream>
#include <algorithm>

using namespace GeneTrail;

std::string convertInt(int number)
{
	std::stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
}

void writeSifFiles(std::vector<std::vector<std::string> > best_paths, std::map<std::string,std::string> regulations, std::string scores)
{
	std::vector<std::string> strs;
	boost::split(strs, scores, boost::is_any_of("."));
	std::cout << strs[0] << std::endl;
	int ki = 2;

	for(std::vector<std::string> path : best_paths)
	{
		std::ofstream myfile;
		std::string file = strs[0] + ".k" + convertInt(ki) + ".sif";
		myfile.open (file.c_str());
		++ki;

		for(unsigned int i=0; i < path.size()-1; ++i)
		{
			auto res = regulations.find(path[i]+path[i+1]);

			if(res != regulations.end() && res->second != "")
			{
				myfile << path[i] << "\t" << res->second << "\t" << path[i+1] << std::endl;
			}
			else
			{
				myfile << path[i] << "\tpp\t" << path[i+1] << std::endl;
			}
		}

		myfile.close();
	}
}

void computeDeregulatedPaths(std::string kegg, std::string scores, int pathlength, bool descending, bool absolute)
{
	GraphType graph;
	GraphParser graph_parser;
	GraphProcessor graph_processor;
	ScoringFileParser scoring_file_parser;

	Pathfinder path_finder;

	graph_parser.readCytoscapeFile<GraphType>(kegg, graph);
	scoring_file_parser.readScoringFile(scores);

	if(absolute)
	{
		scoring_file_parser.sortScoringFileAbsolute();
	}
	else
	{
		scoring_file_parser.sortScoringFile(descending);
	}

	std::set<std::string> vertex_set =  graph_processor.getVertexSet(graph);

	std::vector<std::string> sorted_gene_list = scoring_file_parser.getAllInSet(vertex_set);

	std::set<std::string> vertex_set_with_scores;

	for(std::string s : sorted_gene_list)
	{
		vertex_set_with_scores.insert(s);
	}

	graph_processor.adeptGraph(graph, vertex_set_with_scores);

	std::vector<std::vector<std::string> > best_paths;
	std::map<std::string,std::string> regulations;

	pathlength = (pathlength < (signed)sorted_gene_list.size()-1) ? pathlength : (signed)sorted_gene_list.size()-1;

	path_finder.computeDeregulatedPath(graph, sorted_gene_list, pathlength, best_paths, regulations);

	writeSifFiles(best_paths, regulations, scores);
}

int main(int argc, char* argv[])
{
	int pathlength = 0;
	std::string kegg = "";
	std::string scores = "";


	CommandLineParser p("\nFiDePa - Finding Deregulated Paths (Keller et al. 2009) \n\nhttp://bioinformatics.oxfordjournals.org/content/25/21/2787.full \n\nUSAGE");
	p.addOption("help,h", "Display this message");
	p.addTypedOption<int>("path_length,l", "Maximal length of deregulated the paths");
	p.addTypedOption<std::string>("kegg,k", "[.sif] file containing the network information");
	p.addTypedOption<std::string>("scores,s", "[.txt] file containing gene scores");
	p.addOption("up_regulated,up","Specify to compute up regulated paths");
	p.addOption("down_regulated,down","Specify to compute down regulated paths");
	p.addOption("absolute,abs","Specify to use absolute scores");

	p.parse(argc, argv);

	p.getParameter("path_length",pathlength);
	p.getParameter("kegg",kegg);
	p.getParameter("scores",scores);
	bool up_regulated = p.checkParameter("up_regulated");
	bool down_regulated = p.checkParameter("down_regulated");
	bool absolute = p.checkParameter("absolute");

	if(pathlength > 0 && (up_regulated || down_regulated || absolute) && kegg != "" && scores != "" && !p.checkParameter("help"))
	{
		bool descending = down_regulated ? false : true;
		computeDeregulatedPaths(kegg, scores, pathlength, descending, absolute);
	}
	else
	{
		p.printHelp();
	}

	return 0;
}
