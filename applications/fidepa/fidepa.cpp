#include "../../libraries/core/src/BoostGraphParser.h"
#include "../../libraries/core/src/BoostGraphProcessor.h"
#include "../../libraries/core/src/Pathfinder.h"
#include "../../libraries/core/src/CommandLineParser.h"
#include "../../libraries/core/src/FiDePaRunner.h"

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

int main(int argc, char* argv[])
{
	int pathlength = 0;
	std::string kegg = "";
	std::string scores = "";


	CommandLineParser p("\nFiDePa - Finding Deregulated Paths (Keller et al. 2009) \n\nhttp://bioinformatics.oxfordjournals.org/content/25/21/2787.full \n\nUSAGE");
	p.addOption("help,h", "Display this message");
	p.addTypedOption<int>("path_length,l", "Maximal length of the deregulated paths");
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
		bool increasing = down_regulated ? false : true;
		FiDePaRunner f;
		f.computeDeregulatedPaths(kegg, scores, pathlength, increasing, absolute);
	}
	else
	{
		p.printHelp();
	}

	return 0;
}
