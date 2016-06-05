#include <genetrail2/core/BoostGraphParser.h>
#include <genetrail2/core/BoostGraphProcessor.h>
#include <genetrail2/core/Pathfinder.h>
#include <genetrail2/core/FiDePaRunner.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

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
namespace bpo = boost::program_options;

int main(int argc, char* argv[])
{
	int pathlength = 0;
	std::string kegg = "";
	std::string scores = "";
	bool up_regulated = false, down_regulated = false, absolute = false;

	bpo::variables_map vm;
	bpo::options_description desc("\nFiDePa - Finding Deregulated Paths (Keller et al. 2009)\n\nhttp://bioinformatics.oxfordjournals.org/content/25/21/2787.full\n\nUSAGE");

	desc.add_options()
		("help,h", "Display this message")
		("path_length,l", bpo::value(&pathlength), "Maximal length of the deregulated paths")
		("kegg,k", bpo::value(&kegg), "[.sif] file containing the network information")
		("scores,s", bpo::value(&scores), "[.txt] file containing gene scores")
		("up_regulated,up", bpo::value(&up_regulated)->zero_tokens(), "Specify to compute up regulated paths")
		("down_regulated,down", bpo::value(&down_regulated)->zero_tokens(), "Specify to compute down regulated paths")
		("absolute,abs", bpo::value(&absolute)->zero_tokens(), "Specify to use absolute scores")
	;

	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),
		           vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "ERROR: " << e.what() << "\n";
		desc.print(std::cerr);
		return -1;
	}

	if(!vm["help"].empty()) {
		desc.print(std::cout);
		return 0;
	}

	if(pathlength > 0 && (up_regulated || down_regulated || absolute) && kegg != "" && scores != "")
	{
		bool increasing = down_regulated ? false : true;
		FiDePaRunner f;
		f.computeDeregulatedPaths(kegg, scores, pathlength, increasing, absolute);
	}
	else
	{
		desc.print(std::cerr);
		return -2;
	}

	return 0;
}
