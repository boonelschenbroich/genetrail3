#include "common.h"

bool parseArguments(int argc, char* argv[], Params& p)
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("geo,g", bpo::value<std::string>(&p.geo)->required(), "Name of the GEO file to be parsed.")
		("geo_dir,l", bpo::value<std::string>(&p.geo_dir)->required(), "Path to the geo directory.")
		("duplicates,d", bpo::value<std::string>(&p.methodToHandleDuplicates)->required(), "Method to handle duplicates.")
		("output,o", bpo::value<std::string>(&p.output_file)->required(), "Name of the output file.")
		("gds,s", bpo::value(&p.gds)->zero_tokens(), "Flag indicating if the input is a GDS file. (default is GSE)");

	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);
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


void annotateAndWriteGEO(GeneTrail::GEOMap& geo, const Params& p)
{
	GeneTrail::GPL_Parser gpl_;
	std::map<std::string, std::string> mappings = gpl_.annotate(p.geo_dir, geo.platform);
	geo.gene2exprs = gpl_.mapAndRemoveDuplicates(geo.gene2exprs, mappings, p.methodToHandleDuplicates);
	gpl_.writeGEOMap(p.output_file, geo);
}
