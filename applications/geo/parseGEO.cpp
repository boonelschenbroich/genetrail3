#include <genetrail2/core/GEO.h>
#include <genetrail2/core/GEOGDSParser.h>
#include <genetrail2/core/GEOGSEParser.h>

#include <iostream>
#include <map>
#include <vector>

#include "common.h"

using namespace GeneTrail;

Params p;

int main(int argc, char* argv[])
{
	GEOGDSParser gds;
	GEOGSEParser gse;

	if(!parseArguments(argc, argv, p))
	{
		return -1;
	}

	std::string path = p.geo_dir +  "/" + p.geo;
	GEOMap geo;
	if(p.gds)
	{
		geo = gds.readGDSFile(path);
	}
	else
	{
		geo = gse.readGSEFile(path);
	}

	annotateAndWriteGEO(geo,p);

	return 0;
}
