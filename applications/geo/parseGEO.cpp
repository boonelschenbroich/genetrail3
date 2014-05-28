#include <genetrail2/core/GEO.h>
#include <genetrail2/core/GEOGDSParser.h>
#include <genetrail2/core/GEOGSEParser.h>

#include <iostream>
#include <map>
#include <vector>

#include "common.h"

using namespace GeneTrail;

Params p;

GEOGDSParser gds_;
GEOGSEParser gse_;

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv, p))
	{
		return -1;
	}

	GEOMap geo;
	if(p.gds)
	{
		geo = gds_.readGDSFile(p.geo_dir +  "/" + p.geo);
	}
	else
	{
		geo = gse_.readGSEFile(p.geo_dir +  "/" + p.geo);
	}

	annotateAndWriteGEO(geo,p);

	return 0;
}
