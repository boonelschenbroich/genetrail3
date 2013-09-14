#include <gtest/gtest.h>

#include "../src/graph.h"
#include "../src/scoring_file_parser.h"
#include "../src/graph_parser.h"
#include <config.h>

#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace GeneTrail;

TEST(GraphParser, readAndWrite)
{
	GraphParser parser;
	GraphType g;
	parser.readCytoscapeFile(TEST_DATA_PATH("test_kegg.sif"),g);
	parser.writeCytoscapeFile(TEST_DATA_PATH("test"),g);

	std::ifstream input_sif;
    std::string current = "";
	std::set<std::string> dups;

	//open sif file
    input_sif.open ( TEST_DATA_PATH("test_kegg.sif") );
    while ( std::getline( input_sif, current ) )
    {
    	if ( current != "" )
        {
			boost::erase_all(current, " ");
			boost::erase_all(current, "\t");
			dups.insert(current);
		}
	}

	input_sif.close();

	bool equal = true;

	input_sif.open ( TEST_DATA_PATH("test.sif") );
	while ( std::getline( input_sif, current ) )
    {
        if ( current != "" )
        {
            boost::erase_all(current, " ");
            boost::erase_all(current, "\t");
            std::set<std::string>::iterator it = dups.find(current);
			if(it == dups.end())
			{	
				equal = false;
			}
        }
    }


    EXPECT_TRUE(equal);
}

