#include <gtest/gtest.h>

#include <set>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "../src/scoring_file_parser.h"
#include "../src/graph_parser.h"
#include "../src/graph_processor.h"
#include "../src/pathfinder.h"

#include <config.h>

using namespace GeneTrail;

TEST(Pathfinder, length4)
{

	//Execute FiDePa

	std::system("/home/student/tkehl/workspace/FiDePa/build/fidepa -k /home/student/tkehl/workspace/FiDePa/data/pathfinder_test.sif -s /home/student/tkehl/workspace/FiDePa/data/pathfinder_test_scores.txt -l 4 -up");

	//TEST
	std::ifstream input_sif;
	std::string current = "";
	std::set<std::string> dups,dups2,dups3;

	//K2
	input_sif.open ( TEST_DATA_PATH("pathfinder_test_scores.k2.sif") );

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
	EXPECT_EQ(dups.size(),1);
	std::set<std::string>::iterator it = dups.find("AppB");
	EXPECT_TRUE(it != dups.end());

	//K3
	input_sif.open ( TEST_DATA_PATH("pathfinder_test_scores.k3.sif") );

	while ( std::getline( input_sif, current ) )
	{
		if ( current != "" )
		{
			boost::erase_all(current, " ");
			boost::erase_all(current, "\t");
			dups2.insert(current);
		}
	}

	input_sif.close();
	EXPECT_EQ(2,dups2.size());
	it = dups.find("BppC");
	EXPECT_TRUE(it != dups2.end());
	it = dups.find("CppA");
	EXPECT_TRUE(it != dups2.end());

	//K4
	input_sif.open ( TEST_DATA_PATH("pathfinder_test_scores.k4.sif") );

	while ( std::getline( input_sif, current ) )
	{
		if ( current != "" )
		{
			boost::erase_all(current, " ");
			boost::erase_all(current, "\t");
			dups3.insert(current);
		}
	}

	input_sif.close();
	EXPECT_EQ(3,dups3.size());
	it = dups.find("BppC");
	EXPECT_TRUE(it != dups3.end());
	it = dups.find("CppA");
	EXPECT_TRUE(it != dups3.end());
	it = dups.find("AppD");
	EXPECT_TRUE(it != dups3.end());
}
