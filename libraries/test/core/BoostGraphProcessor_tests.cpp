#include <gtest/gtest.h>

#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <genetrail2/core/BoostGraph.h>
#include <genetrail2/core/BoostGraphParser.h>
#include <genetrail2/core/BoostGraphProcessor.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/GeneSet.h>

#include <config.h>

using namespace GeneTrail;

TEST(GraphProcessor, vertexSet)
{
	BoostGraphParser parser;
	BoostGraphProcessor processor;
	GraphType g;

	parser.readCytoscapeFile(TEST_DATA_PATH("test_kegg.sif"),g);
	std::set<std::string> vertexSet =  processor.getVertexSet(g);

	EXPECT_EQ(vertexSet.size(),10);

	std::set<std::string>::iterator it;

	for(it=vertexSet.begin(); it!=vertexSet.end(); ++it)
	{
		EXPECT_TRUE(vertexSet.find(*it) != vertexSet.end());
	}
}

TEST(GraphProcessor, graphProcessing)
{
	GeneSetReader<double> reader;
	BoostGraphParser parser;
	BoostGraphProcessor processor;
	GraphType g;

	parser.readCytoscapeFile(TEST_DATA_PATH("graph_processor_test.sif"),g);
	std::set<std::string> vertex_set =  processor.getVertexSet(g);

	GeneSet<double> scores = reader.readScoringFile(TEST_DATA_PATH("graph_processor_test_scores.txt"));
	std::vector<std::string> gene_list = scores.getIdentifier(scores.intersect(vertex_set));

	std::set<std::string> vertex_set_with_scores;

	for(std::string s : gene_list)
	{
		vertex_set_with_scores.insert(s);
	}

	processor.adeptGraph(g, vertex_set_with_scores);

	parser.writeCytoscapeFile(TEST_DATA_PATH("graph_processor_test_output"),g);

	std::ifstream input_sif;
	std::string current = "";
	std::set<std::string> dups;

	//open sif file
	input_sif.open ( TEST_DATA_PATH("graph_processor_adepted.sif") );

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

	input_sif.open ( TEST_DATA_PATH("graph_processor_test_output.sif") );

	while ( std::getline( input_sif, current ) )
	{
		if ( current != "" )
		{
			boost::erase_all(current, " ");
			boost::erase_all(current, "\t");
			auto it = dups.find(current);

			if(it == dups.end())
			{
				equal = false;
			}
		}
	}

	input_sif.close();

	EXPECT_TRUE(equal);
}
