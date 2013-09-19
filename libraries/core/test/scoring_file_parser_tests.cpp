#include <gtest/gtest.h>

#include <vector>
#include <tuple>
#include <algorithm>
#include <stdlib.h>
#include <iostream>

#include "../src/scoring_file_parser.h"

#include <config.h>

using namespace GeneTrail;

TEST(Parsing, sortDescending)
{
	ScoringFileParser parser;
	parser.readScoringFile(TEST_DATA_PATH("test_scores.txt"));
	std::vector<std::tuple<std::string,double> > scores = parser.getDescendinglySortedScores();
	bool sorted = true;
	double tmp = 0;

	for(unsigned int i=0; i<scores.size(); ++i)
	{
		if(i == 0)
		{
			tmp = std::get<1>(scores[i]);
		}
		else
		{
			if(tmp < std::get<1>(scores[i]))
			{
				sorted=false;
				break;
			}
		}
	}

	EXPECT_TRUE(sorted);
}

TEST(Parsing, sortDecreasing)
{
	ScoringFileParser parser;
	parser.readScoringFile(TEST_DATA_PATH("test_scores.txt"));
	std::vector<std::tuple<std::string,double> > scores = parser.getDecreasinglySortedScores();
	bool sorted = true;
	double tmp = 0;

	for(unsigned int i=0; i<scores.size(); ++i)
	{
		if(i == 0)
		{
			tmp = std::get<1>(scores[i]);
		}
		else
		{
			if(tmp > std::get<1>(scores[i]))
			{
				sorted=false;
				break;
			}
		}
	}

	EXPECT_TRUE(sorted);
}

TEST(Parsing, sortAbsolute)
{
	ScoringFileParser parser;
	parser.readScoringFile(TEST_DATA_PATH("test_scores.txt"));
	parser.sortScoringFileAbsolute();
	std::vector<std::tuple<std::string,double> > scores = parser.getAbsoluteSortedScores();
	bool sorted = true;
	double tmp = 0;

	for(unsigned int i=0; i<scores.size(); ++i)
	{
		if(i == 0)
		{
			tmp = std::get<1>(scores[i]);
		}
		else
		{
			if(abs(tmp) < abs(std::get<1>(scores[i])))
			{
				sorted=false;
				break;
			}
		}
	}

	EXPECT_TRUE(sorted);
}
