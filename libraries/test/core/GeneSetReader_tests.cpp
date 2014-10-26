#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <utility>

#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>

#include <config.h>

using namespace GeneTrail;

TEST(Parsing, readScoreFileSimple) {
	GeneSetReader<double> parser;
	GeneSet<double> file = parser.readScoringFile(TEST_DATA_PATH("test_scores.txt"));

	auto it = file.begin();
	ASSERT_NE(file.end(), it);
	EXPECT_EQ("90123", it->first); EXPECT_FLOAT_EQ(4.7813, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("01234", it->first); EXPECT_FLOAT_EQ(3.6381, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("12345", it->first); EXPECT_FLOAT_EQ(-6.42397, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("23456", it->first); EXPECT_FLOAT_EQ(8.81721, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("34567", it->first); EXPECT_FLOAT_EQ(-5.91418, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("45678", it->first); EXPECT_FLOAT_EQ(-2.71999, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("56789", it->first); EXPECT_FLOAT_EQ(3.01346, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("67890", it->first); EXPECT_FLOAT_EQ(-1.74416, it->second);
	ASSERT_EQ(file.end(), ++it);
}

TEST(Parsing, readGeneListSimple) {
	GeneSetReader<double> parser;
	auto file = parser.readGeneList(TEST_DATA_PATH("test_genes.txt"));

	auto it = file.begin();
	ASSERT_NE(file.end(), it);
	EXPECT_EQ("dflgjknfjg", it->first); EXPECT_EQ(0.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("dfg213", it->first); EXPECT_EQ(0.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("132fvwr", it->first); EXPECT_EQ(0.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("324fdsvgf", it->first); EXPECT_EQ(0.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("DTGG", it->first); EXPECT_EQ(0.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("gret", it->first); EXPECT_EQ(0.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("gtluh32", it->first); EXPECT_EQ(0.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("§\"$542", it->first); EXPECT_EQ(0.0, it->second);
	EXPECT_EQ(file.end(), ++it);
}

TEST(Parsing, readGeneListWhitespace) {
	GeneSetReader<double> parser;
	auto file = parser.readGeneList(TEST_DATA_PATH("test_genes2.txt"));

	auto it = file.begin();
	ASSERT_NE(file.end(), it);
	EXPECT_EQ("sdfdfdsf", it->first); EXPECT_EQ(0.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("asdad", it->first); EXPECT_EQ(0.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("dsfd", it->first); EXPECT_EQ(0.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("234", it->first); EXPECT_EQ(0.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("123", it->first); EXPECT_EQ(0.0, it->second);
	EXPECT_EQ(file.end(), ++it);
}

TEST(Parsing, readScoreFileWhitespace) {
	GeneSetReader<double> parser;
	auto file = parser.readScoringFile(TEST_DATA_PATH("test_scores2.txt"));

	auto it = file.begin();
	ASSERT_NE(file.end(), it);
	EXPECT_EQ("asadas", it->first);  EXPECT_EQ(2, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("sdfdsf", it->first);  EXPECT_EQ(3.4, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("dfsdfdf", it->first); EXPECT_EQ(-15, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("dsfsd", it->first);   EXPECT_EQ(-7, it->second);
	EXPECT_EQ(file.end(), ++it);
}

TEST(Parsing, readNAFileSimple) {
	GeneSetReader<double> parser;
	auto file = parser.readNAFile(TEST_DATA_PATH("test_scores.na"));

	auto it = file.begin();
	ASSERT_NE(file.end(), it);
	EXPECT_EQ("a", it->first);  EXPECT_FLOAT_EQ(1.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("desf", it->first);  EXPECT_FLOAT_EQ(2.9, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("sgr", it->first);  EXPECT_FLOAT_EQ(-1.5, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("weffe", it->first);  EXPECT_FLOAT_EQ(3212312.1, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("frewf", it->first);  EXPECT_FLOAT_EQ(213, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("dsfr", it->first);  EXPECT_FLOAT_EQ(123, it->second);
	EXPECT_EQ(file.end(), ++it);
}

TEST(Parsing, readNAFileNasty) {
	GeneSetReader<double> parser;
	auto file = parser.readNAFile(TEST_DATA_PATH("test_scores2.na"));

	auto it = file.begin();
	ASSERT_NE(file.end(), it);
	EXPECT_EQ("sdf",       it->first);  EXPECT_FLOAT_EQ(3.2, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("dsdf",      it->first);  EXPECT_FLOAT_EQ(-1.24, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("1",         it->first);  EXPECT_FLOAT_EQ(1.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("dfsf",      it->first);  EXPECT_FLOAT_EQ(132213.2, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("fvju.dsf",  it->first);  EXPECT_FLOAT_EQ(21321.0, it->second);
	ASSERT_NE(file.end(), ++it);
	EXPECT_EQ("()§%$%\"§", it->first);  EXPECT_FLOAT_EQ(234234.0, it->second);
	EXPECT_EQ(file.end(), ++it);
}

TEST(Parsing, sortIncreasingly) {
    GeneSetReader<double> parser;
    GeneSet<double> file = parser.readScoringFile(TEST_DATA_PATH("test_scores.txt"));
    std::vector<std::pair<std::string, double> > scores = file.getIncreasinglySortedScores();
    bool sorted = true;
    for (unsigned int i = 1; i < scores.size(); ++i) {
        if (scores[i-1].second > scores[i].second) {
            sorted = false;
            break;
        }
    }

    EXPECT_TRUE(sorted);
}

TEST(Parsing, sortDecreasingly) {
    GeneSetReader<double> parser;
    GeneSet<double> file = parser.readScoringFile(TEST_DATA_PATH("test_scores.txt"));
    std::vector<std::pair<std::string, double> > scores = file.getDecreasinglySortedScores();
    bool sorted = true;
    for (unsigned int i = 1; i < scores.size(); ++i) {
        if (scores[i-1].second < scores[i].second) {
            sorted = false;
            break;
        }
    }

    EXPECT_TRUE(sorted);
}

TEST(Parsing, sortAbsolute) {
    GeneSetReader<double> parser;
    GeneSet<double> file = parser.readScoringFile(TEST_DATA_PATH("test_scores.txt"));
    std::vector<std::pair<std::string, double> > scores = file.getAbsoluteSortedScores();
    bool sorted = true;
    for (unsigned int i = 1; i < scores.size(); ++i) {
        if (std::abs(scores[i-1].second) < std::abs(scores[i].second)) {
            sorted = false;
            break;
        }
    }

    EXPECT_TRUE(sorted);
}
