#include <gtest/gtest.h>

#include <vector>
#include <tuple>
#include <algorithm>
#include <stdlib.h>
#include <iostream>

#include "../src/ScoringFileParser.h"

#include <config.h>

using namespace GeneTrail;

TEST(Parsing, sortIncreasingly) {
    ScoringFileParser parser;
    parser.readScoringFile(TEST_DATA_PATH("test_scores.txt"));
    std::vector<std::tuple<std::string, double> > scores = parser.getIncreasinglySortedScores();
    bool sorted = true;
    for (unsigned int i = 1; i < scores.size(); ++i) {
        if (std::get<1>(scores[i-1]) > std::get<1>(scores[i])) {
            sorted = false;
            break;
        }
    }

    EXPECT_TRUE(sorted);
}

TEST(Parsing, sortDecreasingly) {
    ScoringFileParser parser;
    parser.readScoringFile(TEST_DATA_PATH("test_scores.txt"));
    std::vector<std::tuple<std::string, double> > scores = parser.getDecreasinglySortedScores();
    bool sorted = true;
    for (unsigned int i = 1; i < scores.size(); ++i) {
        if (std::get<1>(scores[i-1]) < std::get<1>(scores[i])) {
            sorted = false;
            break;
        }
    }

    EXPECT_TRUE(sorted);
}

TEST(Parsing, sortAbsolute) {
    ScoringFileParser parser;
    parser.readScoringFile(TEST_DATA_PATH("test_scores.txt"));
    parser.sortScoringFileAbsolute();
    std::vector<std::tuple<std::string, double> > scores = parser.getAbsoluteSortedScores();
    bool sorted = true;
    for (unsigned int i = 1; i < scores.size(); ++i) {
        if (std::abs(std::get<1>(scores[i-1])) < std::abs(std::get<1>(scores[i]))) {
            sorted = false;
            break;
        }
    }

    EXPECT_TRUE(sorted);
}
