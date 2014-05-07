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
