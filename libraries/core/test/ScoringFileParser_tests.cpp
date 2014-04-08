#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <utility>

#include "../src/ScoringFile.h"
#include "../src/GeneSetReader.h"

#include <config.h>

using namespace GeneTrail;

TEST(Parsing, sortIncreasingly) {
    GeneSetReader parser;
    ScoringFile<double> file = parser.readScoringFile<double>(TEST_DATA_PATH("test_scores.txt"));
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
    GeneSetReader parser;
    ScoringFile<double> file =parser.readScoringFile<double>(TEST_DATA_PATH("test_scores.txt"));
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
    GeneSetReader parser;
    ScoringFile<double> file = parser.readScoringFile<double>(TEST_DATA_PATH("test_scores.txt"));
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
