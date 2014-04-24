#include <gtest/gtest.h>

#include <genetrail2/core/ScoringFile.h>
#include <genetrail2/core/BoostGraphParser.h>
#include <genetrail2/core/BoostGraphProcessor.h>
#include <genetrail2/core/Pathfinder.h>
#include <genetrail2/core/FiDePaRunner.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <string>
#include <cstring>
#include <tuple>
#include <fstream>
#include <set>
#include <vector>
#include <iostream>
#include <ostream>
#include <sstream>
#include <algorithm>
#include <config.h>

using namespace GeneTrail;

std::set<std::string> checkFile(std::string file) {
    std::ifstream input_sif;
    std::string current = "";
    std::set<std::string> dups;

    input_sif.open(file);

    while (std::getline(input_sif, current)) {
        if (current != "") {
            boost::erase_all(current, " ");
            boost::erase_all(current, "\t");
            dups.insert(current);
        }
    }

    input_sif.close();
    return dups;
}

TEST(FiDePaRunner, length_5) {
    //Execute FiDePa
    FiDePaRunner* f = new FiDePaRunner();
    f->computeDeregulatedPaths(TEST_DATA_PATH("FiDePaRunner_test.sif"), TEST_DATA_PATH("FiDePaRunner_test_scores.txt"), 5, true, true);

    //K2
    std::set<std::string> dups = checkFile(TEST_DATA_PATH("FiDePaRunner_test_scores.k2.sif"));
    EXPECT_EQ(dups.size(), 1);
    auto it = dups.find("AppB");
    EXPECT_TRUE(it != dups.end());

    //K3
    dups = checkFile(TEST_DATA_PATH("FiDePaRunner_test_scores.k3.sif"));
    EXPECT_EQ(2, dups.size());
    it = dups.find("BppC");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("CppA");
    EXPECT_TRUE(it != dups.end());

    //K4
    dups = checkFile(TEST_DATA_PATH("FiDePaRunner_test_scores.k4.sif"));
    EXPECT_EQ(3, dups.size());
    it = dups.find("BppC");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("CppA");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("AppD");
    EXPECT_TRUE(it != dups.end());

    //K5
    dups = checkFile(TEST_DATA_PATH("FiDePaRunner_test_scores.k5.sif"));
    EXPECT_EQ(4, dups.size());
    it = dups.find("BppC");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("CppA");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("AppD");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("DppE");
    EXPECT_TRUE(it != dups.end());

    delete f;
}

TEST(FiDePaRunner, pathLengthToHigh) {
    //Execute FiDePa
    FiDePaRunner* f = new FiDePaRunner();
    f->computeDeregulatedPaths(TEST_DATA_PATH("FiDePaRunner_test_pathLengthToHigh.sif"), TEST_DATA_PATH("FiDePaRunner_test_pathLengthToHigh_scores.txt"), 10, true, true);

    //K2
    std::set<std::string> dups = checkFile(TEST_DATA_PATH("FiDePaRunner_test_pathLengthToHigh_scores.k2.sif"));
    EXPECT_EQ(dups.size(), 1);
    auto it = dups.find("AppB");
    EXPECT_TRUE(it != dups.end());

    //K3
    dups = checkFile(TEST_DATA_PATH("FiDePaRunner_test_pathLengthToHigh_scores.k3.sif"));
    EXPECT_EQ(2, dups.size());
    it = dups.find("BppC");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("CppA");
    EXPECT_TRUE(it != dups.end());

    //K4
    dups = checkFile(TEST_DATA_PATH("FiDePaRunner_test_pathLengthToHigh_scores.k4.sif"));
    EXPECT_EQ(3, dups.size());
    it = dups.find("BppC");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("CppA");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("AppD");
    EXPECT_TRUE(it != dups.end());

    //K5
    dups = checkFile(TEST_DATA_PATH("FiDePaRunner_test_pathLengthToHigh_scores.k5.sif"));
    EXPECT_EQ(4, dups.size());
    it = dups.find("BppC");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("CppA");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("AppD");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("DppE");
    EXPECT_TRUE(it != dups.end());

    delete f;
}

TEST(FiDePaRunner, circle) {
    FiDePaRunner* f = new FiDePaRunner();
    f->computeDeregulatedPaths(TEST_DATA_PATH("FiDePaRunner_test_circle.sif"), TEST_DATA_PATH("FiDePaRunner_test_circle_scores.txt"), 10, true, true);
    
    //K2
    std::set<std::string> dups = checkFile(TEST_DATA_PATH("FiDePaRunner_test_circle_scores.k2.sif"));
    EXPECT_EQ(dups.size(), 1);
    auto it = dups.find("AppB");
    EXPECT_TRUE(it != dups.end());

    //K3
    dups = checkFile(TEST_DATA_PATH("FiDePaRunner_test_circle_scores.k3.sif"));
    EXPECT_EQ(2, dups.size());
    it = dups.find("AppB");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("BppC");
    EXPECT_TRUE(it != dups.end());

    //K4
    dups = checkFile(TEST_DATA_PATH("FiDePaRunner_test_circle_scores.k4.sif"));
    EXPECT_EQ(3, dups.size());
    it = dups.find("BppC");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("CppD");
    EXPECT_TRUE(it != dups.end());
    it = dups.find("DppA");
    EXPECT_TRUE(it != dups.end());
    
    delete f;
}
