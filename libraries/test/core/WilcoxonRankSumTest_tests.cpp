#include <gtest/gtest.h>

#include <genetrail2/core/WilcoxonRankSumTest.h>

#include <initializer_list>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace GeneTrail;

double TOLERANCE = 0.00001;

//Indices
constexpr std::initializer_list<size_t> ai = {0, 2, 3};
constexpr std::initializer_list<size_t> bi = {1, 4, 5, 6, 7, 8, 9};

const std::vector<size_t> a(ai);
const std::vector<size_t> b(bi);


//Scores
constexpr std::initializer_list<double> ai2 = {10.0, 8.0, 7.0};
constexpr std::initializer_list<double> bi2 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 9.0};

const std::vector<double> a2(ai2);
const std::vector<double> b2(bi2);

TEST(WilcoxonRankSumTest, Simple1)
{
	WilcoxonRankSumTest<double> test;
    EXPECT_NEAR(test.computeZScore(10, a.begin(), a.end()), test.computeZScore_(25, 3, 7), TOLERANCE);
}

TEST(WilcoxonRankSumTest, Simple2)
{
	WilcoxonRankSumTest<double> test;
    EXPECT_NEAR(test.test(a2.begin(), a2.end(), b2.begin(), b2.end()), test.computeZScore_(25, 3, 7), TOLERANCE);
}