#include <gtest/gtest.h>

#include "../src/HTest.h"
#include "../src/FTest.h"
#include "../src/IndependentTTest.h"
#include "../src/DependentTTest.h"
#include "../src/WilcoxonMannWhitneyTest.h"

#include <config.h>

#include <iostream>
#include <initializer_list>
#include <list>

using namespace GeneTrail;

double TOLERANCE = 0.001;


std::initializer_list<double> ai = {35.5,31.7,31.2,36.6,22.8,28.0,24.6,26.1,34.5,27.7};
std::initializer_list<double> bi = {45.3,36.0,38.6,44.7,31.4,33.5,28.8,35.8,42.9,35.0};

std::list<double> a(ai);
std::list<double> b(bi);

TEST(HTEST, FTest)
{
	FTest<double,std::list<double>::iterator,std::list<double>::iterator> f;
	EXPECT_NEAR(HTest::test(f, a.begin(), a.end(), b.begin(), b.end()), 0.7212895, TOLERANCE);
	EXPECT_NEAR(HTest::twoSidedPValue(f, 0.7212895), 0.6343, TOLERANCE);
	EXPECT_NEAR(HTest::lowerTailedPValue(f, 0.7212895), 0.3172, TOLERANCE);
	EXPECT_NEAR(HTest::upperTailedPValue(f, 0.7212895), 0.6828, TOLERANCE);
	EXPECT_NEAR(HTest::confidenceInterval(f, 0.95).first, 0.1791581, TOLERANCE);
	EXPECT_NEAR(HTest::confidenceInterval(f, 0.95).second, 2.9039072, TOLERANCE);
}

TEST(HTEST, IndependentTest)
{
    IndependentTTest<double,std::list<double>::iterator,std::list<double>::iterator> f;
    EXPECT_NEAR(HTest::test(f, a.begin(), a.end(), b.begin(), b.end()), -3.1538, TOLERANCE);
    EXPECT_NEAR(HTest::twoSidedPValue(f, -3.1538), 0.005626, TOLERANCE);
    EXPECT_NEAR(HTest::lowerTailedPValue(f, -3.1538), 0.002813, TOLERANCE);
    EXPECT_NEAR(HTest::upperTailedPValue(f, -3.1538), 0.9972, TOLERANCE);
    EXPECT_NEAR(HTest::confidenceInterval(f,0.95).first, -12.222103, TOLERANCE);
    EXPECT_NEAR(HTest::confidenceInterval(f,0.95).second, -2.437897, TOLERANCE);
}

TEST(HTEST, DependentTest)
{
    DependentTTest<double,std::list<double>::iterator,std::list<double>::iterator> f;
    EXPECT_NEAR(HTest::test(f, a.begin(), a.end(), b.begin(), b.end()), -11.3706, TOLERANCE);
    EXPECT_NEAR(HTest::twoSidedPValue(f, -11.3706), 1.217e-06, TOLERANCE);
    EXPECT_NEAR(HTest::lowerTailedPValue(f, -11.3706), 6.084e-07, TOLERANCE);
    EXPECT_NEAR(HTest::upperTailedPValue(f, -11.3706), 1, TOLERANCE);
    EXPECT_NEAR(HTest::confidenceInterval(f,0.95).first, -8.788287, TOLERANCE);
    EXPECT_NEAR(HTest::confidenceInterval(f,0.95).second, -5.871713, TOLERANCE);
}

TEST(HTEST, WilcoxonMannWhitneyTest)
{
	WilcoxonMannWhitneyTest<double,std::list<double>::iterator,std::list<double>::iterator> f;
	EXPECT_NEAR(HTest::test(f, a.begin(), a.end(), b.begin(), b.end()), (16 - 10)/std::sqrt((10 * 10 * (10 + 10 + 1)) / 12.0), TOLERANCE);
	std::cout << "WARNING: Build a new example for this Test, this sucks" << std::endl;
	//EXPECT_NEAR(HTest::twoSidedPValue(f, (16 - 10)/std::sqrt((10 * 10 * (10 + 10 + 1)) / 12.0)), 0.008931, TOLERANCE);
	//EXPECT_NEAR(HTest::lowerTailedPValue(f, (16 - 10)/std::sqrt((10 * 10 * (10 + 10 + 1)) / 12.0)), 0.004465, TOLERANCE);
	//EXPECT_NEAR(HTest::upperTailedPValue(f, (16 - 10)/std::sqrt((10 * 10 * (10 + 10 + 1)) / 12.0)), 0.9966, TOLERANCE);
	//EXPECT_NEAR(HTest::confidenceInterval(f,0.95).first, , TOLERANCE);
	//EXPECT_NEAR(HTest::confidenceInterval(f,0.95).second, , TOLERANCE);
}
