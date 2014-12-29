#include <gtest/gtest.h>

#include <genetrail2/core/HTest.h>
#include <genetrail2/core/FTest.h>
#include <genetrail2/core/IndependentTTest.h>
#include <genetrail2/core/DependentTTest.h>
#include <genetrail2/core/OneSampleWilcoxonSignedRankTest.h>
#include <genetrail2/core/IndependentShrinkageTTest.h>

#include <config.h>

#include <iostream>
#include <initializer_list>
#include <list>

using namespace GeneTrail;

double TOLERANCE = 0.0001;


std::initializer_list<double> ai = {35.5,31.7,31.2,36.6,22.8,28.0,24.6,26.1,34.5,27.7};
std::initializer_list<double> bi = {45.3,36.0,38.6,44.7,31.4,33.5,28.8,35.8,42.9,35.0};

std::list<double> a(ai);
std::list<double> b(bi);

TEST(HTEST, FTest)
{
	FTest<double> f;
	EXPECT_NEAR(HTest::test(f, a.begin(), a.end(), b.begin(), b.end()), 0.7212895, TOLERANCE);
	EXPECT_NEAR(HTest::twoSidedPValue(f, 0.7212895), 0.6343, TOLERANCE);
	EXPECT_NEAR(HTest::lowerTailedPValue(f, 0.7212895), 0.3172, TOLERANCE);
	EXPECT_NEAR(HTest::upperTailedPValue(f, 0.7212895), 0.6828, TOLERANCE);
	EXPECT_NEAR(HTest::confidenceInterval(f, 0.95).first, 0.1791581, TOLERANCE);
	EXPECT_NEAR(HTest::confidenceInterval(f, 0.95).second, 2.9039072, TOLERANCE);
}

TEST(HTEST, IndependentTest)
{
    IndependentTTest<double> f;
    EXPECT_NEAR(HTest::test(f, a.begin(), a.end(), b.begin(), b.end()), -3.1538, TOLERANCE);
    EXPECT_NEAR(HTest::twoSidedPValue(f, -3.1538), 0.005626, TOLERANCE);
    EXPECT_NEAR(HTest::lowerTailedPValue(f, -3.1538), 0.002813, TOLERANCE);
    EXPECT_NEAR(HTest::upperTailedPValue(f, -3.1538), 0.997186, TOLERANCE);
    EXPECT_NEAR(HTest::confidenceInterval(f,0.95).first, -12.222103, TOLERANCE);
    EXPECT_NEAR(HTest::confidenceInterval(f,0.95).second, -2.437897, TOLERANCE);
}

TEST(HTEST, DependentTest)
{
    DependentTTest<double> f;
    EXPECT_NEAR(HTest::test(f, a.begin(), a.end(), b.begin(), b.end()), -11.3706, TOLERANCE);
    EXPECT_NEAR(HTest::twoSidedPValue(f, -11.3706), 1.217e-06, TOLERANCE);
    EXPECT_NEAR(HTest::lowerTailedPValue(f, -11.3706), 6.084e-07, TOLERANCE);
    EXPECT_NEAR(HTest::upperTailedPValue(f, -11.3706), 1, TOLERANCE);
    EXPECT_NEAR(HTest::confidenceInterval(f,0.95).first, -8.788287, TOLERANCE);
    EXPECT_NEAR(HTest::confidenceInterval(f,0.95).second, -5.871713, TOLERANCE);
}

TEST(HTEST, DependentTTest_throws)
{
	std::vector<double> a { 1.0, 2.0 };
	std::vector<double> b { 1.0, 2.0, 3.0};

	DependentTTest<double> f;
	EXPECT_THROW(HTest::test(f, a.begin(), a.end(), b.begin(), b.end()), std::out_of_range);
}

TEST(HTEST, IndependentShrinkageTTest)
{
	IndependentShrinkageTTest<double> t;
	std::list<std::list<double>> aa;
	aa.push_back(a);
	std::list<std::list<double>> bb;
	bb.push_back(b);
	std::vector<double> result = t.test(aa.begin(), aa.end(), bb.begin(), bb.end());
	//This should be equal to the normal t-test.
	EXPECT_NEAR(result[0], -3.1538, TOLERANCE);
}

TEST(HTEST, IndependentShrinkageTTest2)
{
	IndependentShrinkageTTest<double> t;

	std::list<std::list<double>> aa;
	std::list<double> a1({35.5,31.7,31.2,36.6,22.8});
	aa.push_back(a1);
	std::list<double> a2({28.0,24.6,26.1,34.5,27.7});
	aa.push_back(a2);
	std::list<double> a3({31.7,26.1,44.7,35.0,45.2});
	aa.push_back(a3);

	std::list<std::list<double>> bb;
	std::list<double> b1({45.3,36.0,38.6,44.7,31.4});
	bb.push_back(b1);
	std::list<double> b2({33.5,28.8,35.8,42.9,35.0});
	bb.push_back(b2);
	std::list<double> b3({35.8,42.9,27.7,38.6,45.3});
	bb.push_back(b3);

	std::vector<double> result = t.test(aa.begin(), aa.end(), bb.begin(), bb.end());

	//These values are computed by the 'st' R package.
	//As they round the shrinkage factor I had to use TOLERANCE=0.1,
	//because we want our values as exact as possible.
	EXPECT_NEAR(result[0], -2.1324815, TOLERANCE);
	EXPECT_NEAR(result[1], -2.1551803, 0.1);
	EXPECT_NEAR(result[2], -0.3615726, 0.1);
}
