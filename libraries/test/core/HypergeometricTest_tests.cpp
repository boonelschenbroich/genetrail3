#include <gtest/gtest.h>

#include <genetrail2/core/HypergeometricTest.h>
#include <genetrail2/core/multiprecision.h>

using namespace GeneTrail;

const double TOLERANCE = 0.00001;

/*
 * Assume that there are 5 green and 45 red marbles in the urn.
 * Draw 10 marbles without replacement. 
 * What is the probability that exactly 4 of the 10 are green? 
 * What is the probability that exactly 5 of the 10 are green?
 */

/*
 * All results are computed via the R implementation of the hypergeometric test.
 * choose(l,k)*choose(m - l, n - k)/choose(m, n)
 */

TEST(HypergeometricTest, bla) {
    HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.lowerTailedPValue(7442,34,699,1).convert_to<double>(), 0.15753959, TOLERANCE);
}

TEST(HypergeometricTest, blubb) {
    HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.lowerTailedPValue(7442,53,699,2).convert_to<double>(), 0.113443211, TOLERANCE);
}

TEST(HypergeometricTest, compareToR) {
    HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.compute(50,5,10,0).convert_to<double>(), 0.3105628, TOLERANCE);
}

TEST(HypergeometricTest, compareToR1) {
	HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.compute(50,5,10,1).convert_to<double>(), 0.4313372, TOLERANCE);
}

TEST(HypergeometricTest, compareToR2) {
	HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.compute(50,5,10,2).convert_to<double>(), 0.2098397, TOLERANCE);
}

TEST(HypergeometricTest, compareToR3) {
	HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.compute(50,5,10,3).convert_to<double>(), 0.04417678, TOLERANCE);
}

TEST(HypergeometricTest, compareToR4) {
	HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.compute(50,5,10,4).convert_to<double>(), 0.003964583, TOLERANCE);
}

TEST(HypergeometricTest, compareToR5) {
	HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.compute(50,5,10,5).convert_to<double>(), 0.0001189375, TOLERANCE);
}

TEST(HypergeometricTest, lowerTailedPValue) {
	HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.lowerTailedPValue(50,5,10,1).convert_to<double>(), 0.3105628 + 0.4313372, TOLERANCE);
}

TEST(HypergeometricTest, lowerTailedPValue2) {
	HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.lowerTailedPValue(50,5,10,4).convert_to<double>(), 0.3105628 + 0.4313372 + 0.2098397 + 0.04417678 + 0.003964583, TOLERANCE);
}

TEST(HypergeometricTest, lowerTailedPValueBounds) {
	HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.lowerTailedPValue(12,4,10,2).convert_to<double>(), 0.09090909, TOLERANCE);
}

TEST(HypergeometricTest, upperTailedPValue) {
	HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.upperTailedPValue(50,5,10,5).convert_to<double>(), 0.0001189375, TOLERANCE);
}

TEST(HypergeometricTest, upperTailedPValue2) {
 	HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.upperTailedPValue(10,4,6,3).convert_to<double>(), 0.3809524 + 0.07142857, TOLERANCE);
 }

TEST(HypergeometricTest, upperTailedPValue3){
	HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.upperTailedPValue(6,3,3,3).convert_to<double>(), 0.05, TOLERANCE);
}

TEST(HypergeometricTest, upperTailedPValueZero){
	HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.upperTailedPValue(50,5,10,6).convert_to<double>(), 0.0, TOLERANCE);
}

TEST(HypergeometricTest, upperTailedPValueBounds) {
	HypergeometricTest<unsigned int, big_float> h;
	EXPECT_NEAR(h.upperTailedPValue(28, 8, 10, 6).convert_to<double>(), 1.447828e-05 + 0.0006949572 + 0.01033749, TOLERANCE);
}

