#include <gtest/gtest.h>

#include "../src/FishersExactTest.h"

using namespace GeneTrail;

const double TOLERANCE = 0.00001;

/*
 * All results are computed via the R implementation of the hypergeometric test.
 * choose(n,i)*choose(m, l + k - i)/choose(m + n, l + k)
 */

TEST(FishersExactTest, compareToR) {
    FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.compute(14,6,10,6,3), 0.08884103, TOLERANCE);
}

TEST(FishersExactTest, compareToR1) {
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.compute(14,9,10,3,2), 0.01665769, TOLERANCE);
}

TEST(FishersExactTest, compareToR2) {
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.compute(14,9,10,3,1), 0.001346076, TOLERANCE);
}

TEST(FishersExactTest, compareToR3) {
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.compute(14,9,10,3,0), 0.000033652, TOLERANCE);
}

TEST(FishersExactTest, compareToR4) {
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.compute(14,6,10,6,4), 0.2332077, TOLERANCE);
}

TEST(FishersExactTest, compareToR5) {
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.compute(14,6,10,6,5), 0.3198277, TOLERANCE);
}

TEST(FishersExactTest, compareToR6) {
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.compute(14,6,10,6,6), 0.2332077, TOLERANCE);
}

TEST(FishersExactTest, compareToR7) {
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.compute(14,6,10,6,7), 0.08884103, TOLERANCE);
}

TEST(FishersExactTest, compareToR8) {
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.compute(14,6,10,6,8), 0.01665769, TOLERANCE);
}

TEST(FishersExactTest, compareToR9) {
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.compute(14,6,10,6,9), 0.001346076, TOLERANCE);
}

TEST(FishersExactTest, compareToR10) {
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.compute(14,6,10,6,10), 3.36519e-05, TOLERANCE);
}

TEST(FishersExactTest, lowerTailedPValue) {
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.lowerTailedPValue(14,9,10,1), 0.0005103872 + 0.01020774, TOLERANCE);
}

TEST(FishersExactTest, lowerTailedPValue2) {
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.lowerTailedPValue(14,9,10,3), 0.08884103 + 0.01665769 + 0.001346076 + 0.000033652, TOLERANCE);
}

TEST(FishersExactTest, upperTailedPValue){
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.upperTailedPValue(14,9,10,7), 0.3266478 + 0.1837394 + 0.04666397 + 0.004083098, TOLERANCE);
}
