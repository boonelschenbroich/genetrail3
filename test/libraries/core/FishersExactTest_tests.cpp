#include <gtest/gtest.h>

#include <genetrail2/core/FishersExactTest.h>

using namespace GeneTrail;

const double TOLERANCE = 0.00001;

/*
 * All results are computed via the R implementation of the hypergeometric test.
 * choose(n,i)*choose(m, l + k - i)/choose(m + n, l + k)
 */

TEST(FishersExactTest, compareToR6) {
	FishersExactTest<unsigned int, double> fet;
	EXPECT_NEAR(fet.compute(14,6,10,6), 0.2332077, TOLERANCE);
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
