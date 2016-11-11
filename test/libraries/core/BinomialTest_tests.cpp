#include <gtest/gtest.h>

#include <genetrail2/core/BinomialTest.h>
#include <genetrail2/core/multiprecision.h>

using namespace GeneTrail;

const double TOLERANCE = 0.00001;

TEST(BinomialTest, lowerTailedPValue) {
    BinomialTest<big_float> b;
	EXPECT_NEAR(b.computeLowerTailedPValue(1000, 12, 0.01).convert_to<double>(), 0.79251, TOLERANCE);
}

TEST(BinomialTest, upperTailedPValue) {
    BinomialTest<big_float> b;
	EXPECT_NEAR(b.computeUpperTailedPValue(1000, 12, 0.005).convert_to<double>(), 0.00533, TOLERANCE);
}