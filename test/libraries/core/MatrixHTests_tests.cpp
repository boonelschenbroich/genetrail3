#include <gtest/gtest.h>

#include <genetrail2/core/MatrixHTest.h>
#include <genetrail2/core/DenseMatrix.h>

#include <config.h>

#include <vector>

using namespace GeneTrail;

double TOLERANCE = 0.0001;

TEST(MatrixHTest, Zscore)
{
	DenseMatrix sample(3, 1);
	sample(0, 0) = 6.47237;
	sample(1, 0) = 6.06996;
	sample(2, 0) = 9.42535;

	sample.setRowName(0, "A");
	sample.setRowName(1, "B");
	sample.setRowName(2, "C");

	sample.setColName(0, "GSM1");

	DenseMatrix reference(3, 3);
	reference(0, 0) = 6.71249;
	reference(0, 1) = 6.22035;
	reference(0, 2) = 5.84363;

	reference(1, 0) = 6.02452;
	reference(1, 1) = 6.40401;
	reference(1, 2) = 6.60834;

	reference(2, 0) = 9.40750;
	reference(2, 1) = 9.26683;
	reference(2, 2) = 8.22949;

	reference.setRowName(0, "A");
	reference.setRowName(1, "B");
	reference.setRowName(2, "C");

	reference.setColName(0, "GSM2");
	reference.setColName(1, "GSM3");
	reference.setColName(2, "GSM4");

	MatrixHTest htest;
	auto scores = htest.test("z-score", sample, reference);
	EXPECT_EQ(scores[0].name(*scores.db()), "A");
	EXPECT_EQ(scores[1].name(*scores.db()), "B");
	EXPECT_EQ(scores[2].name(*scores.db()), "C");

	EXPECT_NEAR(scores[0].score(), 0.4901166, TOLERANCE);
	EXPECT_NEAR(scores[1].score(), -0.9304872, TOLERANCE);
	EXPECT_NEAR(scores[2].score(), 0.7109566, TOLERANCE);
}

