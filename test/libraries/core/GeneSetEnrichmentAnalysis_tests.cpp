#include <cstdint>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include <gtest/gtest.h>

#include <genetrail2/core/Category.h>
#include <genetrail2/core/GeneSetEnrichmentAnalysis.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/multiprecision.h>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>

using namespace GeneTrail;

const double TOLERANCE = 0.00001;

TEST(GeneSetEnrichmentAnalysis, twoSidedPValue) {
	GeneSetEnrichmentAnalysis<double, int> gsea;
	EXPECT_NEAR(gsea.computeTwoSidedPValue(8,4,12), 0.228571, TOLERANCE);
}

TEST(GeneSetEnrichmentAnalysis, twoSidedPValue2) {
	EntityDatabase db;
	std::vector<std::string> test;
	for(int i=0; i<8; ++i){
		test.push_back(boost::lexical_cast<std::string>(i));
	}
	Category cat(&db);
	cat.insert(boost::lexical_cast<std::string>(0));
	cat.insert(boost::lexical_cast<std::string>(1));
	cat.insert(boost::lexical_cast<std::string>(2));
	cat.insert(boost::lexical_cast<std::string>(6));
	GeneSetEnrichmentAnalysis<double, int> gsea;
	EXPECT_NEAR(gsea.computeRunningSum(cat,test.begin(), test.end()), 12, TOLERANCE);
	EXPECT_NEAR(gsea.computeTwoSidedPValue(cat,test), 0.228571, TOLERANCE);
}

TEST(GeneSetEnrichmentAnalysis, littleExample1) {
	GeneSetEnrichmentAnalysis<double, int> gsea;
	EXPECT_NEAR(gsea.computeRightPValue(8,4,12), 0.11428571, TOLERANCE);
	EXPECT_NEAR(gsea.computeOneSidedPValueD(8,4,12),gsea.computeRightPValue(8,4,12), TOLERANCE);
}

TEST(GeneSetEnrichmentAnalysis, littleExample2) {
	GeneSetEnrichmentAnalysis<double, int> gsea;
	EXPECT_NEAR(gsea.computeLeftPValue(8,4,-16), 0.01428, TOLERANCE);
	EXPECT_NEAR(gsea.computeOneSidedPValueD(8,4,-16), gsea.computeLeftPValue(8,4,-16), TOLERANCE);
}

TEST(GeneSetEnrichmentAnalysis, littleExample3) {
        GeneSetEnrichmentAnalysis<double, int> gsea;
        EXPECT_NEAR(gsea.computeRightPValueD(8,4,-16), gsea.computeRightPValue(8,4,-16), TOLERANCE);
}

TEST(GeneSetEnrichmentAnalysis, smallExampleIndices) {
	/*
	 * For this test we assume a list of ten genes,
	 * where three belong the category we want to examine:
	 *
	 * 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9
	 * 0 | 1 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0
	 *
	 * Here '1' indicates, that the genes was found in the
	 * category.
	 */
	std::vector<size_t> indices { 1, 3, 4 };

	GeneSetEnrichmentAnalysis<double, int> gsea;
	auto max_RSc = gsea.computeRunningSum(10, indices.begin(), indices.end());

	EXPECT_EQ(15, max_RSc);
}

TEST(GeneSetEnrichmentAnalysis, smallExampleIndices2) {
	/*
	 * For this test we assume a list of fifteen genes,
	 * where six belong to the category we want to examine:
	 *
	 * 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14
	 * 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |  1 |  1 |  1 |  0 |  1
	 *
	 * Here '1' indicates, that the genes was found in the
	 * category.
	 */
	std::vector<size_t> indices { 0, 7, 10, 11, 12, 14 };

	GeneSetEnrichmentAnalysis<double, int> gsea;
	auto max_RSc = gsea.computeRunningSum(15, indices.begin(), indices.end());

	EXPECT_EQ(-30, max_RSc);
}
