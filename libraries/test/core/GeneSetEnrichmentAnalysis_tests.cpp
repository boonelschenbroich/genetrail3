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
	std::vector<std::string> test;
	for(int i=0; i<8; ++i){
		test.push_back(boost::lexical_cast<std::string>(i));
	}
	Category cat("ca");
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
	EXPECT_NEAR(gsea.computeRightPValue(8,4,-16), 0.01428, TOLERANCE);
	EXPECT_NEAR(gsea.computeOneSidedPValueD(8,4,-16), gsea.computeRightPValue(8,4,-16), TOLERANCE);
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

/*TEST(GeneSetEnrichmentAnalysis, bigExampleCompareOldAndNewImplementation) {
	GeneSetReader<double> reader;
	GeneSet<double> scores = reader.readScoringFile("/home/student/tkehl/GT_TEST/bigScores.txt");
	std::vector<std::string> test_set = scores.getIdentifier(scores.getScores());

	GMTFile input2("/home/student/tkehl/GT_TEST/bigCategories.gmt");
	while(input2) {
		Category c = input2.read();
		GeneSetEnrichmentAnalysis<big_float, int64_t> gsea;
		EXPECT_NEAR(gsea.computeOneSidedPValue(c, test_set).convert_to<double>(), gsea.computeOneSidedPValueD(c, test_set).convert_to<double>(), TOLERANCE);
	}
}*/

TEST(GeneSetEnrichmentAnalysis, bigExampleCompareWebserverAndNewImplementation)
{
	/*GeneSetReader<double> reader;
	GeneSet<double> scores = reader.readScoringFile<double>(
	    "/home/student/tkehl/GT_TEST/bigScores.txt");
	std::vector<std::string> test_set =
	    scores.getIdentifier(scores.getDecreasinglySortedScores());*/

	std::vector<std::string> test_set;
	std::ifstream inputS("/home/student/tkehl/GT_TEST/bigList.txt");
	for(std::string line; getline(inputS, line);) {
		boost::trim(line);
		test_set.push_back(line);
	}
	inputS.close();

	std::map<std::string, double> results;
	std::ifstream input("/home/student/tkehl/GT_TEST/bigGTResults.txt");
	for(std::string line; getline(input, line);) {
		std::vector<std::string> sline;
		boost::split(sline, line, boost::is_any_of(" \t"));
		if(sline.size() == 2) {
			boost::trim(sline[0]);
			boost::trim(sline[1]);
			results[sline[0]] = boost::lexical_cast<double>(sline[1]);
		}else{
			std::cout << line << std::endl;
		}
	}

	GMTFile input2("/home/student/tkehl/GT_TEST/bigCategories.gmt");
	while(input2) {
		Category c = input2.read();
		GeneSetEnrichmentAnalysis<big_float, int64_t> gsea;
		EXPECT_NEAR(
		    gsea.computeOneSidedPValue(c, test_set).convert_to<double>(),
		    results[c.name()],
		    TOLERANCE);
	}
}

/*TEST(GeneSetEnrichmentAnalysis, big1) {
	 GeneSetEnrichmentAnalysis<big_float, int64_t> gsea;
	 int n = 20000;
	 int c = 0;
	 for(int l=10; l<=100; ++l){
		 for(int r=0; r + (n-l) < n; r += (n-l)){
			++c;
	 		EXPECT_NEAR(gsea.computeOneSidedPValueD(n, l, r).convert_to<double>(), gsea.computeRightPValue(n, l, r).convert_to<double>(), TOLERANCE);
		 }
		 for(int r=0; r-l > -n; r-=l){
			++c;
		 	EXPECT_NEAR(gsea.computeOneSidedPValueD(n, l, r).convert_to<double>(), gsea.computeRightPValue(n, l, r).convert_to<double>(), TOLERANCE);
		 }
		 std::cout << c << std::endl;
	 }
}*/
