#include <cstdint>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include <gtest/gtest.h>

#include "../src/Category.h"
#include "../src/GeneSetEnrichmentAnalysis.h"
#include "../src/GeneSetReader.h"
#include "../src/ScoringFile.h"
#include "../src/GMTFile.h"

#include <boost/lexical_cast.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace GeneTrail;
using namespace boost::multiprecision;

const double TOLERANCE = 0.00001;

TEST(GeneSetEnrichmentAnalysis, twoSidedPValue) {
	GeneSetEnrichmentAnalysis<cpp_dec_float_50, int128_t> gsea;
	EXPECT_NEAR(gsea.computeTwoSidedPValue(8,4,12).convert_to<double>(), 0.228571, TOLERANCE);
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
	GeneSetEnrichmentAnalysis<cpp_dec_float_50, int128_t> gsea;
	EXPECT_NEAR(gsea.computeRunningSum(cat,test), 12, TOLERANCE);
	EXPECT_NEAR(gsea.computeTwoSidedPValue(cat,test).convert_to<double>(), 0.228571, TOLERANCE);
}

TEST(GeneSetEnrichmentAnalysis, littleExample1) {
	GeneSetEnrichmentAnalysis<cpp_dec_float_50, int128_t> gsea;
	EXPECT_NEAR(gsea.computeRightPValue(8,4,12).convert_to<double>(), 0.11428571, TOLERANCE);
	EXPECT_NEAR(gsea.computeOneSidedPValueD(8,4,12).convert_to<double>(),gsea.computeRightPValue(8,4,12).convert_to<double>(), TOLERANCE);
}

TEST(GeneSetEnrichmentAnalysis, littleExample2) {
	GeneSetEnrichmentAnalysis<cpp_dec_float_50, int128_t> gsea;
	EXPECT_NEAR(gsea.computeRightPValue(8,4,-16).convert_to<double>(), 0.01428, TOLERANCE);
	EXPECT_NEAR(gsea.computeOneSidedPValueD(8,4,-16).convert_to<double>(),gsea.computeRightPValue(8,4,-16).convert_to<double>(), TOLERANCE);
}

/*TEST(GeneSetEnrichmentAnalysis, bigExampleCompareOldAndNewImplementation) {
	GeneSetReader reader;
	ScoringFile<double> scores = reader.readScoringFile<double>("/home/student/tkehl/GT_TEST/bigScores.txt");
	std::vector<std::string> test_set = scores.getIdentifier(scores.getScores());

	GMTFile input2("/home/student/tkehl/GT_TEST/bigCategories.gmt");
	while(input2) {
		Category c = input2.read();
		GeneSetEnrichmentAnalysis<cpp_dec_float_50, int64_t> gsea;
		EXPECT_NEAR(gsea.computeOneSidedPValue(c, test_set).convert_to<double>(), gsea.computeOneSidedPValueD(c, test_set).convert_to<double>(), TOLERANCE);
	}
}*/

TEST(GeneSetEnrichmentAnalysis, bigExampleCompareWebserverAndNewImplementation)
{
	/*GeneSetReader reader;
	ScoringFile<double> scores = reader.readScoringFile<double>(
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
		GeneSetEnrichmentAnalysis<cpp_dec_float_50, int64_t> gsea;
		EXPECT_NEAR(
		    gsea.computeOneSidedPValue(c, test_set).convert_to<double>(),
		    results[c.name()],
		    TOLERANCE);
	}
}

/*TEST(GeneSetEnrichmentAnalysis, big1) {
	 GeneSetEnrichmentAnalysis<cpp_dec_float_50, int64_t> gsea;
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
