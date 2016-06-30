#include <vector>
#include <iostream>
#include <utility>

#include <gtest/gtest.h>
#include <genetrail2/core/WeightedGeneSetEnrichmentAnalysis.h>

#include <boost/lexical_cast.hpp>

using namespace GeneTrail;

const double TOLERANCE = 0.00001;

TEST(WeightedGeneSetEnrichmentAnalysis, runningSum) {
	std::vector<std::string> names;
	std::vector<double> values;
	for(int i=7; i>0; --i){
		names.emplace_back(boost::lexical_cast<std::string>(i));
		values.emplace_back(i);
	}
	std::vector<size_t> cat;
	cat.emplace_back(0);
	cat.emplace_back(1);
	cat.emplace_back(2);
	cat.emplace_back(6);
	WeightedGeneSetEnrichmentAnalysis<double> gsea;
	EXPECT_NEAR(gsea.computeRunningSum(cat.begin(), cat.end()), 0.4166666666666, TOLERANCE);
}

TEST(WeightedGeneSetEnrichmentAnalysis, runningSum2) {
	std::vector<std::string> names;
	names.emplace_back("A");
	names.emplace_back("B");
	names.emplace_back("C");
	names.emplace_back("D");
	names.emplace_back("E");
	std::vector<double> values{1.51210870, 0.90155071, -0.08425374, -0.42910173, -1.02161460};
	std::vector<size_t> cat;
	cat.emplace_back(1);
	cat.emplace_back(2);
	WeightedGeneSetEnrichmentAnalysis<double> gsea;
	EXPECT_NEAR(gsea.computeRunningSum(cat.begin(), cat.end()), 0.66666666, TOLERANCE);
}

TEST(WeightedGeneSetEnrichmentAnalysis, runningSum3) {
	std::vector<std::string> names;
	names.emplace_back("A");
	names.emplace_back("B");
	names.emplace_back("C");
	names.emplace_back("D");
	names.emplace_back("E");
	names.emplace_back("F");
	names.emplace_back("G");
	names.emplace_back("H");
	names.emplace_back("I");
	names.emplace_back("J");
	std::vector<double> values;
	values.emplace_back(1.88024972);
	values.emplace_back(1.40113722);
	values.emplace_back(1.29404826);
	values.emplace_back(1.04311852);
	values.emplace_back(0.97672805);
	values.emplace_back(0.09449559);
	values.emplace_back(-0.12832204);
	values.emplace_back(-0.19743087);
	values.emplace_back(-0.41267957);
	values.emplace_back(-0.52441915);
	std::vector<size_t> cat;
	cat.emplace_back(1);
	cat.emplace_back(2);
	cat.emplace_back(4);
	cat.emplace_back(7);
	WeightedGeneSetEnrichmentAnalysis<double> gsea;
	EXPECT_NEAR(gsea.computeRunningSum(cat.begin(), cat.end()), 0.6156423, TOLERANCE);
}
