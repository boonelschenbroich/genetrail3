#include <gtest/gtest.h>

#include <genetrail2/core/GeneLabelPermutationTest.h>
#include <genetrail2/core/Statistic.h>

#include <iostream>
#include <initializer_list>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace GeneTrail;

double TOLERANCE = 0.00001;

std::initializer_list<double> ai = {35.5,-31.7,31.2,36.6,-22.8,28.0,-24.6,26.1,-34.5,27.7,45.3,-36.0,38.6,-44.7,31.4,-33.5,28.8,-35.8,42.9,-35.0};
std::vector<double> a(ai);

TEST(Statistic, Min)
{
	GeneLabelPermutationTest<double, std::vector<double>::iterator> g(a.begin(),a.end(),1000);
	EXPECT_NEAR(g.computeLowerTailedPValue(1,-50.0,statistic::min<double,std::vector<double>::iterator>),0.0,TOLERANCE);
	EXPECT_NEAR(g.computeLowerTailedPValue(1,50.0,statistic::min<double,std::vector<double>::iterator>),1.0,TOLERANCE);
}

