#include <gtest/gtest.h>

#include <genetrail2/core/PermutationTest.h>
#include <genetrail2/core/Statistic.h>

#include <iostream>
#include <initializer_list>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace GeneTrail;

double TOLERANCE = 0.00001;

TEST(Statistic, Min)
{
	TestResult<double> t1("T1",0.0,1);
	TestResult<double> t2("T2",1.0,1);
	TestResult<double> t3("T3",3.0,1);
	TestResult<double> t4("T4",5.0,1);
	TestResult<double> t5("T5",7.0,1);
	TestResult<double> t6("T6",9.0,1);
	TestResult<double> t7("T7",11.0,1);
	std::initializer_list<TestResult<double>> ai = {t1,t2,t3,t4,t5,t6,t7};
	std::vector<TestResult<double>> a(ai);

	PermutationTest<double> g(a.begin(),a.end(),1000);
	EXPECT_NEAR(1.000000001,1.0,TOLERANCE);
}

