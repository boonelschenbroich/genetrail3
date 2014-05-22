#include <gtest/gtest.h>

#include <genetrail2/core/Statistic.h>

#include <iostream>
#include <initializer_list>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace GeneTrail;

double TOLERANCE = 0.00001;


std::initializer_list<double> ai = {35.5,-31.7,31.2,36.6,-22.8,28.0,-24.6,26.1,-34.5,27.7};
std::initializer_list<double> bi = {45.3,-36.0,38.6,-44.7,31.4,-33.5,28.8,-35.8,42.9,-35.0};

std::vector<double> a(ai);
std::vector<double> b(bi);

TEST(Statistic, Max)
{
	double max1 = statistic::max<double,std::vector<double>::iterator>(a.begin(),a.end());
	EXPECT_NEAR(max1, 36.6, TOLERANCE);
	double max2 = statistic::max<double, std::vector<double>::iterator>(b.begin(),b.end());
	EXPECT_NEAR(max2, 45.3, TOLERANCE);
}

TEST(Statistic, Min)
{
	double min1 = statistic::min<double, std::vector<double>::iterator>(a.begin(),a.end());
	EXPECT_NEAR(min1, -34.5, TOLERANCE);
	double min2 = statistic::min<double, std::vector<double>::iterator>(b.begin(),b.end());
	EXPECT_NEAR(min2, -44.7, TOLERANCE);
}

TEST(Statistic, Abs)
{
	auto tmp = a;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	for(int i=0; i<a.size(); ++i){
		EXPECT_EQ(tmp[i],std::abs(a[i]));
	}
	tmp = b;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	for(int i=0; i<b.size(); ++i){
		EXPECT_EQ(tmp[i],std::abs(b[i]));
	}
}

TEST(Statistic, Mean)
{
	double mean1 = statistic::mean<double, std::vector<double>::iterator>(a.begin(),a.end());
	EXPECT_NEAR(mean1, 7.15, TOLERANCE);
	double mean2 = statistic::mean<double, std::vector<double>::iterator>(b.begin(),b.end());
	EXPECT_NEAR(mean2, 0.2, TOLERANCE);
}

TEST(Statistic, Median)
{
	double median1 = statistic::median<double, std::vector<double>::iterator>(a.begin(),a.end());
	EXPECT_NEAR(median1, 26.1, TOLERANCE);
	double median2 = statistic::median<double, std::vector<double>::iterator>(b.begin(),b.end());
	EXPECT_NEAR(median2, -33.5, TOLERANCE);
}

TEST(Statistic, Var)
{
	double var1 = statistic::var<double, std::vector<double>::iterator>(a.begin(),a.end());
	EXPECT_NEAR(var1, 957.185, TOLERANCE);
	double var2 = statistic::var<double, std::vector<double>::iterator>(b.begin(),b.end());
	EXPECT_NEAR(var2, 1568.93778, TOLERANCE);
}

TEST(Statistic, Sd)
{
	double sd1 = statistic::sd<double, std::vector<double>::iterator>(a.begin(),a.end());
	EXPECT_NEAR(sd1, 30.93841, TOLERANCE);
	double sd2 = statistic::sd<double, std::vector<double>::iterator>(b.begin(),b.end());
	EXPECT_NEAR(sd2, 39.60982, TOLERANCE);
}

TEST(Statistic, Log){
	auto tmp = a;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	statistic::log<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	for(int i=0; i<a.size(); ++i){
		EXPECT_EQ(tmp[i], std::log(std::abs(a[i])));
	}
	tmp = b;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	statistic::log<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	for(int i = 0; i < b.size(); ++i) {
		EXPECT_EQ(tmp[i], std::log(std::abs(b[i])));
	}
}

TEST(Statistic, Log2)
{
	auto tmp = a;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	statistic::log2<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	for(int i = 0; i < a.size(); ++i) {
		EXPECT_EQ(tmp[i], std::log2(std::abs(a[i])));
	}
	tmp = b;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	statistic::log2<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	for(int i = 0; i < b.size(); ++i) {
		EXPECT_EQ(tmp[i], std::log2(std::abs(b[i])));
	}
}

TEST(Statistic, Log10)
{
	auto tmp = a;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),
	                                                      tmp.end());
	statistic::log10<double, std::vector<double>::iterator>(tmp.begin(),
	                                                       tmp.end());
	for(int i = 0; i < a.size(); ++i) {
		EXPECT_EQ(tmp[i], std::log10(std::abs(a[i])));
	}
	tmp = b;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),
	                                                      tmp.end());
	statistic::log10<double, std::vector<double>::iterator>(tmp.begin(),
	                                                       tmp.end());
	for(int i = 0; i < b.size(); ++i) {
		EXPECT_EQ(tmp[i], std::log10(std::abs(b[i])));
	}
}

TEST(Statistic, Sqrt)
{
	auto tmp = a;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),
	                                                      tmp.end());
	statistic::sqrt<double, std::vector<double>::iterator>(tmp.begin(),
	                                                       tmp.end());
	for(int i = 0; i < a.size(); ++i) {
		EXPECT_EQ(tmp[i], std::sqrt(std::abs(a[i])));
	}
	tmp = b;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),
	                                                      tmp.end());
	statistic::sqrt<double, std::vector<double>::iterator>(tmp.begin(),
	                                                       tmp.end());
	for(int i = 0; i < b.size(); ++i) {
		EXPECT_EQ(tmp[i], std::sqrt(std::abs(b[i])));
	}
}

TEST(Statistic, Pow)
{
	auto tmp = a;
	statistic::pow<double, std::vector<double>::iterator>(tmp.begin(), tmp.end(),2);
	for(int i = 0; i < a.size(); ++i) {
		EXPECT_EQ(tmp[i], a[i] * a[i]);
	}
	tmp = b;
	statistic::pow<double, std::vector<double>::iterator>(tmp.begin(), tmp.end(),2);
	for(int i = 0; i < b.size(); ++i) {
		EXPECT_EQ(tmp[i], b[i] * b[i]);
	}
}

TEST(Statistic, Cov)
{
	auto tmp = a;
	auto tmp2 = b;
	auto covar = statistic::cov<double, std::vector<double>::iterator>(tmp.begin(), tmp.end(), tmp2.begin(), tmp2.end());
	EXPECT_NEAR(covar, -382.553333, TOLERANCE);
}

TEST(Statistic, Pearson)
{
	auto tmp = a;
	auto tmp2 = b;
	auto covar = statistic::pearson_correlation<double, std::vector<double>::iterator>(tmp.begin(), tmp.end(), tmp2.begin(), tmp2.end());
	EXPECT_NEAR(covar, -0.31217, TOLERANCE);
}

TEST(Statistic, Spearman)
{
	auto tmp = a;
	auto tmp2 = b;
	auto covar = statistic::spearman_correlation<double, std::vector<double>::iterator>(tmp.begin(), tmp.end(), tmp2.begin(), tmp2.end());
	EXPECT_NEAR(covar, -0.06666667, TOLERANCE);
}
