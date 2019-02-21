#include <gtest/gtest.h>

#include <genetrail2/core/Entropy.h>

#include <iostream>
#include <initializer_list>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace GeneTrail;

double TOLERANCE = 0.00001;


constexpr std::initializer_list<double> ai = {0.0,0.0,0.0,0.0,0.0,0.0};
constexpr std::initializer_list<double> bi = {1.0,2.0,3.0,4.0,5.0};
constexpr std::initializer_list<double> ci = {-1.0,-2.0,-3.0,-4.0,-5.0};
constexpr std::initializer_list<double> di = {-5.0,-1.0,-2.0,-3.0,-4.0};
constexpr std::initializer_list<double> ei = {
-2.4288314855, -2.2200919553, -1.8547257521, -1.7473247141, -1.7465285986,
-1.7406443388, -1.7235993891, -1.6854530173, -1.6028581134, -1.4917048930,
-1.4767443251, -1.4376808383, -1.3521549811, -1.3487522215, -1.3372159217,
-1.0955078920, -1.0784611231, -1.0423773226, -1.0016687228, -0.9766797163,
-0.8866419096, -0.8005035932, -0.7408119479, -0.7164440467, -0.7150722546,
-0.6700113252, -0.6263447809, -0.6011187943, -0.5899019224, -0.5783118283,
-0.5694084336, -0.5399571922, -0.4901075119, -0.4456101964, -0.4266358195,
-0.4212868011, -0.3733556815, -0.3172085898, -0.2588617204, -0.2139585286,
-0.1941537120, -0.1866998222, -0.1185661060, -0.0990408737, -0.0916497560,
-0.0607272265, -0.0463795206, -0.0443787892, -0.0276632985, -0.0054876268,
 0.0006166116,  0.0056597222,  0.0131847490,  0.0403523690,  0.0522796284,
 0.0745620755,  0.0804426849,  0.0822682864,  0.0892373508,  0.1540986151,
 0.2043294988,  0.2071935391,  0.2085606838,  0.2316877827,  0.2996455612,
 0.3047844542,  0.3272821663,  0.3511956136,  0.3721299719,  0.4091298106,
 0.4210526396,  0.4783000451,  0.5412750386,  0.5622796280,  0.5927251657,
 0.6645578034,  0.6693057720,  0.6791562539,  0.7820694005,  0.8015172524,
 0.8392994135,  0.8526843570,  0.9476941142,  0.9778736730,  1.0108399089,
 1.1700898454,  1.1834549656,  1.1877803545,  1.3179208176,  1.3697384613,
 1.4255614385,  1.5799665186,  1.6269383089,  1.7725670430,  1.9113730397,
 2.0787986791,  2.2991867011,  2.4874321107,  2.9620978930,  3.5846890041
};

std::vector<double> a(ai);
std::vector<double> b(bi);
std::vector<double> c(ci);
std::vector<double> d(di);
std::vector<double> e(ei);

TEST(CummulativeEntropy, Null)
{
	double e = Entropy::cummulative_entropy<double, std::vector<double>::iterator>(a.begin(),a.end());
	EXPECT_NEAR(e, 0.0, TOLERANCE);
}

TEST(CummulativeEntropy, Easy)
{
	double e = Entropy::cummulative_entropy<double, std::vector<double>::iterator>(b.begin(),b.end());
	EXPECT_NEAR(e, 1.173414, TOLERANCE);
}

TEST(CummulativeEntropy, EasyNegative)
{
	double e = Entropy::cummulative_entropy<double, std::vector<double>::iterator>(c.begin(),c.end());
	EXPECT_NEAR(e, 1.173414, TOLERANCE);
}

TEST(CummulativeEntropy, EasyUnordered)
{
	double e = Entropy::cummulative_entropy<double, std::vector<double>::iterator>(d.begin(),d.end());
	EXPECT_NEAR(e, 1.173414, TOLERANCE);
}

TEST(CummulativeEntropy, Random100)
{
	double en = Entropy::cummulative_entropy<double, std::vector<double>::iterator>(e.begin(),e.end());
	EXPECT_NEAR(en, 0.9527841, TOLERANCE);
}

TEST(ConditionalCummulativeEntropy, IndividualElements)
{
	std::vector<std::vector<size_t>> bins;
	for(size_t i=0; i<5; ++i) {
		std::vector<size_t> v;
		v.emplace_back(i);
		bins.emplace_back(v);
	}

	double e = Entropy::conditional_cummulative_entropy_for_sorted_vectors(b, bins);
	EXPECT_NEAR(e, 0.0, TOLERANCE);
}

TEST(ConditionalCummulativeEntropy, OneBin)
{
	std::vector<std::vector<size_t>> bins;
	std::vector<size_t> v;
	for(size_t i=0; i<5; ++i) {
		v.emplace_back(i);
	}
	bins.emplace_back(v);

	double conditional_entropy = Entropy::conditional_cummulative_entropy_for_sorted_vectors(b, bins);
	double entropy = Entropy::cummulative_entropy<double, std::vector<double>::iterator>(b.begin(),b.end());
	EXPECT_NEAR(conditional_entropy, entropy, TOLERANCE);
}

TEST(ConditionalCummulativeEntropy, Random100)
{
	std::vector<size_t> v;
	for(size_t i=0; i<10; ++i) {
		v.emplace_back(i);
	}
	std::sort(e.begin(), e.end());
	double conditional_entropy = Entropy::conditional_cummulative_entropy_for_sorted_vectors(e, v);
	EXPECT_NEAR(conditional_entropy, 0.2468272, TOLERANCE);
}

/*TEST(ConditionalCummulativeEntropyEstimator, Null)
{
	std::vector<std::vector<double>> v;
	v.emplace_back(a);
	v.emplace_back(a);
	std::vector<double> ve = Entropy::conditional_cummulative_entropy_estimator(a,v);
	for(auto e : ve){
		EXPECT_NEAR(e, 0.0, TOLERANCE);
	}
}

TEST(ConditionalCummulativeEntropyEstimator, SameVector)
{
	std::vector<std::vector<double>> v;
	v.emplace_back(b);
	std::vector<double> ve = Entropy::conditional_cummulative_entropy_estimator(b,v);
	EXPECT_NEAR(ve[0], 0.0, TOLERANCE);
	EXPECT_NEAR(ve[1], 0.138629, TOLERANCE);
	EXPECT_NEAR(ve[2], 0.381909, TOLERANCE);
	EXPECT_NEAR(ve[3], 0.727127, TOLERANCE);
	EXPECT_NEAR(ve[4], 1.17341, TOLERANCE);

	v.emplace_back(b);
	v.emplace_back(b);
	ve = Entropy::conditional_cummulative_entropy_estimator(b,v);
	EXPECT_NEAR(ve[0], 0.0, TOLERANCE);
	EXPECT_NEAR(ve[1], 0.138629, TOLERANCE);
	EXPECT_NEAR(ve[2], 0.381909, TOLERANCE);
	EXPECT_NEAR(ve[3], 0.727127, TOLERANCE);
	EXPECT_NEAR(ve[4], 1.17341, TOLERANCE);
}

constexpr std::initializer_list<double> xi = {0.8810036, 4.0287239, 6.4166063, 0.8866507, 3.6051118, 1.9610657, 6.0285489, 6.2974628, 4.4068446, 5.5483233};
constexpr std::initializer_list<double> yi = {4.445753, 3.908179, 9.341760, 7.865155, 1.235848, 7.081790, 8.202711, 4.414463, 8.042939, 6.259775};
constexpr std::initializer_list<double> zi = {5.182783,  3.525994,  7.182913,  1.043968,  7.352742,  9.669824,  5.313436, 8.290781,  3.662965, 10.484494};
std::vector<double> x(xi);
std::vector<double> y(yi);
std::vector<double> z(zi);

TEST(ConditionalCummulativeEntropyEstimator, ConditionedOnTwoVectors10)
{
	std::vector<std::vector<double>> v;
	v.emplace_back(y);
	v.emplace_back(z);
	std::vector<double> ve = Entropy::conditional_cummulative_entropy_estimator(x,v);
	EXPECT_EQ(ve.size(),10);
	EXPECT_NEAR(ve[0], 0.0, TOLERANCE);
	EXPECT_NEAR(ve[1], 0.2486497, TOLERANCE);
	EXPECT_NEAR(ve[2], 0.3610577, TOLERANCE);
	EXPECT_NEAR(ve[3], 0.5792411, TOLERANCE);
	EXPECT_NEAR(ve[4], 0.6764643, TOLERANCE);
	EXPECT_NEAR(ve[5], 1.2131427, TOLERANCE);
	EXPECT_NEAR(ve[6], 1.4193435, TOLERANCE);
	EXPECT_NEAR(ve[7], 1.5987954, TOLERANCE);
	EXPECT_NEAR(ve[8], 1.7039022, TOLERANCE);
	EXPECT_NEAR(ve[9], 1.7563245, TOLERANCE);
}*/