#include <gtest/gtest.h>

#include <genetrail2/core/PValue.h>

#include<vector>

using namespace GeneTrail;

const double TOLERANCE = 0.00001;

std::initializer_list<std::pair<std::string,double> > il = { std::make_pair ("A",0.05),  std::make_pair ("B",0.01),  std::make_pair ("C",0.07),  std::make_pair ("D",0.03),  std::make_pair ("E",0.1) };
std::vector<std::pair<std::string,double> > pvalues(il);

std::initializer_list<std::pair<std::string,double> > il2 = { std::make_pair ("A",0.01),  std::make_pair ("B",0.01),  std::make_pair ("C",0.01),  std::make_pair ("D",0.05),  std::make_pair ("E",0.05) };
std::vector<std::pair<std::string,double> > pvalues2(il2);

std::initializer_list<std::pair<std::string,double> > il3 = { std::make_pair ("A",2.01),  std::make_pair ("B",3.01),  std::make_pair ("C",5.01),  std::make_pair ("D",6.05),  std::make_pair ("E",3.05) };
std::vector<std::pair<std::string,double> > pvalues3(il3);

std::initializer_list<std::pair<std::string,double> > il4 = { std::make_pair ("A",0.0),  std::make_pair ("B",0.0),  std::make_pair ("C",0.0),  std::make_pair ("D",0.0),  std::make_pair ("E",0.0) };
std::vector<std::pair<std::string,double> > pvalues4(il4);

TEST(PValue, StepUp){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::stepUp(pvalues3));
	EXPECT_NEAR(adj[0].second, 1.0, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 1.0, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 1.0, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 1.0, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 1.0, TOLERANCE);
}

TEST(PValue, StepDown){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::stepDown(pvalues3));
	EXPECT_NEAR(adj[0].second, 1.0, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 1.0, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 1.0, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 1.0, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 1.0, TOLERANCE);
}

TEST(PValue, StepUp2){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::stepUp(pvalues4));
	EXPECT_NEAR(adj[0].second, 0.0, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.0, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.0, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.0, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.0, TOLERANCE);
}

TEST(PValue, StepDown2){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::stepDown(pvalues4));
	EXPECT_NEAR(adj[0].second, 0.0, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.0, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.0, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.0, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.0, TOLERANCE);
}

TEST(PValue, StepUp3){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::stepUp(pvalues));
	EXPECT_NEAR(adj[0].second, 0.01, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.01, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.03, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.03, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.1, TOLERANCE);
}

TEST(PValue, StepDown3){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::stepDown(pvalues));
	EXPECT_NEAR(adj[0].second, 0.05, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.05, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.07, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.07, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.1, TOLERANCE);
}

/*
 * p.adjust(c(0.05,0.01,0.07,0.03,0.1),method="bonferroni")
 */
TEST(PValue, Bonferroni){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::bonferroni(pvalues));
	EXPECT_NEAR(adj[0].second, 0.25, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.05, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.35, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.15, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.5, TOLERANCE);
}

/*
 * p.adjust(c(0.05,0.01,0.07,0.03,0.1),method="fdrBH")
 */
TEST(PValue, FDR){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::benjamini_hochberg(pvalues));
	EXPECT_NEAR(adj[0].second, 0.05000000, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.07500000, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.08333333, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.08750000, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.10000000, TOLERANCE);
}

/*
 * p.adjust(c(0.01,0.01,0.01,0.05,0.05),method="fdrBH")
 */
TEST(PValue, FDR2){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::benjamini_hochberg(pvalues2));
	EXPECT_NEAR(adj[0].second, 0.01666667, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.01666667, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.01666667, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.05000000, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.05000000, TOLERANCE);
}

/*
 * p.adjust(c(0.05,0.01,0.07,0.03,0.1),method="hochberg")
 */
TEST(PValue, Hochberg){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::hochberg(pvalues));
	EXPECT_NEAR(adj[0].second, 0.05, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.1, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.1, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.1, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.1, TOLERANCE);
}

/*
 * p.adjust(c(0.05,0.01,0.07,0.03,0.1),method="holm")
 */
TEST(PValue, Holm){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::holm(pvalues));
	EXPECT_NEAR(adj[0].second, 0.05, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.12, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.15, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.15, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.15, TOLERANCE);
}

/*
 * p.adjust(c(0.05,0.01,0.07,0.03,0.1),method="fdrBY")
 */
TEST(PValue, BenjaminiYekutieli){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::benjamini_yekutieli(pvalues));
	EXPECT_NEAR(adj[0].second, 0.1141667, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.1712500, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.1902778, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.1997917, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.2283333, TOLERANCE);
}

/*
 * p.adjust(c(0.01,0.01,0.01,0.05,0.05),method="fdrBY")
 */
TEST(PValue, BenjaminiYekutieli2){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::benjamini_yekutieli(pvalues2));
	EXPECT_NEAR(adj[0].second, 0.03805556, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.03805556, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.03805556, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.11416667, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.11416667, TOLERANCE);
}

/*
 * sapply(c(1:5),function(i)1-((1-p[i])^5))
 */
TEST(PValue, Sidak){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::sidak(pvalues));
	EXPECT_NEAR(adj[0].second, 0.22621906, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.04900995, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.30431163, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.14126597, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.40951000, TOLERANCE);
}

/*
 * cummin(sapply(c(1:5),function(i)1-((1-sp[i])^(5-i+1))))
 */
TEST(PValue, HolmSidak){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::holm_sidak(pvalues));
	EXPECT_NEAR(adj[0].second, 0.04900995, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.11470719, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.14262500, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.14262500, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.14262500, TOLERANCE);
}

/*
 * cummax(sapply(c(1:5),function(i)1-((1-sp[i])^(5/i))))
 */
TEST(PValue, Finner){
	std::vector<std::pair<std::string,double> > adj(pvalue<double>::finner(pvalues));
	EXPECT_NEAR(adj[0].second, 0.04900995, TOLERANCE);
	EXPECT_NEAR(adj[1].second, 0.07332097, TOLERANCE);
	EXPECT_NEAR(adj[2].second, 0.08193660, TOLERANCE);
	EXPECT_NEAR(adj[3].second, 0.08672055, TOLERANCE);
	EXPECT_NEAR(adj[4].second, 0.10000000, TOLERANCE);
}

TEST(PValue, Fisher){
	double f = pvalue<double>::fisher(pvalues);
	EXPECT_NEAR(f, 3.796847e-04, TOLERANCE);
}

TEST(PValue, Stouffer){
	std::vector<double> weights(5,0.2);
	double f = pvalue<double>::stouffer(pvalues, weights);
	EXPECT_NEAR(f, 0.0, 0.0001);
}
