/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2016 Tim Kehl <tkehl@bioinf.uni-sb.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public
 * License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef GT2_REGULATION_REGULATOR_ENRICHMENT_ALGORITHMS_H
#define GT2_REGULATION_REGULATOR_ENRICHMENT_ALGORITHMS_H

#include "RegulatorGeneAssociationEnrichmentResult.h"

#include <genetrail2/core/multiprecision.h>
#include <genetrail2/core/HTest.h>
#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/GeneSetEnrichmentAnalysis.h>
#include <genetrail2/core/WilcoxonRankSumTest.h>
#include <genetrail2/core/ExactQuantileTest.h>

#include <cmath>
#include <limits>

#include <boost/math/special_functions/sign.hpp>

namespace GeneTrail
{
class WRSTest
{
  public:
	WRSTest() {}

	RegulatorEnrichmentResult compute_p_value(RegulatorEnrichmentResult& result, size_t )
	{
		result.p_value = HTest::upperTailedPValue(test_, result.score);
		return result;
	}
	
	template <typename Iterator>
	double compute_score(size_t n, Iterator begin, Iterator end)
	{
		return test_.computeZScore(n, begin, end);
	}
	
	double normalize_score(size_t, size_t, double score)
	{
		//Nothing to do
		return score;
	}
	
  private:
	WilcoxonRankSumTest<double> test_;
};

class KSTest
{
  public:
  
  	using big_int = int64_t;
  
	KSTest() {}

	RegulatorEnrichmentResult compute_p_value(RegulatorEnrichmentResult& result, size_t n)
	{
		big_float p_value;
		p_value = test_.computeRightPValue(n, result.hits, result.score);
		result.p_value = p_value.convert_to<double>();
		return result;
	}
	
	template <typename Iterator>
	double compute_score(size_t n, Iterator begin, Iterator end)
	{
		return test_.computeRunningSum(n, begin, end);
	}
	
	double normalize_score(size_t n, size_t hits, double score)
	{
		return test_.computeScore(n, hits, score).convert_to<double>();
	}

  private:
	GeneSetEnrichmentAnalysis<big_float, big_int> test_;
};

class EMTest
{
  public:
	EMTest()
	{}

	RegulatorEnrichmentResult compute_p_value(RegulatorEnrichmentResult& result,
	                                          size_t n)
	{
		if(!initialized) {
			initialized = true;
			test_.initialize(n);
		}

		result.p_value = test_.computePValue(n, result.hits, result.score, n / 2).convert_to<double>();

		return result;
	}

	template <typename Iterator>
	double compute_score(size_t, Iterator begin, Iterator end)
	{
		return statistic::median<double>(begin, end);
	}

	double normalize_score(size_t , size_t , double score)
	{
		return score;
	}

  private:
	bool initialized = false;
	ExactQuantileTest<big_float> test_;
};
}

#endif // GT2_REGULATION_REGULATOR_ENRICHMENT_ALGORITHMS_H
