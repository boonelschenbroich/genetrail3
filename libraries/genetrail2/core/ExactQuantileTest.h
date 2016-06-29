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
#ifndef GENETRAIL2_EXACT_QUANTILE_TEST_H
#define GENETRAIL2_EXACT_QUANTILE_TEST_H

#include "macros.h"
#include "Statistic.h"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

namespace GeneTrail
{

template <typename float_type> class GT2_EXPORT ExactQuantileTest
{

  public:
	ExactQuantileTest() {}

	void initialize(size_t n)
	{
		increasingly_sorted_scores.resize(n);
		decreasingly_sorted_scores.resize(n);
		for(size_t i = 0; i < n; ++i) {
			increasingly_sorted_scores[i] = i;
			decreasingly_sorted_scores[n - i - 1] = i;
		}
	}

	float_type computeSLOW(const size_t& n, const size_t& k,
	                       const size_t& start, const size_t& end)
	{
		float_type result = 0;
		for(size_t j = start; j <= end; ++j) {
			result += boost::math::binomial_coefficient<float_type>(j - 1, k) *
			          boost::math::binomial_coefficient<float_type>(n - j, k);
		}
		return result /
		       boost::math::binomial_coefficient<float_type>(n, 2 * k + 1);
	}

	float_type computeUpperTailedPValueSLOW(const size_t& n, const size_t& k,
	                                        const size_t& i)
	{
		float_type result = computeSLOW(n, k, k + 1, i);
		return result;
	}

	float_type computeLowerTailedPValueSLOW(const size_t& n, const size_t& k,
	                                        const size_t& i)
	{
		float_type result = computeSLOW(n, k, i, n - k);
		return result;
	}

	float_type compute(const float_type& n, const float_type& k,
	                   const float_type& start, const float_type& end,
	                   float_type& binomial_coefficient_lower,
	                   float_type& binomial_coefficient_upper)
	{
		float_type result =
		    binomial_coefficient_lower * binomial_coefficient_upper;
		for(float_type j = start; j <= end; ++j) {
			binomial_coefficient_lower *= (j - 1) / (j - 1 - k);
			binomial_coefficient_upper *= (n - j + 1 - k) / (n - j + 1);
			result += binomial_coefficient_lower * binomial_coefficient_upper;
		}
		return result;
	}

	float_type compute(const size_t& n, const size_t& k, const size_t& start,
	                   const size_t& end)
	{
		float_type binomial_coefficient_lower =
		    boost::math::binomial_coefficient<float_type>(start - 1, k);
		float_type binomial_coefficient_upper =
		    boost::math::binomial_coefficient<float_type>(n - start, k);
		return compute(boost::numeric_cast<float_type>(n),
		               boost::numeric_cast<float_type>(k),
		               boost::numeric_cast<float_type>(start + 1),
		               boost::numeric_cast<float_type>(end),
		               binomial_coefficient_lower, binomial_coefficient_upper) /
		       boost::math::binomial_coefficient<float_type>(n, 2 * k + 1);
	}

	float_type computeUpperTailedPValueOdd(const size_t& n, const size_t& k,
	                                       const size_t& i)
	{
		return compute(n, k, k + 1, i);
	}

	float_type computeLowerTailedPValueOdd(const size_t& n, const size_t& k,
	                                       const size_t& i)
	{
		return compute(n, k, i, n - k);
	}

	std::tuple<std::vector<float_type>, std::vector<float_type>>
	initializeBinomialCoefficientsCummulative(const size_t& n, const size_t& k)
	{
		std::vector<float_type> binoms(n, boost::numeric_cast<float_type>(0));
		std::vector<float_type> cumm(n, boost::numeric_cast<float_type>(0));
		float_type fk = boost::numeric_cast<float_type>(k);
		float_type fi = boost::numeric_cast<float_type>(k + 1.0);
		binoms[k] = 1;
		cumm[k] = 1;
		for(size_t i = k + 1; i < n; ++i) {
			binoms[i] = binoms[i - 1] * (fi / (fi - fk));
			cumm[i] = binoms[i] + cumm[i - 1];
			++fi;
		}
		return std::make_tuple(binoms, cumm);
	}

	std::vector<float_type> initializeBinomialCoefficients(const size_t& n,
	                                                       const size_t& k)
	{
		std::vector<float_type> binoms(n, boost::numeric_cast<float_type>(0));
		binoms[k] = 1;
		for(size_t i = k + 1; i < n; ++i) {
			binoms[i] =
			    binoms[i - 1] * (boost::numeric_cast<float_type>(i) /
			                     boost::numeric_cast<float_type>(i - k));
		}
		return binoms;
	}

	template <typename Comperator>
	float_type computeEven(const size_t& n, const size_t& k,
	                       const size_t& quantile,
	                       const std::vector<size_t>& scores, Comperator comp)
	{
		std::vector<float_type> binomsCumm;
		std::vector<float_type> binoms;
		std::tie(binoms, binomsCumm) =
		    initializeBinomialCoefficientsCummulative(n, k - 1);
		float_type result = 0;
		for(size_t j = k + 1; j < n - k; ++j) {
			size_t i = k;
			for(; i < j; ++i) {
				if(comp(boost::numeric_cast<float_type>(scores[j] + scores[i]) /
				            2.0,
				        quantile)) {
					break;
				}
			}
			result += binomsCumm[i - 1] * binoms[n - j];
		}
		return result / boost::math::binomial_coefficient<float_type>(n, 2 * k);
	}

	float_type computeLowerTailedPValueEven(const size_t& n, const size_t& k,
	                                        const size_t& quantile)
	{
		return computeEven(n, k, quantile, increasingly_sorted_scores,
		                   [](float_type a, float_type b) { return a > b; });
	}

	float_type computeUpperTailedPValueEven(const size_t& n, const size_t& k,
	                                        const size_t& quantile)
	{
		return computeEven(n, k, quantile, decreasingly_sorted_scores,
		                   [](float_type a, float_type b) { return a < b; });
	}

	float_type computePValue(const size_t& n, const size_t& k,
	                         const size_t& quantile,
	                         const size_t& globalQuantile)
	{
		bool odd = k % 2 != 0;
		float_type p = 0.0;
		if(quantile < globalQuantile) {
			if(odd) {
				p = computeLowerTailedPValueOdd(n, (k - 1) / 2, quantile - 1);
			} else {
				p = computeLowerTailedPValueEven(n, k / 2, quantile);
			}
		} else {
			if(odd) {
				p = computeUpperTailedPValueOdd(n, k / 2, quantile);
			} else {
				p = computeUpperTailedPValueEven(n, k / 2, quantile);
			}
		}
		return p;
	}

  private:
	std::vector<size_t> increasingly_sorted_scores;
	std::vector<size_t> decreasingly_sorted_scores;
};
}

#endif // GENETRAIL2_EXACT_QUANTILE_TEST_H
