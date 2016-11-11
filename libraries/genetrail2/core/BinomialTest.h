/**
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2016 Tim Kehl tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_BINOMIAL_TEST_H
#define GT2_CORE_BINOMIAL_TEST_H

#include "Category.h"
#include "macros.h"
#include "multiprecision.h"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <cmath>
#include <iostream>

namespace GeneTrail
{

/**
 * Hypergeometric test
 */
template <typename float_type> class GT2_EXPORT BinomialTest
{
  public:
	/**
	 * This method implements the Binomial test.
	 *
	 * @param n Number of trials
	 * @param i Number of successes
	 * @param p Probability of success
	 *
	 * @return The probability of the parameters
	 */
	float_type compute(const size_t& n, const size_t& i,
	                   const float_type& p) const
	{
		return boost::math::binomial_coefficient<float_type>(n, i) * pow(p, i) *
		       pow(1.0 - p, n - i);
	}

	/**
	 * This method implements the Binomial test.
	 *
	 * @param n Number of trials
	 * @param i Number of successes
	 * @param p Probability of success
	 *
	 * @return Lower-tailed p-value
	 */
	float_type computeLowerTailedPValue(const size_t& n, const size_t& i,
	                                    const float_type& p) const
	{
		float_type pval = 0.0;
		for(size_t k = 0; k <= i; ++k) {
			pval += compute(n, k, p);
		}
		return pval;
	}

	/**
	 * This method implements the Binomial test.
	 *
	 * @param n Number of trials
	 * @param i Number of successes
	 * @param p Probability of success
	 *
	 * @return Upper-tailed p-value
	 */
	float_type computeUpperTailedPValue(const size_t& n, const size_t& i,
	                                    const float_type& p) const
	{
		float_type pval = 0.0;
		for(size_t k = i; k <= n; ++k) {
			pval += compute(n, k, p);
		}
		return pval;
	}

	/**
	 * This method implements the Binomial test.
	 *
	 * @param reference Category representing the reference set
	 * @param test Category representing the test set
	 * @param category Category representing the analyzed category
	 *
	 * @return P-value
	 */
	float_type computePValue(const Category& reference, const Category& testSet,
	                         const Category& category) const
	{
		return computePValue_(reference, testSet, category,
		                      [](size_t n, float_type p, size_t i) {
			return boost::numeric_cast<float_type>(i) <
			       boost::numeric_cast<float_type>(n) * p;
		});
	}

	/**
	 * This method implements the Binomial test.
	 *
	 * @param reference Category representing the reference set
	 * @param test Category representing the test set
	 * @param category Category representing the analyzed category
	 *
	 * @return P-value
	 */
	float_type computeLowerTailedPValue(const Category& reference,
	                                    const Category& testSet,
	                                    const Category& category) const
	{
		return computePValue_(reference, testSet, category,
		                      [](size_t, float_type, size_t) { return true; });
	}

	/**
	 * This method implements the Binomial test.
	 *
	 * @param reference Category representing the reference set
	 * @param test Category representing the test set
	 * @param category Category representing the analyzed category
	 *
	 * @return P-value
	 */
	float_type computeUpperTailedPValue(const Category& reference,
	                                    const Category& testSet,
	                                    const Category& category) const
	{
		return computePValue_(reference, testSet, category,
		                      [](size_t, float_type, size_t) { return false; });
	}

	float_type computeScore(const Category& reference, const Category& testSet,
	                        const Category& category)
	{
		size_t n = category.size();
		size_t i = Category::intersect("null", category, testSet).size();
		size_t k = Category::intersect("null", category, reference).size();
		size_t refsize = reference.size();
		float_type p = boost::numeric_cast<float_type>(k) /
		       boost::numeric_cast<float_type>(refsize); 
		return compute(n, i, p);
	}

	size_t numberOfHits(const Category& testSet, const Category& category) const
	{
		return Category::intersect("null", category, testSet).size();
	}

	double expectedNumberOfHits(const Category& reference,
	                            const Category& category) const
	{
		size_t n = category.size();
		size_t k = Category::intersect("null", category, reference).size();
		size_t refsize = reference.size();
		double e = ((double)k / (double)refsize) * (double) n;
		return e;
	}

  private:
	template <typename Comparator>
	float_type computePValue_(const Category& reference,
	                          const Category& testSet, const Category& category,
	                          Comparator comp) const
	{
		float_type pval;
		size_t n = category.size();
		size_t i = Category::intersect("null", category, testSet).size();
		size_t k = Category::intersect("null", category, reference).size();
		size_t refsize = reference.size();
		float_type p = boost::numeric_cast<float_type>(k) /
		               boost::numeric_cast<float_type>(refsize);
		if(comp(n, p, i)) {
			pval = computeLowerTailedPValue(n, i, p);
		} else {
			pval = computeUpperTailedPValue(n, i, p);
		}
		return pval;
	}
};
}

#endif // GT2_CORE_BINOMIAL_TEST_H
