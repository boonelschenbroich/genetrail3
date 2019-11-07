/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2014-2019 Tim Kehl tkehl@bioinf.uni-sb.de>
 *               2014-2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_HYPERGEOMETRIC_TEST_H
#define GT2_CORE_HYPERGEOMETRIC_TEST_H

#include "macros.h"

#include <boost/math/special_functions/binomial.hpp>
#include <iostream>

namespace GeneTrail
{

	/**
	 * Hypergeometric test
	 */
	template <typename uintt, typename return_type>
	class GT2_EXPORT HypergeometricTest
	{
		public:
		/**
		 * This method implements the Hypergeometric test.
		 *
		 * @param m Number of genes in the reference set (population size)
		 * @param l Number of genes in category (success states in the population)
		 * @param n Number of genes in test set (Number of draws)
		 * @param k Number of genes in the test set that belong to the category (Number of successes)
		 *
		 * @return The hypergeometric probability of the parameters
		 */
		return_type compute(const uintt& m, const uintt& l, const uintt& n,
		                    const uintt& k) const
		{
			return (boost::math::binomial_coefficient<return_type>(l, k) *
			        boost::math::binomial_coefficient<return_type>(m - l, n - k)) /
			       boost::math::binomial_coefficient<return_type>(m, n);
		}

		/**
		 * This method computes the lower-tailed p-value for the given event of
		 * the Hypergeometric test.
		 *
		 * @param m Number of genes in the reference set (population size)
		 * @param l Number of genes in category (success states in the population)
		 * @param n Number of genes in test set (Number of draws)
		 * @param k Number of genes in the test set that belong to the category (Number of successes)
		 *
		 * @return p-value
		 */
		return_type lowerTailedPValue(const uintt& m, const uintt& l,
		                              const uintt& n, const uintt& k) const
		{
			return_type p = 0.0;
			// Make sure we do not compute undefined binomial coefficients
			uintt i = (n > m - l) ? (n + l) - m : 0; // Brackets are here to avoid underrun
			for(; i <= k; ++i) {
				p += boost::math::binomial_coefficient<return_type>(l, i) *
				     boost::math::binomial_coefficient<return_type>(m - l, n - i);
			}
			return p / boost::math::binomial_coefficient<return_type>(m, n);
		}

		/**
		 * This method computes the maximum k needed to get a significant lower-tailed 
		 * p-value for the Hypergeometric test.
		 *
		 * @param m Number of genes in the reference set (population size)
		 * @param l Number of genes in category (success states in the population)
		 * @param n Number of genes in test set (Number of draws)
		 *
		 * @return k
		 */
		uintt findSignificantKLowerTailed(const uintt& m, const uintt& l,
		                              const uintt& n, const return_type& threshold) const
		{
			return_type p = 0.0;
			uintt d = std::min(n, l);
			// Make sure we do not compute undefined binomial coefficients
			uintt i = (n > m - l) ? (n + l) - m : 0; // Brackets are here to avoid underrun
			for(; i <= d; ++i) {
				p += boost::math::binomial_coefficient<return_type>(l, i) *
				     boost::math::binomial_coefficient<return_type>(m - l, n - i);
				if (p > threshold) {
					return i;
				}
			}
			return d;
		}

		/**
		 * This method computes the upper-tailed p-value for the given event of
		 * the Hypergeometric test.
		 *
		 * @param m Number of genes in the reference set (population size)
		 * @param l Number of genes in category (success states in the population)
		 * @param n Number of genes in test set (Number of draws)
		 * @param k Number of genes in the test set that belong to the category (Number of successes)
		 *
		 * @return p-value
		 */
		return_type upperTailedPValue(const uintt& m, const uintt& l,
		                              const uintt& n, const uintt& k) const
		{
			return_type p = 0.0;
			uintt d = std::min(n, l);
			for(uintt i = k; i <= d; ++i) {
				p += boost::math::binomial_coefficient<return_type>(l, i) *
				     boost::math::binomial_coefficient<return_type>(m - l, n - i);
			}
			return p / boost::math::binomial_coefficient<return_type>(m, n);
		}

		/**
		 * This method computes the minimum k needed to get a significant upper-tailed 
		 * p-value for the Hypergeometric test.
		 *
		 * @param m Number of genes in the reference set (population size)
		 * @param l Number of genes in category (success states in the population)
		 * @param n Number of genes in test set (Number of draws)
		 *
		 * @return k
		 */
		uintt findSignificantKUpperTailed(const uintt& m, const uintt& l,
		                              const uintt& n, const return_type& threshold) const
		{
			return_type p = 0.0;
			return_type denom = boost::math::binomial_coefficient<return_type>(m, n);
			uintt d = std::min(n, l);
			for(uintt i = d; i > 1; --i) {
				p += boost::math::binomial_coefficient<return_type>(l, i) *
				     boost::math::binomial_coefficient<return_type>(m - l, n - i);
				if (p/denom > threshold) {
					return i;
				}
			}
			return 0;
		}
	};
}

#endif // GT2_CORE_HYPERGEOMETRIC_TEST_H

