/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_FISHERS_EXACT_TEST_H
#define GT2_FISHERS_EXACT_TEST_H

#include "macros.h"

#include <boost/math/special_functions/binomial.hpp>

#include <iostream>

namespace GeneTrail {

    /**
     * Fisher's exact test
     */
    template <typename unsigned_integer_type = unsigned int, typename return_type = double >
    class GT2_EXPORT FishersExactTest {
    public:

        /**
         * This method implements the Fisher's exact test.
         *
         * @param m Number of genes in the reference set (population size)
         * @param l Number of genes in category (success states in the population)
         * @param n Number of genes in test set (Number of draws)
         * @param k Number of genes in the test set that belong to the category (Number of successes)
         * @param i
         *
         * @return
         */
        return_type compute(const unsigned_integer_type& m, const unsigned_integer_type& l, const unsigned_integer_type& n, const unsigned_integer_type& k, const unsigned_integer_type& i) {
            return (boost::math::binomial_coefficient<return_type>(n, i) * boost::math::binomial_coefficient<return_type>(m, l + k - i)) / boost::math::binomial_coefficient<return_type>(m + n, l + k);
        }

        /**
         * This method computes the lower-tailed p-value for the given event of the Fisher's exact test.
         *
         * @param m Number of genes in the reference set (population size)
         * @param l Number of genes in category (success states in the population)
         * @param n Number of genes in test set (Number of draws)
         * @param k Number of genes in the test set that belong to the category (Number of successes)
         * @param i
         *
         * @return
         */
        return_type lowerTailedPValue(const unsigned_integer_type& m, const unsigned_integer_type& l, const unsigned_integer_type& n, const unsigned_integer_type& k) {
			return_type p = 0;
            for (unsigned_integer_type i = 0; i <= k; ++i) {
                p += boost::math::binomial_coefficient<return_type>(n, i) * boost::math::binomial_coefficient<return_type>(m, l + k - i);
            }
            return p / boost::math::binomial_coefficient<return_type>(m + n, l + k);
        }

        /**
         * This method computes the upper-tailed p-value for the given event of the Fisher's exact test.
         *
         * @param m Number of genes in the reference set (population size)
         * @param l Number of genes in category (success states in the population)
         * @param n Number of genes in test set (Number of draws)
         * @param k Number of genes in the test set that belong to the category (Number of successes)
         *
         * @return
         */
        return_type upperTailedPValue(const unsigned_integer_type& m, const unsigned_integer_type& l, const unsigned_integer_type& n, const unsigned_integer_type& k) {
			return_type p = 0;
			unsigned_integer_type d = std::min(n, l + k);
			for (unsigned_integer_type i = k; i <= d; ++i) {
            	p += boost::math::binomial_coefficient<return_type>(n, i) * boost::math::binomial_coefficient<return_type>(m, l + k - i);
			}
            return p/boost::math::binomial_coefficient<return_type>(m + n, l + k);
        }
    };
}

#endif // GT2_FISHERS_EXACT_TEST_H

