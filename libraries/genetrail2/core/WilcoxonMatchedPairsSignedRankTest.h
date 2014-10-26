/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2013 Tim Kehl tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_WILCOXON_MATCHED_PAIRS_SIGNED_RANK_TEST_H
#define GT2_CORE_WILCOXON_MATCHED_PAIRS_SIGNED_RANK_TEST_H


#include <genetrail2/core/macros.h>
#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/OneSampleWilcoxonSignedRankTest.h>

#include <vector>
#include <cmath>
#include <cassert>
#include <iterator>

namespace GeneTrail {

    /**
     * Wilcoxon-Mann-WhitneyTest
     */
    template <typename value_type>
    class GT2_EXPORT WilcoxonMatchedPairsSignedRankTest {

		public:

        WilcoxonMatchedPairsSignedRankTest(value_type tol = 1e-4, value_type mu = 0.0) : tolerance_(tol), mu_(mu), t_(tol, mu) {
		}

        /**
         * This method implements a variant of the Wilcoxon Signed Rank Test
         *
         * @param Iterator
         * @param Iterator
		 * @param Iterator
		 * @param Iterator
         * @return Z-score for the differences between the two groups
         */
        template<typename InputIterator1, typename InputIterator2>
        value_type test(const InputIterator1& first_begin, const InputIterator1& first_end, const InputIterator2& second_begin, const InputIterator2& second_end) {
			assert(std::distance(first_begin,first_end) == std::distance(second_begin,second_end));
			std::vector<value_type> diff(first_begin,first_end);

			int i = 0;
            for (auto it=diff.begin(); it != diff.end(); ++it) {
                *it -= *(second_begin + i);
				++i;
            }

			score_ = t_.test(diff.begin(),diff.end());

            return score_;
        }

		boost::math::normal distribution(){
			return t_.distribution();
		}

		std::pair<value_type, value_type> confidenceInterval(const value_type& alpha) {
			t_.confidenceInterval(alpha);
		}

    	protected:

        value_type tolerance_;
		value_type score_;
		value_type mu_;
		OneSampleWilcoxonSignedRankTest<value_type> t_;
    };
}

#endif // GT2_CORE_WILCOXON_MATCHED_PAIRS_SIGNED_RANK_TEST_H

