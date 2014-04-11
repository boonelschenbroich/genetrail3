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
#ifndef WILCOXON_MANN_WHITNEY_H
#define WILCOXON_MANN_WHITNEY_H

#include "config.h"

#include "Statistic.h"

namespace GeneTrail {

    /**
     * Wilcoxon-Mann-WhitneyTest
     */
    template <typename value_type, typename InputIterator1,typename InputIterator2>
    class GT2_EXPORT WilcoxonMannWhitneyTest {
    public:

        WilcoxonMannWhitneyTest(value_type tol = 1e-4) : tolerance_(tol) {
        }


		value_type computeU(size_t size1, size_t size2, int rank_sum1){
			return (size1 * size2) + (size1 * (size1 + 1)) / 2 - rank_sum1;
		}

		value_type computeZScore(size_t size1, size_t size2, value_type u){
			value_type mean = (size1 + size2) / 2.0;
			value_type sd = std::sqrt((size1 * size2 * (size1 + size2 + 1)) / 12.0);
			return (u - mean) / sd;
		}

		value_type computeZScore(int rank_sum1, size_t size1, int rank_sum2, size_t size2){
			value_type u1 = computeU(size1, size2, rank_sum1);
			value_type u2 = computeU(size2, size1, rank_sum2);
			depleted_ = (u1 < u2);
			value_type u = depleted_ ? u1 : u2;
			return computeZScore(size1, size2, u);
		}

        /**
         * This method implements a variant of the Wilcoxon Mann Whitney Test
         * 
         *
         * @param Iterator
         * @param Iterator
		 * @param Iterator
		 * @param Iterator
         * @return Z-score for the differences between the two groups
         */
        value_type test(const InputIterator1& first_begin, const InputIterator1& first_end, const InputIterator2& second_begin, const InputIterator2& second_end) {
            std::vector< std::pair< value_type, unsigned int> > r;

            int rank_sum1 = 0;
            int rank_sum2 = 0;

            size_t size1 = std::distance(first_begin,first_end);
            size_t size2 = std::distance(second_begin,second_end);

            auto iit = first_begin;
			for (unsigned int i = 0; i < size1; ++i) {
                r.push_back(std::make_pair(*iit, 0));
				++iit;
            }

			iit = second_begin;
            for (unsigned int i = 0; i < size2; ++i) {
                r.push_back(std::make_pair(*iit, 1));
				++iit;
            }

            sort(r.begin(), r.end());

            for (unsigned int i = 0; i < r.size(); ++i) {
                if (r[i].second == 0) {
                    rank_sum1 += i + 1;
                } else {
                    rank_sum2 += i + 1;
                }
            }

			score_ = computeZScore(rank_sum1, size1, rank_sum2, size2);

            return score_;
        }

		boost::math::normal distribution(){
			boost::math::normal dist(0,1);
			return dist;
		}

		std::pair<value_type, value_type> confidenceInterval(const value_type& alpha) {
			value_type conf = boost::math::quantile(boost::math::complement(distribution(), (1 - alpha) / 2.0));
			return std::make_pair(score_ / conf, score_ * conf);
		}

		bool enriched(){
			return !depleted_;
		}

    protected:

        value_type tolerance_;
		value_type score_;
		bool depleted_;
    };
}

#endif // WILCOXON_MANN_WHITNEY_H

