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
#ifndef GT2_CORE_WILCOXON_RANK_SUM_TEST_H
#define GT2_CORE_WILCOXON_RANK_SUM_TEST_H


#include <genetrail2/core/macros.h>
#include <genetrail2/core/Statistic.h>

namespace GeneTrail {

    /**
     * Wilcoxon Rank Sum Test
     */
    template <typename value_type, typename InputIterator1,typename InputIterator2>
    class GT2_EXPORT WilcoxonRankSumTest {
    public:

        WilcoxonRankSumTest(value_type tol = 1e-4) : tolerance_(tol) {
        }

		value_type computeZScore(value_type rank_sum1, value_type size1, value_type rank_sum2, value_type size2){
			value_type mu = (size1 * (size1 + size2 + 1) / ((value_type)2.0));
			value_type sd = sqrt((size1 * size2 * (size1 + size2 + 1)) / ((value_type)12.0));
			enriched_ = rank_sum1 > mu;
			return (rank_sum1 - mu) / sd;
		}

        /**
         * This method implements a variant of the Wilcoxon Rank Sum Test
         *
         * @param Iterator
         * @param Iterator
		 * @param Iterator
		 * @param Iterator
         * @return Z-score for the differences between the two groups
         */
        value_type test(const InputIterator1& first_begin, const InputIterator1& first_end, const InputIterator2& second_begin, const InputIterator2& second_end) {
            std::vector< std::pair< value_type, unsigned int> > r;

            value_type rank_sum1 = 0;
            value_type rank_sum2 = 0;

            value_type size1 = std::distance(first_begin,first_end);
			value_type size2 = std::distance(second_begin,second_end);

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

            sort(r.begin(), r.end(),
				[](const std::pair< value_type, int>& a, const std::pair< value_type, int>& b) {
					return a.first < b.first;
				}
			);

            for (unsigned int i = 0; i < r.size(); ++i) {
                unsigned int k = 0;
				value_type rank = i + 1;
				// Find ranks with same value
				for(unsigned int j = i+1; j < r.size(); ++j){
					if(r[j].first != r[i].first){
						break;
					}else{
						k += 1;
					}
				}

				// Build median rank
				if(k != 0){
					rank = ((value_type)(rank + rank + k)) / ((value_type) 2.0);
				}

				// Add to rank sums
				for(int l=i; l <= i+k; ++l){
					if (r[l].second == 0) {
                   		rank_sum1 += rank;
                	} else {
                    	rank_sum2 += rank;
                	}
				}
            }

			//std::cout << "R1: " << rank_sum1 << " n1: " << size1 << " R2: " << rank_sum2 << " n2: " << size2 << std::endl;
			//std::cout << "mu1: " << size1 * (size1 + size2 + 1) / 2.0 << std::endl;
			//std::cout << "mu2: " << size2 * (size1 + size2 + 1) / 2.0 << std::endl;
			score_ = computeZScore(rank_sum1, size1, rank_sum2, size2);
			//std::cout << "ZScore: " << score_ << std::endl;
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
			return enriched_;
		}

    protected:

        value_type tolerance_;
		value_type score_;
		bool enriched_;
    };
}

#endif // GT2_CORE_WILCOXON_RANK_SUM_TEST_H

