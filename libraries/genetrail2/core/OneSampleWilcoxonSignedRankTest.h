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
#ifndef GT2_CORE_ONE_SAMPLE_WILCOXON_SIGNED_RANK_TEST_H
#define GT2_CORE_ONE_SAMPLE_WILCOXON_SIGNED_RANK_TEST_H


#include <genetrail2/core/macros.h>
#include <genetrail2/core/Statistic.h>

#include <vector>
#include <cmath>

namespace GeneTrail {

    /**
     * Wilcoxon-Mann-WhitneyTest
     */
    template <typename value_type>
    class GT2_EXPORT OneSampleWilcoxonSignedRankTest {
    public:

        OneSampleWilcoxonSignedRankTest(value_type tol = 1e-4, value_type mu = 0.0) : tolerance_(tol), mu_(mu) {
        	enriched_ = false;
		}


		value_type mean(value_type n){
			return (n * (n + 1.0)) / 4.0;
		}

		value_type var(value_type n){
			return std::sqrt((n*(n + 1.0)*(2.0*n+1.0))/24.0);
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
        template<typename InputIterator>
        value_type test(const InputIterator& begin, const InputIterator& end) {

			std::vector<value_type> diff(begin, end);

            for (auto it=diff.begin(); it != diff.end(); ++it) {
                *it -= mu_;
            }

			std::vector< std::pair< value_type, int> > ranks;
			for(unsigned int i=0; i<diff.size(); ++i){
				if(diff[i] >= 0){
					ranks.push_back(std::make_pair(diff[i], 1));
				}else{
					ranks.push_back(std::make_pair(-diff[i], -1));
				}
			}

			sort(diff.begin(), diff.end());

			value_type p_rank_sum = 0.0;
			value_type n_rank_sum = 0.0;


            for (unsigned int i = 0; i < ranks.size(); ++i) {
                if (ranks[i].second == 1) {
                    p_rank_sum += i+1;
                } else {
                    n_rank_sum += i+1;
                }
            }

			value_type rank_sum = 0.0;
			if(p_rank_sum < n_rank_sum){
				rank_sum = p_rank_sum;
			}else{
				enriched_ = true;
				rank_sum = n_rank_sum;
			}

			score_ = (rank_sum - mean(diff.size())) / var(diff.size());

            return score_;
        }

		bool enriched(){
			return enriched_;
		}

		boost::math::normal distribution(){
			boost::math::normal dist(0,1);
			return dist;
		}

		std::pair<value_type, value_type> confidenceInterval(const value_type& alpha) {
			value_type conf = boost::math::quantile(boost::math::complement(distribution(), (1 - alpha) / 2.0));
			return std::make_pair(score_ / conf, score_ * conf);
		}

    protected:
		bool enriched_;
        value_type tolerance_;
		value_type score_;
		value_type mu_;
    };
}

#endif // GT2_CORE_ONE_SAMPLE_WILCOXON_SIGNED_RANK_TEST_H

