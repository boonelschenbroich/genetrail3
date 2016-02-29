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


#include "macros.h"
#include "Scores.h"
#include "Statistic.h"

#include <boost/math/distributions/normal.hpp>

namespace GeneTrail {

    /**
     * Wilcoxon Rank Sum Test
     */
    template <typename ValueType>
    class GT2_EXPORT WilcoxonRankSumTest {
    public:
		using value_type = ValueType;

        WilcoxonRankSumTest(value_type tol = 1e-4) : tolerance_(tol) {
        }

		value_type computeZScore(value_type rank_sum, value_type size1, value_type size2){
			value_type mu = (size1 * (size1 + size2 + 1) / ((value_type)2.0));
			value_type sd = std::sqrt((size1 * size2 * (size1 + size2 + 1)) / ((value_type)12.0));
			enriched_ = rank_sum > mu;
			return (rank_sum - mu) / sd;
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
        template<typename InputIterator1, typename InputIterator2>
        value_type test(const InputIterator1& first_begin, const InputIterator1& first_end, const InputIterator2& second_begin, const InputIterator2& second_end) {
            std::vector< std::pair< value_type, unsigned int> > r;

            value_type rank_sum1 = 0;
            value_type rank_sum2 = 0;

			size_t size1 = 0;
			for(auto iit = first_begin; iit != first_end; ++iit, ++size1) {
				r.emplace_back(*iit, 0);
            }

			size_t size2 = 0;
			for(auto iit = second_begin; iit != second_end; ++iit, ++size2) {
				r.emplace_back(*iit, 1);
            }

            sort(r.begin(), r.end(),
				[](const std::pair< value_type, int>& a, const std::pair< value_type, int>& b) {
					return a.first < b.first;
				}
			);

			size_t i = 0;
            while(i < r.size()) {
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
				for(unsigned int l=i; l <= i+k; ++l){
					if (r[l].second == 0) {
                   		rank_sum1 += rank;
                	} else {
                    	rank_sum2 += rank;
                	}
				}

				i += k + 1;
            }

			//std::cout << "R1: " << rank_sum1 << " n1: " << size1 << " R2: " << rank_sum2 << " n2: " << size2 << std::endl;
			//std::cout << "mu1: " << size1 * (size1 + size2 + 1) / 2.0 << std::endl;
			//std::cout << "mu2: " << size2 * (size1 + size2 + 1) / 2.0 << std::endl;
			score_ = computeZScore(rank_sum1, size1, size2);
			//std::cout << "ZScore: " << score_ << std::endl;
            return score_;
        }

		/**
		 * This method implements a variant of the Wilcoxon Rank Sum Test
		 *
		 * It is specialised for working with a list of scores sorted by score
		 * and a range of input iterators indicating scores belonging to one of
		 * the groups. The indices must also be provided in to order of the
		 * scores.
		 *
		 * @param scores a scores object that was sorted by score.
		 * @param first start of the range of indices into the passed scores
		 * object.
		 * @param last  end of the range of indices.
		 *
		 * @return Z-score for the differences between the two groups
		 *
		 */
		template <typename InputIterator>
		value_type test_sorted(const Scores& scores, InputIterator first,
		                       const InputIterator& last)
		{
		    auto size1 = std::distance(first, last);
		    auto size2 = scores.size() - size1;
		    assert(size1 <= scores.size());

		    value_type rank_sum1 = 0;
		    value_type rank_sum2 = 0;

		    size_t i = 0;
		    while(i < scores.size()) {
			    size_t k = 0;
			    value_type rank = i + 1;
			    // Find ranks with same value
			    for(size_t j = i + 1; j < scores.size(); ++j) {
				    if(scores[j].score() != scores[i].score()) {
					    break;
				    } else {
					    k += 1;
				    }
			    }

			    // Build median rank
			    if(k != 0) {
				    rank = ((value_type)(rank + rank + k)) / ((value_type)2.0);
			    }

			    // Add to rank sums
			    for(auto l = scores.begin() + i;
			        l != scores.begin() + i + k + 1; ++l) {
				    if(first != last && l->index() == first->index()) {
					    rank_sum1 += rank;
					    ++first;
				    } else {
					    rank_sum2 += rank;
				    }
			    }

			    i += k + 1;
		    }

		    return score_ = computeZScore(rank_sum1, size1, size2);
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

