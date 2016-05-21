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
#include "Category.h"
#include "Statistic.h"
#include "multiprecision.h"

#include <boost/math/distributions/normal.hpp>

#include <iostream>

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
		
		template<typename Iterator>
		size_t intersectionSize(const Category& category, Iterator begin, const Iterator& end)
		{
			size_t n = 0;
			for(; begin != end; ++begin) {
				if(category.contains(*begin)) {
					++n;
				}
			}
			return n;
		}

		value_type computeZScore_(value_type rank_sum, value_type size1, value_type size2){
			value_type mu = (size1 * (size1 + size2 + 1) / ((value_type)2.0));
			value_type sd = std::sqrt((size1 * size2 * (size1 + size2 + 1)) / ((value_type)12.0));
			enriched_ = rank_sum > mu;
			return (rank_sum - mu) / sd;
		}
		
		value_type computeZScore(std::vector< std::pair< value_type, unsigned int>>& r, size_t size1, size_t size2){
            value_type rank_sum1 = 0;
            value_type rank_sum2 = 0;

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
				for(unsigned int l=i; l <= i+k; ++l){
					if (r[l].second == 0) {
                   		rank_sum1 += rank;
                	}
					//Not needed as we are only intersted in the test set.
					/* else {
                    	rank_sum2 += rank;
                	}*/
				}
            }

			score_ = computeZScore_(rank_sum1, size1, size2);
            return score_;
		}
		
		/**
		 * This method implements a variant of the Wilcoxon Rank Sum Test
		 * List L is assumed to be sorted decreasingly.
		 *
		 * @param category Category for which the Z-score should be computed.
		 * @param begin Iterator (begin) of entries (ids of EntityDatabase) of list L
		 * @param end Iterator (end) of entries (ids of EntityDatabase) of list L
		 * @return The Z-score for the given category.
		 */
		template<typename Iterator>
		value_type computeZScore(const Category& category, const Iterator& begin, const Iterator& end)
		{
			std::vector< std::pair< value_type, unsigned int> > r;

			size_t size1 = 0;
			size_t size2 = 0;
			
			//The counter is to ensure that we do not have ties.
			size_t counter = 0;
			for(auto iit = end; iit != begin; ) {
				--iit;
				if(category.contains(*iit)) {
					++size1;
					r.emplace_back(counter++, 0);
				} else {
					++size2;
					r.emplace_back(counter++, 1);
				}
            }
			
			return computeZScore(r, size1, size2);
		}
		
		/**
		 * This method implements a variant of the Wilcoxon Rank Sum Test.
		 * List L is assumed to be sorted decreasingly and indices start by 0.
		 *
		 * @param n Length of the analysed list L.
		 * @param begin Iterator (begin) of indices of list L that belong to the test set. 
		 * @param end Iterator (end) of indices of list L that belong to the test set. 
		 * @return The Z-score for the given category.
		 */
		template<typename Iterator>
		value_type computeZScore(size_t n, const Iterator& begin, const Iterator& end)
		{
			size_t rank_sum = 0;
			
			for(auto it = begin; it != end; ++it){
				rank_sum += n - *it;
			}
			
			size_t size = std::distance(begin, end);
			return computeZScore_(rank_sum, size, n - size);
		}

        /**
         * This method implements a variant of the Wilcoxon Rank Sum Test
		 * Test set and reference set are assumed to contain scores that 
		 * indicate importance (bigger values are more important).
         *
         * @param first_begin Iterator (begin) of test set
         * @param first_end Iterator (end) of test set
		 * @param second_begin Iterator  (begin) of reference set
		 * @param second_end Iterator  (end) of reference set
         * @return Z-score for the differences between the two groups
         */
        template<typename InputIterator1, typename InputIterator2>
        value_type test(const InputIterator1& first_begin, const InputIterator1& first_end, const InputIterator2& second_begin, const InputIterator2& second_end) {
            std::vector< std::pair< value_type, unsigned int> > r;

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
			
			return computeZScore(r, size1, size2);
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

