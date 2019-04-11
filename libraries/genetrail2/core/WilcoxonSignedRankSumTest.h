/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2013 Tim Kehl <tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_WILCOXON_SIGNED_RANK_SUM_TEST_H
#define GT2_CORE_WILCOXON_SIGNED_RANK_SUM_TEST_H


#include "macros.h"
#include "Statistic.h"

#include <vector>
#include <cmath>

namespace GeneTrail {

    /**
     * Wilcoxon-SignedRanksumTest
     */
    template <typename value_type>
    class GT2_EXPORT WilcoxonSignedRankSumTest {
    public:

		WilcoxonSignedRankSumTest(value_type tol = 1e-4, value_type mu = 0.0) : tolerance_(tol), mu_(mu) {
		}

		value_type var(value_type n){
			return std::sqrt((n*(n + 1.0)*(2.0*n+1.0))/6);
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
	    std::vector<double> diff(begin, end);
            for (auto it=diff.begin(); it != diff.end(); ++it) {
                *it -= mu_;
            }
	    std::vector< std::pair<double, int> > ranks;
	    for(unsigned int i=0; i<diff.size(); ++i){
		if(diff[i] > 0){
			ranks.push_back(std::make_pair(diff[i], 1));
		}else if(diff[i] < 0){
			ranks.push_back(std::make_pair(-diff[i], -1));
		}
	    }
	  
	    sort(ranks.begin(), ranks.end(),
		[](const std::pair< double, int>& a, const std::pair< double, int>& b) {
			return b.first > a.first;
		});

	    double rank_sum = 0.0;
	    if(ranks.size() != 0){
	      for (unsigned int i = 0; i < ranks.size(); ++i) {
		size_t limit = 0;
		for(unsigned int x = i+1 ; x < ranks.size(); ++x){
		  if(ranks[x].first != ranks[i].first){
		    limit = x;
		    break;
		  }
		}
		if(limit == 0){
		  //std::cout << ranks[i].first <<"  " <<ranks[i].second << "  " << i+1 << std::endl;
		  if (ranks[i].second == 1) {
		      rank_sum += i+1;
		  } else {
		      rank_sum -= i+1;
		  }
		}else{
		  double avg = 0;
		  for(size_t s = i+1; s < limit+1; ++s){
		    avg += s;
		  }
		  avg /= (limit-i);
		  for(size_t u = i; u < limit ; ++u){
		  // std::cout << ranks[u].first <<"  " <<ranks[u].second << "  " <<avg << std::endl;
		    if (ranks[u].second == 1) {
			rank_sum += avg;
		    } else {
			rank_sum -= avg;
		    }
		  }
		  i = limit-1;
		}
	      }
	      score_ = rank_sum / var(ranks.size());
	    }
	    return(score_);
	    
	}

    protected:
	value_type tolerance_;
	value_type score_ = 0.0;
	value_type mu_;
    };
}

#endif 

