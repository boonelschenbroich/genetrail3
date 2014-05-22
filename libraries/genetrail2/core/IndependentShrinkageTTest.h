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
#ifndef GT2_CORE_INDEPENDENT_SHRINKAGE_T_TEST_H
#define GT2_CORE_INDEPENDENT_SHRINKAGE_T_TEST_H

#include <vector>
#include <cmath>
#include <tuple>
#include <utility>

#include "config.h"

#include "Statistic.h"
#include "ShrinkageTTest.h"

namespace GeneTrail {

    /**
     * Shrinkage T-Test
     */
    template <typename value_type, typename InputIterator1,typename InputIterator2>
	class GT2_EXPORT IndependentShrinkageTTest : public ShrinkageTTest<value_type, InputIterator1, InputIterator2> {

		public:

		/**
		 * Computes the shrinkage t-statistic for each gene in the given range.
		 *
		 * @param begin_first InputIterator
		 * @param end_first InputIterator
		 * @param begin_second InputIterator
		 * @param end_second InputIterator
		 * @return Vector of t-scores
		 */
		std::vector<value_type> test(InputIterator1 begin_first, InputIterator1 end_first, InputIterator2 begin_second, InputIterator2 end_second){
			std::vector<value_type> t_scores;
			std::vector<Gene<value_type>> genes_first, genes_second;
			value_type var_median_first, var_median_second;
			std::tie(genes_first, var_median_first) = this->computeAllVariances(begin_first, end_first);
			std::tie(genes_second, var_median_second) = this->computeAllVariances(begin_second, end_second);
			this->computePooledVariances(genes_first, var_median_first);
			this->computePooledVariances(genes_second, var_median_second);

			for(int i=0; i<genes_first.size(); ++i){
				auto* a = &genes_first[i];
				auto* b = &genes_second[i];
				t_scores.emplace_back((a->mean - b->mean) / std::sqrt( a->shrink_var/a->size + b->shrink_var/b->size ));
			}
			return t_scores;
		}
	};
}

#endif // GT2_CORE_INDEPENDENT_SHRINKAGE_T_TEST_H

