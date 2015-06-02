/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2015 Tim Kehl tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_ONE_SAMPLE_SHRINKAGE_T_TEST_H
#define GT2_CORE_ONE_SAMPLE_SHRINKAGE_T_TEST_H

#include <vector>
#include <cmath>
#include <tuple>
#include <utility>

#include "macros.h"

#include "Statistic.h"
#include "ShrinkageTTest.h"

namespace GeneTrail {

	/**
	 * Shrinkage T-Test
	 */
	template <typename value_type>
	class GT2_EXPORT OneSampleShrinkageTTest : public ShrinkageTTest<value_type> {
		public:
        /**
		 * Computes the shrinkage t-statistic for each gene in the given range.
		 *
		 * @param begin InputIterator
		 * @param end InputIterator
		 * @return Vector of t-scores
		 */
		template<typename InputIterator>
		std::vector<value_type> test(InputIterator begin, InputIterator end){
			std::vector<value_type> t_scores;
			std::vector<Gene<value_type>> genes;
			value_type var_median;
			std::tie(genes, var_median) = this->computeAllVariances(begin, end);
			this->computePooledVariances(genes, var_median);

            for(unsigned int i=0; i<genes.size(); ++i){
				auto* a = &genes[i];
				t_scores.emplace_back((a->mean) / std::sqrt(a->shrink_var/a->size));
			}

			return t_scores;
		}

        protected:

        value_type tolerance_;
		value_type df_;
		value_type mu_;
		value_type stdErr_;
		value_type score_;

	};
}

#endif // GT2_CORE_ONE_SAMPLE_SHRINKAGE_T_TEST_H
