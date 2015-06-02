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
#ifndef GT2_CORE_DEPENDENT_SHRINKAGE_T_TEST_H
#define GT2_CORE_DEPENDENT_SHRINKAGE_T_TEST_H

#include <vector>
#include <cmath>
#include <tuple>
#include <utility>

#include "config.h"

#include "Statistic.h"
#include "ShrinkageTTest.h"
#include "OneSampleShrinkageTTest.h"

namespace GeneTrail {

    /**
     * Shrinkage T-Test
     */
    template <typename value_type>
	class GT2_EXPORT DependentShrinkageTTest : public ShrinkageTTest<value_type> {

		public:

		template <typename InputIterator1, typename InputIterator2>
		void checkSize(InputIterator1 begin_first, InputIterator1 end_first, InputIterator2 begin_second, InputIterator2 end_second){
			const size_t num_entries1 = std::distance(begin_first, end_first);
			const size_t num_entries2 = std::distance(begin_second, end_second);

			if(num_entries1 != num_entries2) {
				throw std::out_of_range(
				    "The number of entries in the input ranges do not match. "
				    "Range 1: " +
				    boost::lexical_cast<std::string>(num_entries1) +
				    ", Range 2: " +
				    boost::lexical_cast<std::string>(num_entries2));
			}
		}

		/**
		 * Computes the shrinkage t-statistic for each gene in the given range.
		 *
		 * @param begin_first InputIterator
		 * @param end_first InputIterator
		 * @param begin_second InputIterator
		 * @param end_second InputIterator
		 * @return Vector of t-scores
		 */
		template <typename InputIterator1, typename InputIterator2>
		std::vector<value_type> test(InputIterator1 begin_first, InputIterator1 end_first, InputIterator2 begin_second, InputIterator2 end_second){
			checkSize(begin_first, end_first, begin_second, end_second);
			size_t distance = std::distance(begin_first, end_first);
			std::vector<std::vector<value_type>> diffv;
			diffv.reserve(distance);

			auto it_first = begin_first;
			auto it_second = begin_second;
			for(; (it_first != end_first) || (it_second != end_second); ++it_first, ++it_second){
				std::vector<value_type> diff(std::distance(it_first->begin(), it_first->end()));
				std::transform(it_first->begin(), it_first->end(), it_second->begin(), diff.begin(),
			               [](value_type a, value_type b) { return a - b; });
				diffv.push_back(diff);
			}

			OneSampleShrinkageTTest<value_type> t;
			return t.test(diffv.begin(), diffv.end());
		}
	};
}

#endif // GT2_CORE_DEPENDENT_SHRINKAGE_T_TEST_H

