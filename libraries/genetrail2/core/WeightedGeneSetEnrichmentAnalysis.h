/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_CORE_WEIGHTED_GENE_SET_ENRICHMENT_ANALYSIS_H
#define GT2_CORE_WEIGHTED_GENE_SET_ENRICHMENT_ANALYSIS_H

#include "macros.h"

#include "Category.h"
#include "EnrichmentResult.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <utility>
#include <tuple>
#include <functional>
#include <cassert>

namespace GeneTrail
{
	template <typename float_type>
	class GT2_EXPORT WeightedGeneSetEnrichmentAnalysis
	{
		constexpr static bool debug = false;

		using Names = std::vector<std::string>;
		using Values = std::vector<double>;
		using Indices = std::vector<size_t>;
		using Pair = std::pair<Names, Values>;

		private:
		const Names& names_;
		Values values_;

		public:
		/**
		 * Constructor
		 */
		WeightedGeneSetEnrichmentAnalysis(const Names& names,
		                                  const Values& values)
		    : names_(names)
		{
			values_.reserve(values.size());
			for(auto v : values) {
				values_.emplace_back(std::abs(v));
			}
		}

		/**
		 * This method computes the interaction between a category and a
		 *testSet.
		 *
		 * @param category Category
		 * @param testSet Container
		 */
		Indices intersection(const Category& category, const Names& names) const
		{
			Indices result;
			result.reserve(std::min(category.size(), names.size()));

			for(int i = 0; i < names.size(); ++i) {
				if(category.contains(names[i])) {
					result.emplace_back(i);
				}
			}

			return result;
		}

		/**
		 * This method computes the sum of absolute values in the container.
		 *
		 * @param S Container
		 */
		float_type sum(Indices::const_iterator it, Indices::const_iterator end) const
		{
			using namespace std;
			float_type result = 0.0;
			for(; it != end; ++it) {
				result += values_[*it];
			}
			return result;
		}

		/**
		 * This method computes the absolute maximum of the two values.
		 *
		 * @param a Floating point number
		 * @param b Floating point number
		 */
		float_type absMax(const float_type a, const float_type b) const
		{
			using namespace std;
			return abs(a) < abs(b) ? b : a;
		}

		/**
		 * This method computes the running sum statistic, based on the given
		 *categories.
		 *
		 * @param category Category for which the RSc should be computed.
		 * @return The RSc for the given categories.
		 */
		float_type computeRunningSum(const Category& category) const
		{
			Indices S = intersection(category, names_);
			return computeRunningSum(S.begin(), S.end());
		}

		/**
		 * This method computes the running sum statistic, based on the given
		 *categories.
		 *
		 * @param category Category for which the RSc should be computed.
		 * @return The RSc for the given categories.
		 */
		float_type computeRunningSum(Indices::const_iterator begin,
		                             Indices::const_iterator end) const
		{
			assert(names_.size() == values_.size());

			const auto category_size = std::distance(begin, end);

			if(category_size == 0) {
				return float_type();
			}

			const float_type NR_inv = 1.0 / sum(begin, end);
			const float_type missv = 1.0 / (names_.size() - category_size);

			float_type maxRS = 0.0;

			float_type RS = -(*begin * missv);
			float_type minRS = RS;
			RS += NR_inv * values_[*begin];
			maxRS = (maxRS > RS) ? maxRS : RS;
			size_t lastIndex = *begin;
			for(auto it = begin + 1; it != end; ++it) {
				RS -= (*it - lastIndex - 1) * missv;
				minRS = (minRS < RS) ? minRS : RS;
				RS += NR_inv * values_[*it];
				maxRS = (maxRS > RS) ? maxRS : RS;

				lastIndex = *it;
			}

			return (maxRS > -minRS) ? maxRS : minRS;
		}
	};
}

#endif // GT2_CORE_WEIGHTED_GENE_SET_ENRICHMENT_ANALYSIS_H

