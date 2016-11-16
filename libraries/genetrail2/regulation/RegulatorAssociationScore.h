/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2016 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_REGULATOR_ASSOCIATION_SCORE_H
#define GT2_REGULATOR_ASSOCIATION_SCORE_H

#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/Exception.h>

#include <cassert>
#include <string>
#include <unordered_map>

namespace GeneTrail
{

class GT2_EXPORT PearsonCorrelation
{
  public:
	PearsonCorrelation() {}

	template <typename Iterator1, typename Iterator2>
	typename std::iterator_traits<Iterator1>::value_type
	compute(Iterator1 begin1, Iterator1 end1, Iterator2 begin2,
	        Iterator2 end2) const
	{
		assert(std::distance(begin1, end1) == std::distance(begin2, end2));
		return statistic::pearson_correlation<
		    typename std::iterator_traits<Iterator1>::value_type>(begin1, end1,
		                                                          begin2, end2);
	}
};

class GT2_EXPORT SpearmanCorrelation
{
  public:
	SpearmanCorrelation() {}

	template <typename Iterator1, typename Iterator2>
	typename std::iterator_traits<Iterator1>::value_type
	compute(Iterator1 begin1, Iterator1 end1, Iterator2 begin2,
	        Iterator2 end2) const
	{
		assert(std::distance(begin1, end1) == std::distance(begin2, end2));
		return statistic::spearman_correlation<
		    typename std::iterator_traits<Iterator1>::value_type>(begin1, end1,
		                                                          begin2, end2);
	}
};

class GT2_EXPORT KendallCorrelation
{
  public:
	KendallCorrelation() {}

	template <typename Iterator1, typename Iterator2>
	typename std::iterator_traits<Iterator1>::value_type
	compute(Iterator1 begin1, Iterator1 end1, Iterator2 begin2,
	        Iterator2 end2) const
	{
		assert(std::distance(begin1, end1) == std::distance(begin2, end2));
		return statistic::spearman_correlation<
		    typename std::iterator_traits<Iterator1>::value_type>(begin1, end1,
		                                                          begin2, end2);
	}
};
}

#endif // GT2_REGULATOR_ASSOCIATION_SCORE_H
