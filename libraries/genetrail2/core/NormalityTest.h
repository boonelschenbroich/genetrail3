/*
 * GeneTrail2 - An efficient library for interpreting genetic data
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
#ifndef GT2_CORE_NORMALITY_TEST_H
#define GT2_CORE_NORMALITY_TEST_H

#include <boost/math/distributions.hpp>

#include "AndersonDarlingTest.h"

namespace GeneTrail
{
	/**
	 * A collection of normality tests based on the Anderson-Darling Test.
	 */
	namespace NormalityTest
	{
		template <typename value_type, typename InputIterator,
		          typename Distribution>
		bool test(Distribution dist, const InputIterator& begin,
		          const InputIterator& end, value_type alpha)
		{
			AndersonDarlingTest<value_type, Distribution> ADTest(dist);
			ADTest.test(begin, end);
			return ADTest.pValue() > alpha;
		}

		template <typename value_type, typename InputIterator>
		bool testForNormalDistribution(const InputIterator& begin,
		                               const InputIterator& end,
		                               value_type alpha)
		{
			auto mean = statistic::mean<value_type>(begin, end);
			auto sd = statistic::sd<value_type>(begin, end);
			boost::math::normal dist(mean, sd);
			return test(dist, begin, end, alpha);
		}

		template <typename value_type, typename InputIterator>
		bool testForStandardNormalDistribution(const InputIterator& begin,
		                                       const InputIterator& end,
		                                       value_type alpha)
		{
			boost::math::normal dist(0, 1);
			return test(dist, begin, end, alpha);
		}
	}
}

#endif // GT2_CORE_NORMALITY_TEST_H

