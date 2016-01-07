/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2014-2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *               2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_HTEST_H
#define GT2_CORE_HTEST_H

#include <boost/math/distributions.hpp>

namespace GeneTrail
{
	namespace HTest
	{
		template <typename Test, typename InputIterator1,
		          typename InputIterator2>
		typename Test::value_type test(Test& t,
		                               const InputIterator1& first_begin,
		                               const InputIterator1& first_end,
		                               const InputIterator2& second_begin,
		                               const InputIterator2& second_end)
		{
			return t.test(first_begin, first_end, second_begin, second_end);
		}

		template <typename Test, typename InputIterator>
		typename Test::value_type test(Test& t, const InputIterator& begin,
		                               const InputIterator& end)
		{
			return t.test(begin, end);
		}

		template <typename Test>
		typename Test::value_type
		twoSidedPValue(Test& t, const typename Test::value_type& score)
		{
			return 2.0 * boost::math::cdf(boost::math::complement(
			                 t.distribution(), fabs(score)));
		}

		template <typename Test>
		std::pair<typename Test::value_type, typename Test::value_type>
		confidenceInterval(Test& t, const typename Test::value_type& alpha)
		{
			return t.confidenceInterval(alpha);
		}

		template <typename Test>
		typename Test::value_type
		lowerTailedPValue(Test& t, const typename Test::value_type& score)
		{
			return boost::math::cdf(t.distribution(), score);
		}

		template <typename Test>
		typename Test::value_type
		upperTailedPValue(Test& t, const typename Test::value_type& score)
		{
			return boost::math::cdf(
			    boost::math::complement(t.distribution(), score));
		}
	}
}

#endif // GT2_CORE_HTEST_H
