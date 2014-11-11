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

