/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
 *
 * This program is free software: you can redist_ribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is dist_ributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public
 * License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef GT2_CORE_ANDERSON_DARLING_TEST_H
#define GT2_CORE_ANDERSON_DARLING_TEST_H

#include <genetrail2/core/macros.h>
#include <genetrail2/core/Statistic.h>

#include <boost/math/distributions.hpp>

#include <cmath>

namespace GeneTrail
{

	/**
	 * Anderson-Darling Test for normality.
	 */
	template <typename value_type, typename Distribution>
	class GT2_EXPORT AndersonDarlingTest
	{
		public:
		/**
		 * Default constructor.
		 *
		 * @param d
		 * @todo Clarify what the distribution is for...
		 */
		AndersonDarlingTest(Distribution d)
		    : dist_(d)
		{
		}

		/**
		 * This method implements a variant of the Anderson-Darling Test.
		 * This variant assumes that mean and variance of the data are both
		 * unknown.
		 * The test checks if the given values follow a standard normal
		 * distribution.
		 *
		 * @param begin The begining of the range to which the test should be applied.
		 * @param end The end of the range to which the test should be applied.
		 * @return A2* - score
		 */
		template <typename InputIterator>
		value_type test(const InputIterator& begin, const InputIterator& end)
		{
			using namespace boost::math;

			std::vector<value_type> scores(begin, end);

			// Sort Z-scores
			std::sort(scores.begin(), scores.end());

			value_type sum = 0.0;
			const value_type n = scores.size();
			for(unsigned int i = 0; i < n; ++i) {
				sum += (2 * (i + 1) - 1) *
				       (std::log(cdf(dist_, scores[i])) +
				        std::log(cdf(complement(dist_, scores[n - i - 1]))));
			}

			value_type A2 = -n - (1 / n) * sum;
			score_ = A2 * (1 + 0.75 / n + 2.25 / (n * n));
			return score_;
		}

		/**
		 * The p-value is computed according to Table 4.9 in Stephens (1986).
		 * See R package "nortest"
		 *
		 * @return The computed p-value.
		 */
		value_type pValue()
		{
			value_type p;
			if(score_ < 0.2) {
				p = 1.0 - std::exp(-13.436 + 101.14 * score_ -
				                   223.73 * std::pow(score_, 2.0));
			} else if(score_ < 0.34) {
				p = 1.0 - std::exp(-8.318 + 42.796 * score_ -
				                   59.938 * std::pow(score_, 2.0));
			} else if(score_ < 0.6) {
				p = std::exp(0.9177 - 4.279 * score_ -
				             1.38 * std::pow(score_, 2.0));
			} else {
				p = std::exp(1.2937 - 5.709 * score_ +
				             0.0186 * std::pow(score_, 2.0));
			}

			return p;
		}

		protected:
		value_type score_;
		Distribution dist_;
	};
}

#endif // GT2_CORE_ANDERSON_DARLING_TEST_H

