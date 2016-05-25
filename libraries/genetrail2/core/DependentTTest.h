/*
 * GeneTrail2 - An efficient library for interpreting genetic data
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
#ifndef GT2_CORE_DEPENDENT_TTEST_H
#define GT2_CORE_DEPENDENT_TTEST_H

#include <list>
#include <cmath>

#include "macros.h"

#include "Statistic.h"

namespace GeneTrail {

    /**
     * Dependent Student's T-Test
     */
    template <typename ValueType>
    class GT2_EXPORT DependentTTest {
    public:
		using value_type = ValueType;

        DependentTTest(value_type tol = 1e-5, value_type mu = 0.0) : tolerance_(tol), mu_(mu) {
        }

		/**
		 * Dependent Student's T-Test
		 *
		 * @param first_begin Iterator to the begin of the first range of sample
		 *        values.
		 * @param first_end Iterator to the end of the first range of sample
		 *        values.
		 * @param second_begin Iterator to the begin of the second range of
		 *        sample values.
		 * @param second_end Iterator to the end of the second range of sample
		 *        values.
		 *
		 * @return T-score for the differences between the two groups
		 *
		 * @throws std::out_of_range If the number of elements in the input
		 *         ranges do not match
		 */
		template<typename InputIterator1, typename InputIterator2>
        value_type test(const InputIterator1& first_begin, const InputIterator1& first_end, const InputIterator2& second_begin, const InputIterator2& second_end) {
			const size_t num_entries1 = std::distance(first_begin, first_end);
			const size_t num_entries2 = std::distance(second_begin, second_end);

			if(num_entries1 != num_entries2) {
				throw std::out_of_range(
				    "The number of entries in the input ranges do not match. "
				    "Range 1: " +
				    boost::lexical_cast<std::string>(num_entries1) +
				    ", Range 2: " +
				    boost::lexical_cast<std::string>(num_entries2));
			}


			std::vector<value_type> diff(num_entries1);
			std::transform(first_begin, first_end, second_begin, diff.begin(),
			               [](value_type a, value_type b) { return a - b; });

            auto var = statistic::var<value_type>(diff.begin(), diff.end());
			stdErr_ = std::sqrt(var / num_entries1);

            if (std::fabs(var / num_entries1) < tolerance_) {
                score_ = mu_;
            } else {
				score_ = (statistic::mean<value_type>(diff.begin(),diff.end()) - mu_) / stdErr_;
            }

			df_ = num_entries1 - 1;
            return score_;
        }

		/**
		 * Get the used students_t distribution.
		 *
		 * @warning The behavior ofthis function is undefined until test has been called.
		 * @returns An instance of boost::math::students_t.
		 */
		boost::math::students_t distribution() const {
			boost::math::students_t dist(this->df_);
			return dist;
		}

		/**
		 * Compute a confidence interval of level alpha.
		 *
		 * @param alpha The confidence level to be used for computing the
		 *              interval.
		 *
		 * @returns A pair containing the left and right bounds of the confidence interval.
		 */
		std::pair<value_type, value_type> confidenceInterval(const value_type& alpha) {
			boost::math::students_t dist(df_);
			value_type conf = boost::math::quantile(boost::math::complement(dist, (1 - alpha) / 2.0));
			return std::make_pair((score_ - conf) * stdErr_, (score_ + conf) * stdErr_);
		}

    protected:
        value_type tolerance_;
		value_type df_;
		value_type stdErr_;
		value_type score_;
		value_type mu_;
    };
}

#endif // GT2_CORE_DEPENDENT_TTEST_H

