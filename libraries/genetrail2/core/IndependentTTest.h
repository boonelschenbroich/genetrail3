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
#ifndef GT2_CORE_INDEPENDENT_TTEST_H
#define GT2_CORE_INDEPENDENT_TTEST_H

#include <cmath>

#include "macros.h"

#include "Statistic.h"

#include <boost/math/distributions/students_t.hpp>

namespace GeneTrail {

    /**
     * Independent Student's T-Test
     */
    template <typename ValueType>
    class GT2_EXPORT IndependentTTest {
    public:
		using value_type = ValueType;

        IndependentTTest(value_type tol = 1e-5, value_type mu = 0.0) : tolerance_(tol), mu_(mu) {
        }

        /**
         * This method implements the "Independent Student's T-Test".
         * This tests checks if two distributions have the same mean.
         * The data arrays are allowed to be drawn from population with unequal variances.
         *
         * NUMERICAL RECIPES page 729
         *
         * @param first_begin Iterator
         * @param first_end Iterator
		 * @param second_begin Iterator
		 * @param second_end Iterator
         * @return t-score for the differences between the two groups
         */
        template<typename InputIterator1, typename InputIterator2>
        value_type test(const InputIterator1& first_begin, const InputIterator1& first_end, const InputIterator2& second_begin, const InputIterator2& second_end) {
			value_type mean1, mean2, var1, var2;
			size_t size1, size2;

			std::tie(mean1, var1, size1) = statistic::mean_var_size<value_type>(first_begin, first_end);
			std::tie(mean2, var2, size2) = statistic::mean_var_size<value_type>(second_begin, second_end);

			auto stdErr1_2 = var1 / size1;
			auto stdErr2_2 = var2 / size2;
			auto stdErr_2 = stdErr1_2 + stdErr2_2;
			this->stdErr_ = std::sqrt(stdErr_2);
			this->df_ = (stdErr_2 * stdErr_2) / (stdErr1_2*stdErr1_2 / (size1 - 1) + stdErr2_2*stdErr2_2 / (size2 - 1));

            if (stdErr_ < tolerance_) {
                score_ = 0;
            } else {
                score_ = (mean1 - mean2 - mu_) / stdErr_;
            }
			return score_;
        }

		boost::math::students_t distribution() {
			boost::math::students_t dist(this->df_);
			return dist;
		}

		std::pair<value_type, value_type> confidenceInterval(const value_type& alpha) {
			boost::math::students_t dist(df_);
			value_type conf = boost::math::quantile(boost::math::complement(dist, (1 - alpha) / 2.0));
			return std::make_pair((score_ - conf) * stdErr_, (score_ + conf) * stdErr_);
		}

    protected:

        value_type tolerance_;
		value_type df_;
		value_type mu_;
		value_type stdErr_;
		value_type score_;
    };
}

#endif // GT2_CORE_INDEPENDENT_TTEST_H

