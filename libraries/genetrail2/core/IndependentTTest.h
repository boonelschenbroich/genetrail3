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
#ifndef INDEPENDENT_TTEST_H
#define INDEPENDENT_TTEST_H

#include <cmath>

#include "macros.h"

#include "Statistic.h"

namespace GeneTrail {

    /**
     * Independent Student's T-Test
     */
    template <typename value_type, typename InputIterator1,typename InputIterator2>
    class GT2_EXPORT IndependentTTest {
    public:

        IndependentTTest(value_type tol = 1e-5) : tolerance_(tol) {
        }

        /**
         * This method implements the "Independent Student's T-Test".
         * This tests checks if two distributions have the same mean.
         * The data arrays are allowed to be drawn from population with unequal variances.
         *
         * NUMERICAL RECIPES page 729
         *
         * @param Iterator
         * @param Iterator
		 * @param Iterator
		 * @param Iterator
         * @return A PValue object containing a t-score for the differences between the two groups
         */
        value_type test(const InputIterator1& first_begin, const InputIterator1& first_end, const InputIterator2& second_begin, const InputIterator2& second_end) {
            value_type mean1 = statistic::mean<value_type, InputIterator1>(first_begin,first_end);
            value_type mean2 = statistic::mean<value_type, InputIterator2>(second_begin,second_end);
            value_type var1 = statistic::var<value_type, InputIterator1>(first_begin,first_end);
            value_type var2 = statistic::var<value_type, InputIterator2>(second_begin,second_end);

            size_t size1 = std::distance(first_begin,first_end);
            size_t size2 = std::distance(second_begin,second_end);

            value_type stdErr1 = std::sqrt(var1 / size1);
            value_type stdErr2 = std::sqrt(var2 / size2);
            this->stdErr_ = std::sqrt(std::pow(stdErr1, 2) + std::pow(stdErr2, 2));
			this->df_ = std::pow(stdErr_, 4) / (std::pow(stdErr1, 4) / (size1 - 1) + std::pow(stdErr2, 4) / (size2 - 1));

            if (stdErr_ < tolerance_) {
                score_ = 0;
            } else {
                score_ = (mean1 - mean2) / stdErr_;
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
		value_type stdErr_;
		value_type score_;
    };
}

#endif // INDEPENDENT_TTEST_H

