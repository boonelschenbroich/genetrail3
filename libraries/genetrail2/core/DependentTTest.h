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
    template <typename value_type>
    class GT2_EXPORT DependentTTest {
    public:

        DependentTTest(value_type tol = 1e-5, value_type mu = 0.0) : tolerance_(tol), mu_(mu) {
        }

        /**
         *
         * Dependent Student's T-Test
         *
         * @param Iterator
         * @param Iterator
		 * @param Iterator
		 * @param Iterator
         * @return T-score for the differences between the two groups
         */
		template<typename InputIterator1, typename InputIterator2>
        value_type test(const InputIterator1& first_begin, const InputIterator1& first_end, const InputIterator2& second_begin, const InputIterator2& second_end) {

			std::list<value_type> diff(first_begin,first_end);
            auto it2 = second_begin;

            for (auto it=diff.begin(); it != diff.end(); ++it) {
                *it -= *it2;
                ++it2;
            }

            auto var = statistic::var<value_type>(diff.begin(), diff.end());
            size_t size = diff.size();
			stdErr_ = std::sqrt(var / size);

            if (std::fabs(var / size) < tolerance_) {
                score_ = mu_;
            } else {
				score_ = (statistic::mean<value_type>(diff.begin(),diff.end()) - mu_) / stdErr_;
            }

			df_ = size - 1;
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
		value_type mu_;
    };
}

#endif // GT2_CORE_DEPENDENT_TTEST_H

