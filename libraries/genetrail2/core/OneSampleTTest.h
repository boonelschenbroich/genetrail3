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
#ifndef GT2_CORE_ONE_SAMPLE_TTEST_H
#define GT2_CORE_ONE_SAMPLE_TTEST_H

#include "macros.h"

#include "Statistic.h"

#include <boost/math/distributions/students_t.hpp>

#include <cmath>
#include <iterator>
#include <iostream>

namespace GeneTrail {

    /**
     * One Sample Student's T-Test
     */
    template <typename ValueType>
    class GT2_EXPORT OneSampleTTest {
    public:
		using value_type = ValueType;

        OneSampleTTest(value_type tol = 1e-5, value_type mu = 0.0) : tolerance_(tol), mu_(mu) {
        }

		template<typename InputIterator1, typename InputIterator2>
		value_type test(const InputIterator1& begin, const InputIterator1& end, const InputIterator2&, const InputIterator2&) {
			return test(begin, end);
		}

        /**
         * This method implements the "One Sample Student's T-Test".
         *
         * @param begin Iterator indicating where to begin
         * @param end Iterator indicating where to end
         * @return t-score
         */
        template<typename InputIterator>
        value_type test(const InputIterator& begin, const InputIterator& end) {
            auto mean = statistic::mean<value_type>(begin, end);
            auto var = statistic::var<value_type>(begin, end);
			auto size = std::distance(begin, end);

            this->stdErr_ = std::sqrt(var / size);
			this->df_ = size - 1;

            if (stdErr_ < tolerance_) {
                score_ = mu_;
            } else {
                score_ = (mean - mu_) / stdErr_;
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

#endif // GT2_CORE_ONE_SAMPLE_TTEST_H

