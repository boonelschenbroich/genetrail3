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

#include <cmath>
#include <iterator>

#include "macros.h"

#include "Statistic.h"

namespace GeneTrail {

    /**
     * One Sample Student's T-Test
     */
    template <typename value_type>
    class GT2_EXPORT OneSampleTTest {
    public:

        OneSampleTTest(value_type tol = 1e-5, value_type mu = 0.0) : tolerance_(tol), mu_(mu) {
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
            value_type mean = statistic::mean<value_type, InputIterator>(begin, end);
            value_type var = statistic::var<value_type, InputIterator>(begin, end);
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

