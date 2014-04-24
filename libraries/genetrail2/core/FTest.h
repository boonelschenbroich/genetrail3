/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2013 Tim Kehl tkehl@bioinf.uni-sb.de>
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
#ifndef F_TTEST_H
#define F_TTEST_H

#include "config.h"

#include "Statistic.h"

#include <boost/math/distributions.hpp>

namespace GeneTrail {

    /**
     * F-Test
     */
    template <typename value_type, typename InputIterator1,typename InputIterator2>
    class GT2_EXPORT FTest {
    public:

        FTest(value_type tol = 1e-4) : tolerance_(tol),size1_(0),size2_(0),score_(0) {
        }

        /**
         * This method implements the standerd F-Test.
         *
         * @param Iterator
         * @param Iterator
		 * @param Iterator
		 * @param Iterator
         * @return F-score for the differences between the two groups
         */
        value_type test(const InputIterator1& first_begin, const InputIterator1& first_end, const InputIterator2& second_begin, const InputIterator2& second_end) {
			this->size1_ = std::distance(first_begin,first_end);
			this->size2_ = std::distance(second_begin,second_end);
			value_type var1 = statistic::var<value_type, InputIterator1>(first_begin,first_end);
            value_type var2 = statistic::var<value_type, InputIterator2>(second_begin,second_end);
            score_ = var1 / var2;
            return score_;
        }

		boost::math::fisher_f distribution(){
			boost::math::fisher_f dist(size1_ - 1, size2_ - 1);
			return dist;
		}

		value_type score(){
			score_;
		}

		std::pair<value_type, value_type> confidenceInterval(const value_type& alpha) {
			value_type conf = boost::math::quantile(boost::math::complement(distribution(), (1 - alpha) / 2.0));
			return std::make_pair(score_ / conf, score_ * conf);
		}

    protected:

        value_type tolerance_;
		size_t size1_;
		size_t size2_;
		value_type score_;
    };

}

#endif // F_TEST_H

