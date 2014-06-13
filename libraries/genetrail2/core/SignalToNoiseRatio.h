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
#ifndef GT2_CORE_SIGNAL_TO_NOISE_RATIO_H
#define GT2_CORE_SIGNAL_TO_NOISE_RATIO_H

#include "macros.h"

#include "Statistic.h"

#include <boost/math/distributions.hpp>

namespace GeneTrail {

    /**
     * SignalToNoiseRatio
     */
    template <typename value_type, typename InputIterator1,typename InputIterator2>
    class GT2_EXPORT SignalToNoiseRatio {
    public:

        SignalToNoiseRatio(value_type tol = 1e-4) : tolerance_(tol),score_(0) {
        }

        /**
         * This method implements the standerd F-Test.
         *
         * @param Iterator
         * @param Iterator
		 * @param Iterator
		 * @param Iterator
         * @return Signal to noise ratio between the two groups
         */
        value_type test(const InputIterator1& first_begin, const InputIterator1& first_end, const InputIterator2& second_begin, const InputIterator2& second_end) {
			value_type mean1 = statistic::mean<value_type, InputIterator1>(first_begin,first_end);
			value_type mean2 = statistic::mean<value_type, InputIterator1>(first_begin,first_end);
			value_type sd1 = statistic::sd<value_type, InputIterator1>(first_begin,first_end);
            value_type sd2 = statistic::sd<value_type, InputIterator2>(second_begin,second_end);
            score_ = (mean1 - mean2) / (sd1 + sd2);
            return score_;
        }

		boost::math::fisher_f distribution(){
			boost::math::normal dist(0,1);
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
		value_type score_;
    };

}

#endif // GT2_CORE_SIGNAL_TO_NOISE_RATIO_H

