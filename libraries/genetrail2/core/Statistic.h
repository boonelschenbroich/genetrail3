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
#ifndef STATISTIC_NEW_H
#define STATISTIC_NEW_H

#include <algorithm>

namespace GeneTrail {

    /**
     * A collection of mathematical operations.
     */
    namespace statistic {

        /**
         * This method calculates the maximal value of a given container.
         *
         * @param c Standard C++ container holding values of type Type
         * @return Maximum of container c
         */
        template <typename value_type, typename InputIterator>
        value_type max(InputIterator begin, InputIterator end) {
            return *std::max_element(begin, end);
        }

        /**
         * This method calculates the minimal value of a given container.
         *
         * @param c Standard C++ container holding values of type Type
         * @return Minimum of container c
         */
		template <typename value_type, typename InputIterator>
        value_type min(InputIterator begin, InputIterator end) {
            return *std::min_element(begin, end);
        }

        /**
         * This method calculates the mean of a given container.
         *
         * @param c Standard C++ container holding values of type Type
         * @return Mean of container c
         */
        template <typename value_type, typename InputIterator>
        double mean(InputIterator begin, InputIterator end) {
            return (std::accumulate(begin, end,  0.0) / ((value_type)std::distance(begin,end)));
        }

        /**
         * This method calculates the median of a given container.
         *
         * @param c Standard C++ container holding values of type Type
         * @return Median of container c
         */
        template <typename value_type, typename InputIterator>
        value_type median(InputIterator begin, InputIterator end) {
            std::nth_element(begin, begin + (((std::distance(begin,end) + 1) / 2) - 1), end);
            return *(begin + (((std::distance(begin,end) + 1) / 2) - 1));
        }

        /**
         * This method calculates the pooled variance of a given container.
         *
         * @param c Standard C++ container holding values of type Type
         * @return Variance of container c
         */
        template <typename value_type, typename InputIterator>
        value_type var(InputIterator begin, InputIterator end) {
            int n = std::distance(begin,end);
            if (n == 1) {
                return ((value_type) 0);
            } else {
                value_type s = std::accumulate(begin, end, ((value_type) 0));
                value_type mean = s / n;
                value_type v = (value_type) 0;
				value_type ep = (value_type) 0;
                for (auto it = begin; it != end; ++it) {
                    s = *it - mean;
                    ep += s;
                    v += s*s;
                }
                return (v - ep * ep / n) / (n - 1);
            }
        }

        /**
         * This method calculates the standard deviation of a given container.
         *
         * @param c Standard C++ container holding values of type Type
         * @return Standard deviation of container c
         */
        template <typename value_type, typename InputIterator>
        value_type sd(InputIterator begin, InputIterator end) {
            return std::sqrt(var(begin,end));
        }
    };
}

#endif // STATISTIC_H

