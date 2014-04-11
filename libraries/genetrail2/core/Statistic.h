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

#include <cmath>
#include <algorithm>
#include <vector>

namespace GeneTrail {

    /**
     * A collection of mathematical operations.
     */
    namespace statistic {

		/**
		 * This method calculates the absolute value for each entry of a given container.
		 *
		 * @param begin InputIterator
		 * @param end InputIterator
		 */
		template <typename value_type, typename InputIterator>
		void abs(InputIterator begin, InputIterator end) {
			std::transform (begin, end , begin, static_cast< value_type(*)(value_type) >( std::abs ));
		}

		/**
		 * This method calculates the square root for each entry of a given container.
		 *
		 * @param begin InputIterator
		 * @param end InputIterator
		 */
		template <typename value_type, typename InputIterator>
		void sqrt(InputIterator begin, InputIterator end){
			std::transform (begin, end , begin, static_cast< value_type(*)(value_type) >( std::sqrt ));
		}

		/**
		 * This method calculates the natural logarithm for each entry of a given container.
		 *
		 * @param begin InputIterator
		 * @param end InputIterator
		 */
		template <typename value_type, typename InputIterator>
		void log(InputIterator begin, InputIterator end){
			std::transform (begin, end , begin, static_cast< value_type(*)(value_type) >( std::log ));
		}

		/**
		 * This method calculates the logarithm (base 10) for each entry of a given container.
		 *
		 * @param begin InputIterator
		 * @param end InputIterator
		 */
		template <typename value_type, typename InputIterator>
		void log10(InputIterator begin, InputIterator end){
			std::transform (begin, end , begin, static_cast< value_type(*)(value_type) >( std::log10 ));
		}

		/**
		 * This method calculates the logarithm (base 2) for each entry of a given container.
		 *
		 * @param begin InputIterator
		 * @param end InputIterator
		 */
		template <typename value_type, typename InputIterator>
		void log2(InputIterator begin, InputIterator end){
			std::transform (begin, end , begin, static_cast< value_type(*)(value_type) >( std::log2 ));
		}

		/**
		 * This method calculates the n-th power for each entry of a given container.
		 * @param begin InputIterator
		 * @param end InputIterator
		 */
		template <typename value_type, typename InputIterator>
		void pow(InputIterator begin, InputIterator end, int n){
			std::transform (begin, end , begin, [n](value_type x){ return std::pow(x,n);});
		}

        /**
         * This method calculates the maximal value of a given container.
         *
         * @param begin InputIterator
		 * @param end InputIterator
         * @return Maximum of the given range
         */
        template <typename value_type, typename InputIterator>
        value_type max(InputIterator begin, InputIterator end) {
            return *std::max_element(begin, end);
        }

        /**
         * This method calculates the minimal value of a given container.
         *
         * @param begin InputIterator
		 * @param end InputIterator
         * @return Minimum of the given range
         */
		template <typename value_type, typename InputIterator>
        value_type min(InputIterator begin, InputIterator end) {
            return *std::min_element(begin, end);
        }

        /**
         * This method calculates the mean of a given container.
         *
         * @param begin InputIterator
		 * @param end InputIterator
         * @return Mean of the given range
         */
        template <typename value_type, typename InputIterator>
        double mean(InputIterator begin, InputIterator end) {
            return (std::accumulate(begin, end,  0.0) / ((value_type)std::distance(begin,end)));
        }

		/**
		 * This method calculates the max-mean statistic of a given container.
		 *
		 * @param begin InputIterator
		 * @param end InputIterator
		 * @return Max-mean statistic of the given range
		 */
		template <typename value_type, typename InputIterator>
		value_type max_mean(InputIterator begin, InputIterator end){
			value_type positive_sum = 0;
			value_type negative_sum = 0;
			for(auto iter = begin; iter != end; ++iter) {
				if(*iter < 0) {
					negative_sum -= *iter;
				} else {
					positive_sum += *iter;
				}
			}
			value_type n = std::distance(begin,end);
			return std::max(negative_sum / n, positive_sum / n);
		}

		/**
         * This method calculates the median of a given container.
         *
         * @param begin InputIterator
		 * @param end InputIterator
         * @return Median of the given range
         */
        template <typename value_type, typename InputIterator>
        value_type median(InputIterator begin, InputIterator end) {
			std::vector<value_type> tmp(begin,end);
			unsigned int median_position = (((std::distance(tmp.begin(),tmp.end()) + 1) / 2) - 1);
            std::nth_element(tmp.begin(), tmp.begin() + median_position, tmp.end());
            return *(tmp.begin() + median_position);
        }

        /**
         * This method calculates the pooled variance of a given container.
         *
         * @param begin InputIterator
		 * @param end InputIterator
         * @return Variance of the given range
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
         * @param begin InputIterator
		 * @param end InputIterator
         * @return Standard deviation of the given range
         */
        template <typename value_type, typename InputIterator>
        value_type sd(InputIterator begin, InputIterator end) {
            return std::sqrt(var<value_type,InputIterator>(begin,end));
        }
    }
}

#endif // STATISTIC_H

