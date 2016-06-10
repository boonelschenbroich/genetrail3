/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_CORE_GENE_SET_ENRICHMENT_ANALYSIS_H
#define GT2_CORE_GENE_SET_ENRICHMENT_ANALYSIS_H

#include "macros.h"

#include "Category.h"

#include <boost/math/special_functions/binomial.hpp>

#include <functional>
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <utility>
#include <tuple>

namespace GeneTrail
{
	template <typename float_type, typename big_int_type>
	class GT2_EXPORT GeneSetEnrichmentAnalysis
	{
		constexpr static bool debug = false;

		public:
		template<typename Iterator>
		size_t intersectionSize(const Category& category, Iterator begin, const Iterator& end)
		{
			size_t n = 0;
			for(; begin != end; ++begin) {
				if(category.contains(*begin)) {
					++n;
				}
			}
			return n;
		}

		big_int_type absMax(big_int_type a, big_int_type b)
		{
			return std::abs(a) < std::abs(b) ? b : a;
		}

		/**
		 * This method computes the running sum statistic, based on the given
		 * categories.
		 *
		 * @param category Category for which the RSc should be computed.
		 * @param begin Iterator pointing to the beginning of a sorted list of genes,
		 *              based on which the RSc should be computed.
		 * @param end   Iterator pointing to the end of the gene list.
		 *
		 * @return The maximum value of the running sum for the given category.
		 */
		template<typename Iterator>
		big_int_type computeRunningSum(const Category& category,
									   Iterator begin, const Iterator& end)
		{
			size_t n = std::distance(begin, end);
			size_t l = intersectionSize(category, begin, end);
			size_t nl = n - l;
			big_int_type RSc = 0;
			big_int_type rs = 0;

			for(; begin != end; ++begin) {
				if(category.contains(*begin)) {
					rs += nl;
					RSc = absMax(RSc, rs);
				} else {
					rs -= l;
					RSc = absMax(RSc, rs);
				}
			}

			return RSc;
		}

		template<typename Iterator>
		big_int_type computeRunningSum(size_t n, Iterator begin, const Iterator& end)
		{
			if(begin == end) {
				return 0;
			}

			size_t l = std::distance(begin, end);
			size_t nl = n - l;

			big_int_type RSc = *begin * l;
			big_int_type rs = -RSc;

			rs += nl;

			RSc = absMax(RSc, rs);

			size_t lastIndex = *begin;

			for(++begin; begin != end; lastIndex = *begin, ++begin) {
				rs -= (*begin - lastIndex - 1) * l;
				RSc = absMax(RSc, rs);
				rs += nl;
				RSc = absMax(RSc, rs);
			}

			rs -= (n - lastIndex - 1) * l;
			RSc = absMax(RSc, rs);

			return RSc;
		}

		/**
		 * This method computes the running sum statistic and the corresponding
		 * p-value, based on the given categories.
		 *
		 * @param category Category for which the p-value should be computed.
		 * @param testSet Sorted list of genes, based on which the p-value
		 *                should be computed.
		 *
		 * @return The p-value for the given categories.
		 */
		float_type
		computeTwoSidedPValue(const Category& category,
		                      const std::vector<std::string>& testSet)
		{
			big_int_type RSc = computeRunningSum(category, testSet.begin(), testSet.end());
			return computeTwoSidedPValue(
			    testSet.size(), intersectionSize(category, testSet.begin(), testSet.end()), RSc);
		}

		/**
		 * This method computes a two-sided p-value for the given running sum
		 * statistic via dynamic programming.
		 *
		 * @param n The number of genes in the test set
		 * @param l The number of genes in the category
		 * @param RSc The running sum statistic
		 * @return p-value
		 */
		float_type computeTwoSidedPValue(const size_t& n, const size_t& l,
		                                 const big_int_type& RSc)
		{
			return computePValue_(
			    n, l, abs(RSc),
			    [](big_int_type v, big_int_type RSc) { return abs(v) < RSc; });
		}

		/**
		 * This method computes a lower-tailed p-value for the given running sum
		 * statistic via dynamic programming.
		 *
		 * @param n The number of genes in the test set
		 * @param l The number of genes in the category
		 * @param RSc The running sum statistic
		 * @return p-value
		 */
		float_type computeLeftPValue(const size_t& n, const size_t& l,
		                             const big_int_type& RSc)
		{
			assert(n >= l);
			return computePValue_(
			    n, l, std::abs(RSc),
			    [](big_int_type v, big_int_type RSc) { return -RSc < v; });
		}

		/**
		 * This method computes a upper-tailed p-value for the given running
		 * sum statistic via dynamic programming.
		 *
		 * @param n The number of genes in the test set
		 * @param l The number of genes in the category
		 * @param RSc The running sum statistic
		 * @return p-value
		 */
		float_type computeRightPValue(const size_t& n, const size_t& l,
		                              const big_int_type& RSc)
		{
			assert(n >= l);
			return computePValue_(
			    n, l, std::abs(RSc),
			    [](big_int_type v, big_int_type RSc) { return v < RSc; });
		}

		/**
		 * Internal routine for the pvalue computation.
		 * It implements the dynamic programming recursion.
		 *
		 * M[k,i] = M[k-1,i] + M[k,i-1] if comp(i*(n-l) + k*l) == true
		 * M[k,i] = 0 otherwise
		 *
		 * We exploit the fact that we can reuse the column i-1 in the
		 * computation of column i. This way we need less storage and can
		 * perform an in-place update.
		 *
		 * @todo A possible improvement would be to compute the complement
		 *       of this value for small running sums, as small input RScs
		 *       lead to a high number of additions.
		 */
		template <typename Comparator>
		float_type computePValue_(const size_t n, const size_t l,
		                          const big_int_type& RSc, Comparator comp)
		{
			const size_t nl = n - l;

			std::vector<float_type> M(nl + 1, 0);
			M[0] = 1;

			if(debug) {
				std::cout << M[0] << " ";
			}

			// Initialize
			for(size_t k = 1; k <= nl; ++k) {
				if(comp(-k * l, RSc)) {
					M[k] = 1;
				}
				if(debug) {
					std::cout << M[k] << " ";
				}
			}

			if(debug) {
				std::cout << std::endl;
			}

			for(size_t i = 1; i <= l; ++i) {
				big_int_type inl = i * nl;
				M[0] = comp(inl, RSc) ? 1 : 0;
				if(debug) {
					std::cout << M[0] << " ";
				}
				for(size_t k = 1; k <= nl; ++k) {
					if(comp(inl - k * l, RSc)) {
						M[k] += M[k - 1]; // This is the bottleneck
					} else {
						M[k] = 0;
					}
					if(debug) {
						std::cout << M[k] << " ";
					}
				}
				if(debug) {
					std::cout << std::endl;
				}
			}

			float_type binom =
			    boost::math::binomial_coefficient<float_type>(n, l);

			if(M[nl] >= binom) {
				return 0.0;
			}
			return 1.0 - M[nl] / binom;
		}

		/**
		 * This method computes the running sum statistic and the corresponding
		 *p-value, based on the given categories.
		 *
		 * @param category Category for which the p-value should be computed.
		 * @param testSet Sorted list of genes, based on which the p-value
		 *should be computed.
		 * @return The p-value for the given categories.
		 */
		float_type
		computeOneSidedPValue(const Category& category,
		                      const std::vector<std::string>& testSet)
		{
			big_int_type RSc = computeRunningSum(category, testSet.begin(), testSet.end());
			if(RSc < 0) {
				return computeLeftPValue(
				    testSet.size(), intersectionSize(category, testSet.begin(), testSet.end()), RSc);
			} else {
				return computeRightPValue(
				    testSet.size(), intersectionSize(category, testSet.begin(), testSet.end()), RSc);
			}
		}

		/**
		 * DEPRECATED - Please use the new implementation
		 * This method computes the running sum statistic and the corresponding
		 *p-value, based on the given categories.
		 *
		 * @param category Category for which the p-value should be computed.
		 * @param testSet Sorted list of genes, based on which the p-value
		 *should be computed.
		 * @return The p-value for the given categories.
		 */
		float_type
		computeOneSidedPValueD(const Category& category,
		                       const std::vector<std::string>& testSet)
		{
			big_int_type RSc = computeRunningSum(category, testSet);
			return computeOneSidedPValueD(
			    testSet.size(), intersectionSize(category, testSet), RSc);
		}

		/**
		 * DEPRECATED - Please use the new implementation
		 * This method computes a one-sided p-value for the given running sum
		 *statistic.
		 * (This is not the original implementation of the Paper)
		 *
		 * @param n The number of genes in the test set
		 * @param l The number of genes in the category
		 * @param RSc The running sum statistic
		 * @return p-value
		 */
		float_type computeOneSidedPValueD(const size_t& n, const size_t& l,
		                                  const big_int_type& RSc)
		{
			std::map<big_int_type, float_type> pathways;
			pathways[0] = 1;

			bool enriched = RSc >= 0;
			for(size_t i = 0; i < n; i++) {
				iterate(pathways, i, RSc, n, l, enriched);
			}

			auto iter = pathways.find(0);

			if(iter != pathways.end()) {
				float_type binom =
				    boost::math::binomial_coefficient<float_type>(n, l);
				return 1.0 - iter->second / binom;
			}
			return 1.0;
		}

		/**
		 * This method is needed to fill the map with the number of possible
		 *paths.
		 *
		 * @param newpathway Temporary copy of the pathway map
		 * @param value The number of possible paths at this position.
		 * @param newpos The new position in the thought dynamic programming
		 *matrix
		 */
		void fillNewPathway(std::map<big_int_type, float_type>& newpathway,
		                    float_type value, const big_int_type& newpos)
		{
			auto iter1 = newpathway.find(newpos);
			if(iter1 == newpathway.end()) {
				newpathway[newpos] = value;
			} else {
				newpathway[newpos] = iter1->second + value;
			}
		}

		/**
		 * This method iterates over the previous layer of the thought dynamic
		 *programming matrix and fills the current one.
		 *
		 * @param pathway The previous layer of the thought dynamic programming
		 *matrix
		 * @param i The index of the current layer
		 * @param RSc The value of the running sum statistic
		 * @param n The number of genes in the test set
		 * @param l The number of genes in the category
		 * @param enriched Boolean flag indicating if we consider a enriched or
		 *depleted
		 */
		void iterate(std::map<big_int_type, float_type>& pathway,
		             const big_int_type& i, const big_int_type& RSc,
		             const size_t& n, const size_t& l, const bool& enriched)
		{
			auto iter = pathway.begin();
			std::map<big_int_type, float_type> newpathway;

			while(iter != pathway.end()) {
				// where could we be in the next step
				// which balls are available
				big_int_type black = (iter->first + i * l) / n;
				big_int_type white = i - black;

				// having drawn a black ball
				if(black < l) {
					big_int_type newpos = iter->first + n - l;
					// only relevant if we are smaller RSc
					if(!enriched || (newpos < RSc)) {
						fillNewPathway(newpathway, iter->second, newpos);
					}
				}

				// having drawn a white ball
				if(white < (n - l)) {
					big_int_type newpos = iter->first - l;
					if(enriched) {
						// checking whether we can still reach zero
						big_int_type black = (newpos + (i + 1) * l) / n;
						if((newpos + (l - black) * (n - l) >= 0) &&
						   (black <= l)) {
							fillNewPathway(newpathway, iter->second, newpos);
						}
					} else {
						if(newpos > RSc) {
							fillNewPathway(newpathway, iter->second, newpos);
						}
					}
				}
				++iter;
			}
			pathway = newpathway;
			newpathway.clear();
		}
	};
}

#endif // GT2_CORE_GENE_SET_ENRICHMENT_ANALYSIS_H

