/*
 * GeneTrail2 - An efficient library for interpreting genetic data
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
#ifndef GT2_CORE_SHRINKAGE_T_TEST_H
#define GT2_CORE_SHRINKAGE_T_TEST_H

#include <vector>
#include <cmath>
#include <tuple>
#include <utility>

#include "macros.h"

#include "Statistic.h"

#include <boost/iterator/filter_iterator.hpp>

namespace GeneTrail {

	template <typename value_type>
	struct Gene
	{
		value_type mean;
		value_type var;
		value_type Var;
		value_type shrink_var;
		size_t size;
	};

	/**
	 * Shrinkage T-Test
	 *
	 * Implements the regularized t-test as Opgen-Rhein and Strimmer (2007)
	 * "Accurate ranking of differentially expressed genes by a
	 * distribution-free shrinkage approach"
	 *
	 * Given samples from a multivariate distribution, the variance used in the
	 * t-test is shrunk towards the median variance accross all variables.
	 *
	 * The shrinkage factor is estimated from the data.
	 */
	template <typename value_type>
	class GT2_EXPORT ShrinkageTTest {

		public:
		/**
		 * Computes the different variances for a single gene.
		 *
		 * @todo Use a compensated algorithm for computing the variance.
		 *
		 * @param begin InputIterator
		 * @param end InputIterator
		 * @return Gene struct holding all needed fields
		 */
		template <typename InputIterator>
		Gene<value_type> computeVariances(InputIterator begin,
		                                  InputIterator end)
		{
			const size_t size = std::distance(begin, end);
			wik_.resize(size);

			value_type mean = 0.0;
			for (size_t i = 0; i < size; ++i, ++begin) {
				mean += wik_[i] = *begin;
			}
			mean /= size;

			value_type sum_wik = 0.0;
			for(auto& item : wik_) {
				item -= mean;
				item *= item;
				sum_wik += item;
			}

			const auto wk = sum_wik / size;

			//TODO: Use compensated variance algorithm here
			value_type Var = (value_type) 0;
			for (const auto& item : wik_) {
				auto diff = item - wk;
				Var += diff*diff;
			}

			const size_t size_1 = size - 1;
			Var *= ((value_type)size) / (size_1 * size_1 * size_1);

			return Gene<value_type>{mean, sum_wik / size_1, Var, 0.0, size};
		}

		/**
		 * Computes variances for all genes.
		 *
		 * @param begin InputIterator
		 * @param end InputIterator
		 * @return Vector of genes
		 */
		template <typename InputIterator>
		std::tuple<std::vector<Gene<value_type>>,value_type> computeAllVariances(InputIterator begin, InputIterator end){
			const auto n = std::distance(begin, end);

			std::vector<Gene<value_type>> variances;
			std::vector<value_type> vars;

			variances.reserve(n);
			vars.reserve(n);

			for (auto iter = begin; iter != end; ++iter) {
				variances.emplace_back(computeVariances(iter->begin(), iter->end()));
				vars.emplace_back(variances.back().var);
			}

			return std::make_pair(variances, statistic::median<value_type>(vars.begin(), vars.end()));
		}

		/**
		 * Computes variances for all genes.
		 *
		 * @param begin InputIterator
		 * @param end InputIterator
		 * @return Vector of genes
		 */
		template <typename InputIterator>
		std::tuple<std::vector<Gene<value_type>>,value_type> computeAllVariancesRemoveNaN(InputIterator begin, InputIterator end) {
			const auto n = std::distance(begin, end);

			std::vector<Gene<value_type>> variances;
			std::vector<value_type> vars;

			variances.reserve(n);
			vars.reserve(n);

			auto is_not_nan = [](value_type x) { return !std::isnan(x); };
			for (auto iter = begin; iter != end; ++iter) {
				auto fit_begin = boost::make_filter_iterator(is_not_nan, iter->begin(), iter->end());
				auto fit_end   = boost::make_filter_iterator(is_not_nan, iter->end(), iter->end());

				variances.emplace_back(computeVariances(fit_begin, fit_end));
				vars.emplace_back(variances.back().var);
			}

			return std::make_pair(variances, statistic::median<value_type>(vars.begin(), vars.end()));
		}

		/**
		 * Computes the pooling factor needed for the shrinking step.
		 *
		 * @param genes Vector of genes
		 * @param var_median Median of all empirical variances.
		 * @return Pooling factor
		 */
		value_type estimatePoolingFactor(const std::vector<Gene<value_type>>& genes, value_type var_median){
			value_type sum_Var = (value_type) 0;
			value_type sum = (value_type) 0;
			for(Gene<value_type> g : genes){
				sum_Var += g.Var;
				sum += std::pow(g.var - var_median, 2);
			}
			return std::min((value_type)1, sum_Var / sum);
		}

		/**
		 * Computes all pooled variances.
		 *
		 * @param[in,out] genes Vector of genes
		 * @param[in] var_median Median of all empirical variances
		 */
		void computePooledVariances(std::vector<Gene<value_type>>& genes, value_type var_median){
			value_type lambda = estimatePoolingFactor(genes,var_median);
			value_type lambda_comp = ((value_type)1) - lambda;
			value_type tmp = lambda * var_median;
			for(Gene<value_type>& g : genes){
				g.shrink_var = tmp + lambda_comp * g.var;
			}
		}

		private:
			std::vector<value_type> wik_;
	};
}

#endif // GT2_CORE_SHRINKAGE_T_TEST_H

