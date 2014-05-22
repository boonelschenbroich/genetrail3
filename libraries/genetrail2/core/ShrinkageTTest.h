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
#ifndef GT2_CORE_SHRINKAGE_T_TEST_H
#define GT2_CORE_SHRINKAGE_T_TEST_H

#include <vector>
#include <cmath>
#include <tuple>
#include <utility>

#include "config.h"

#include "Statistic.h"

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
     */
    template <typename value_type, typename InputIterator1,typename InputIterator2>
	class GT2_EXPORT ShrinkageTTest {

		public:

		/**
		 * Computes the different variances for a single gene.
		 *
		 * @param begin InputIterator
		 * @param end InputIterator
		 * @return Gene struct holding all needed fields
		 */
    	template <typename InputIterator>
      	Gene<value_type> computeVariances(InputIterator begin, InputIterator end) {
			Gene<value_type> g;

			value_type mean = statistic::mean<value_type, InputIterator>(begin, end);
        	g.mean = mean;
			size_t size = std::distance(begin, end);
			g.size = size;

			std::vector<value_type> wik;
        	for (auto iter = begin; iter != end; ++iter) {
          		wik.emplace_back(std::pow(*iter - mean, 2));
        	}

			auto sum_wik = std::accumulate(wik.begin(), wik.end(),  0.0);
			auto wk = sum_wik / ((value_type)size);
			value_type vk = sum_wik / ((value_type)(size - 1));
			g.var = vk;

			value_type Var = (value_type) 0;
			for (auto iter = wik.begin(); iter != wik.end(); ++iter) {
				Var += std::pow(*iter - wk, 2);
			}
			Var *= ((value_type)size) / std::pow((value_type)(size -1),3);
			g.Var = Var;

			return g;
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
			std::vector<Gene<value_type>> variances;
			std::vector<value_type> vars;
			for (auto iter = begin; iter != end; ++iter) {
				variances.emplace_back(computeVariances(iter->begin(), iter->end()));
				vars.emplace_back(variances.back().var);
			}
			return std::make_pair(variances, statistic::median<value_type, decltype(vars.begin())>(vars.begin(), vars.end()));
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
				sum += std::pow( g.var - var_median ,2);
			}
			return std::min((value_type)1, sum_Var / sum);
		}

		/**
		 * Computes all pooled variances.
		 *
		 * @param genes Vector of genes
		 * @param var_median Median of all empirical variances
		 * @return Vector of pooled variances
		 */
		void computePooledVariances(std::vector<Gene<value_type>>& genes, value_type var_median){
			value_type lambda = estimatePoolingFactor(genes,var_median);
			value_type lambda_comp = ((value_type)1) - lambda;
			value_type tmp = lambda * var_median;
			for(Gene<value_type>& g : genes){
				g.shrink_var = tmp + lambda_comp * g.var;
			}
		}
	};
}

#endif // GT2_CORE_SHRINKAGE_T_TEST_H

