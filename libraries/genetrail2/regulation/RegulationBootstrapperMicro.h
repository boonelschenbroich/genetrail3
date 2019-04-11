
/**
* GeneTrail2 - An efficent library for interpreting genetic data
* Copyright (C) 2016 Tim Kehl tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_REGULATION_REGULATION_BOOTSTRAPPER_MICRO_H
#define GT2_REGULATION_REGULATION_BOOTSTRAPPER_MICRO_H

#include <genetrail2/core/DenseColumnSubset.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/Matrix.h>
#include <genetrail2/core/MatrixIterator.h>
#include <genetrail2/core/Statistic.h>

#include <iostream>
#include <random>
#include <vector>
#include <tuple>
#include <cmath>

#include "RegulationFile.h"

namespace GeneTrail
{
template <typename ValueType> class GT2_EXPORT RegulationBootstrapperMicro
{
  public:
	using value_type = ValueType;
	using Regulation = std::tuple<size_t, size_t, value_type>;
	using dist_type = std::uniform_int_distribution<size_t>;

	RegulationBootstrapperMicro(DenseMatrix* matrix, DenseMatrix* matrix_micro, unsigned seed, size_t firstControl)
	    : matrix_(matrix),
	      matrix_micro_(matrix_micro),
	      bootstrap_sample_matrix_(matrix_->cols()),
	      subset_matrix_(matrix, bootstrap_sample_matrix_.begin(), bootstrap_sample_matrix_.end()),
	      subset_micro_(matrix_micro, bootstrap_sample_matrix_.begin(), bootstrap_sample_matrix_.end()),
	      generator_(seed),
	      distribution_dis_(0, firstControl-1),
	      distribution_con_(firstControl,matrix_->cols() - 1),
	      firstControl_(firstControl)
	     
	{
	      std::iota(bootstrap_sample_matrix_.begin(), bootstrap_sample_matrix_.end(), 0u);
	}
	

	/**
	 * This method creates a random bootstrap sample (with replacements)
	 * that can then be used as indices for the matrix.
	 */
	void create_bootstrap_sample()
	{		
		for(size_t i = 0; i < firstControl_; ++i) {
		      bootstrap_sample_matrix_[i] = distribution_dis_(generator_);
		}
		for(size_t x = firstControl_; x < matrix_->cols();++x){
			bootstrap_sample_matrix_[x] = distribution_con_(generator_);
		}
		std::sort(bootstrap_sample_matrix_.begin(), bootstrap_sample_matrix_.end(),
		          [](const size_t& a, const size_t& b) { return a < b; });
		          
	}
	void set_bootstrap_sample(std::vector<size_t>& vec){
	  bootstrap_sample_matrix_ = vec;
	}

	/**
	 * This method creates a matrix subset using the current bootstrap sample
	 * and recomputes the correlations between all regulator-target pairs.
	 *
	 * @param regulations Vector of regulations (std::tuple<size_t, size_t,
	 *value_type>)
	 * @param use_absolute_values Flag indicating if correlations should be
	 *sorted absolute or relative
	 */
	template <typename RegulatorImpactScore>
	void perform_bootstrapping_run(RegulationFile<double>& rfile,
								   std::vector<Regulation>& regulations,
								   bool normalize_scores,
	                               bool use_absolute_values,
								   bool sort_decreasingly,
	                               RegulatorImpactScore score)
	{
		
		subset_matrix_.assign(bootstrap_sample_matrix_.begin(), bootstrap_sample_matrix_.end());
		subset_micro_.assign(bootstrap_sample_matrix_.begin(), bootstrap_sample_matrix_.end());
		
		size_t target_idx = std::get<1>(regulations[0]);
		for(Regulation& r : regulations) {

			size_t regulator_idx = std::get<0>(r);
			RowMajorMatrixIterator<Matrix> regulator_it(&subset_micro_,
			                                            regulator_idx),
			    target_it(&subset_matrix_, target_idx);
			
			value_type result =
			    score.compute(regulator_it->begin(), regulator_it->end(),
			                  target_it->begin(), target_it->end());
			    
			if(normalize_scores) {
				size_t tmp = rfile.regulator2regulations(regulator_idx).size();
				std::get<2>(r) = result / tmp;
			} else {
				std::get<2>(r) = result;
			}
		}

		sort_(regulations, use_absolute_values, sort_decreasingly);

	}

  private:
	/**
	 * This method sorts the given regutions vector.
	 *
	 * @param regulations Vector of regulations (std::tuple<size_t, size_t,
	 *value_type>)
	 * @param use_absolute_values Flag indicating if correlations should be
	 *sorted absolute or relative
	 */
	void sort_(std::vector<Regulation>& regulations, bool use_absolute_values, bool sort_decreasingly)
	{
		// Sort values decreasingly
		if(use_absolute_values) {
			if (sort_decreasingly) {
				std::sort(regulations.begin(), regulations.end(),
						[](const Regulation& a, const Regulation& b) {
					return std::abs(std::get<2>(a)) > std::abs(std::get<2>(b));
				});
			} else {
				std::sort(regulations.begin(), regulations.end(),
						[](const Regulation& a, const Regulation& b) {
					return std::abs(std::get<2>(a)) < std::abs(std::get<2>(b));
				});
			}
		} else {
			if (sort_decreasingly) {
				std::sort(regulations.begin(), regulations.end(),
						[](const Regulation& a, const Regulation& b) {
					return std::get<2>(a) > std::get<2>(b);
				});
			} else {
				std::sort(regulations.begin(), regulations.end(),
						[](const Regulation& a, const Regulation& b) {
					return std::get<2>(a) < std::get<2>(b);
				});
			}
		}
	}

	DenseMatrix* matrix_;
	DenseMatrix* matrix_micro_;
	std::vector<size_t> bootstrap_sample_matrix_;
	DenseColumnSubset subset_matrix_;
	DenseColumnSubset subset_micro_;
	size_t samples_;
	std::mt19937 generator_;
	dist_type distribution_dis_;
	dist_type distribution_con_;
	size_t firstControl_;

};
}

#endif // GT2_REGULATION_REGULATION_BOOTSTRAPPER_H
