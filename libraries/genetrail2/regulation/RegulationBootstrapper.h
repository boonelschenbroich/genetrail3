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

#ifndef GT2_REGULATION_REGULATION_BOOTSTRAPPER_H
#define GT2_REGULATION_REGULATION_BOOTSTRAPPER_H

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

namespace GeneTrail
{
template <typename ValueType> class GT2_EXPORT RegulationBootstrapper
{
  public:
	using value_type = ValueType;
	using Regulation = std::tuple<size_t, size_t, value_type>;
	using dist_type = std::uniform_int_distribution<size_t>;

	RegulationBootstrapper(DenseMatrix* matrix, unsigned seed)
	    : matrix_(matrix),
	      bootstrap_sample_(matrix_->cols()),
	      subset_(matrix, bootstrap_sample_.begin(), bootstrap_sample_.end()),
	      generator_(seed),
	      distribution_(0, matrix_->cols() - 1)
	{
		// Use entire matrix
		std::iota(bootstrap_sample_.begin(), bootstrap_sample_.end(), 0u);
	}

	/**
	 * This method creates a random bootstrap sample (with replacements)
	 * that can then be used as indices for the matrix.
	 */
	void create_bootstrap_sample()
	{
		for(size_t i = 0; i < matrix_->cols(); ++i) {
			bootstrap_sample_[i] = distribution_(generator_);
		}
		std::sort(bootstrap_sample_.begin(), bootstrap_sample_.end(),
		          [](const size_t& a, const size_t& b) { return a < b; });
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
	void perform_bootstrapping_run(std::vector<Regulation>& regulations,
	                               bool use_absolute_values,
	                               RegulatorImpactScore score)
	{
		subset_.assign(bootstrap_sample_.begin(), bootstrap_sample_.end());

		size_t target_idx = std::get<1>(regulations[0]);
		for(Regulation& r : regulations) {
			size_t regulator_idx = std::get<0>(r);

			RowMajorMatrixIterator<Matrix> regulator_it(&subset_,
			                                            regulator_idx),
			    target_it(&subset_, target_idx);

			value_type result =
			    score.compute(regulator_it->begin(), regulator_it->end(),
			                  target_it->begin(), target_it->end());

			std::get<2>(r) = result;
		}

		sort_(regulations, use_absolute_values);
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
	void sort_(std::vector<Regulation>& regulations, bool use_absolute_values)
	{
		// Sort values decreasingly
		if(use_absolute_values) {
			std::sort(regulations.begin(), regulations.end(),
			          [](const Regulation& a, const Regulation& b) {
				return std::abs(std::get<2>(a)) > std::abs(std::get<2>(b));
			});
		} else {
			std::sort(regulations.begin(), regulations.end(),
			          [](const Regulation& a, const Regulation& b) {
				return std::get<2>(a) > std::get<2>(b);
			});
		}
	}

	DenseMatrix* matrix_;
	std::vector<size_t> bootstrap_sample_;
	DenseColumnSubset subset_;
	size_t samples_;
	std::mt19937 generator_;
	dist_type distribution_;
};
}

#endif // GT2_REGULATION_REGULATION_BOOTSTRAPPER_H
