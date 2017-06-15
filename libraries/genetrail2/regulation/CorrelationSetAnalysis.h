/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2016 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_CORRELATION_SET_ANALYSIS_H
#define GT2_CORRELATION_SET_ANALYSIS_H

#include <boost/numeric/conversion/cast.hpp>

#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseColumnSubset.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/Exception.h>
#include <genetrail2/core/HTest.h>
#include <genetrail2/core/macros.h>
#include <genetrail2/core/MatrixIterator.h>
#include <genetrail2/core/NameDatabases.h>
#include <genetrail2/core/Statistic.h>

#include <genetrail2/regulation/RegulationFile.h>
#include <genetrail2/regulation/RegulatorEffectResult.h>
#include <genetrail2/regulation/RegulatoryImpactFactors.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <tuple>
#include <random>
#include <chrono>
#include <set>

namespace GeneTrail
{
template<typename ValueType>
class GT2_EXPORT CorrelationSetAnalysis
{
  public:
	using value_type = ValueType;
	using Regulation = std::tuple<size_t, size_t, value_type>;

	CorrelationSetAnalysis(	DenseMatrix* matrix,
							MatrixNameDatabase& name_database,
							std::vector<size_t>& sorted_targets,
							RegulationFile<value_type>& regulationFile)
	    : matrix_(matrix),
		  correlations_(matrix->rows(), matrix->rows()),
		  name_database_(name_database),
		  sorted_targets_(sorted_targets),
	      regulationFile_(regulationFile),
		  regulators_(regulationFile.regulators()),
		  max_regulator_(*std::max_element(regulators_.begin(), regulators_.end()))
	{
		results_.resize(max_regulator_ + 1);
		init_();
	}

	template <typename CorrelationCoefficinent> 
	std::vector<RegulatorEffectResult> run(CorrelationCoefficinent func, size_t seed, size_t runs)
	{
		std::cout << "INFO: Calculating correlations" << std::endl;
		create_correlations(func);

		std::cout << "INFO: Calculating scores" << std::endl;
		for(size_t regulator : regulators_) {
			compute_score_(regulator);
		}

		std::cout << "INFO: Calculating p-values" << std::endl;
		std::mt19937 twister(seed);
		std::uniform_int_distribution<int> distribution(0, matrix_->rows()-1);
		for(size_t i=0; i<runs; ++i){
			std::cout << "INFO: Performing permutation " << i << std::endl;
			perform_permutation_(twister, distribution);
		}

		std::vector<RegulatorEffectResult> results;
		results.reserve(regulators_.size());
		for(auto& regulator : regulators_){
			results_[regulator].p_value = ((double)(results_[regulator].number_of_extremer_scores + 1)) / ((double)runs);
			results.emplace_back(results_[regulator]);
		}
		return results_;
	}

  private:
	void init_() {
		std::set<size_t> number_of_targets;
		for(size_t regulator : regulators_) {
			number_of_targets.emplace(regulationFile_.regulator2regulations(regulator).size());
		}
		std::copy(number_of_targets.begin(), number_of_targets.end(), std::back_inserter(number_of_targets_));
		std::sort(number_of_targets_.begin(), number_of_targets_.end());
		max_number_of_targets_ = *std::max_element(number_of_targets_.begin(), number_of_targets_.end());
		number_of_targets_to_regulator_.resize(max_number_of_targets_ + 1);
		for(size_t regulator : regulators_) {
			number_of_targets_to_regulator_[regulationFile_.regulator2regulations(regulator).size()].emplace_back(regulator);
		}
	}

	template <typename CorrelationCoefficient>
	void create_correlations(CorrelationCoefficient func){
		for (size_t i=0; i<matrix_->rows(); ++i) {
			correlations_.set(i, i, 1.0);
			for (size_t j=i+1; j<matrix_->rows(); ++j) {
				RowMajorMatrixIterator<Matrix> regulator_it(matrix_, i),
			    							   target_it(matrix_, j);
				auto corr = func.compute(regulator_it->begin(), regulator_it->end(),
			                  			 target_it->begin(), target_it->end());
				correlations_.set(i, j, corr);
				correlations_.set(j, i, corr);
			}
		}
	}

	void compute_score_(size_t regulator){
		value_type sum = 0.0;
		auto regulations = regulationFile_.regulator2regulations(regulator);
		size_t size = regulations.size();
		for (size_t i=0; i < size; ++i){
			for (size_t j=i+1; j < size; ++j){
				auto t1 = std::get<1>(regulations[i]);
				auto t2 = std::get<1>(regulations[j]);
				sum += correlations_(t1, t2);
			}
		}
		value_type n = (value_type)size;
		results_[regulator].name = name_database_(regulator);
		results_[regulator].hits = size;
		results_[regulator].score = (2.0 / (n*(n-1.0))) * sum;
	}

	void perform_permutation_(std::mt19937 twister,	std::uniform_int_distribution<int> distribution){
		auto next = number_of_targets_.begin();
		double sum = 0.0;
		std::vector<size_t> regs;
		regs.reserve(max_number_of_targets_);
		for (size_t i=0; i <= max_number_of_targets_; ++i){
			regs.emplace_back(distribution(twister));
		}
		for (size_t i=0; i <= max_number_of_targets_; ++i){
			for (size_t j=i+1; j <= max_number_of_targets_; ++j){
				sum += correlations_(regs[i], regs[j]);
			}
			if(i == (*next)) {
				double n = (double)i;
				double mean = (2.0 / (n*(n-1.0))) * sum;
				for (size_t regulator : number_of_targets_to_regulator_[*next]){
					if(mean > results_[regulator].score) {
						results_[regulator].number_of_extremer_scores += 1;
					}
				}
				++next;
			}
		}
	}

  	DenseMatrix* matrix_;
	DenseMatrix correlations_;
	MatrixNameDatabase& name_database_;
	std::vector<size_t>& sorted_targets_;
	RegulationFile<double>& regulationFile_;

	std::vector<size_t> regulators_;
	size_t max_regulator_;
	std::vector<RegulatorEffectResult> results_;

	std::vector<size_t> number_of_targets_;
	size_t max_number_of_targets_;
	std::vector<std::vector<size_t>> number_of_targets_to_regulator_;
};
}

#endif // GT2_CORRELATION_SET_ANALYSIS_H
