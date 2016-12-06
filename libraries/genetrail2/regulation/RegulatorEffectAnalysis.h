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

#ifndef GT2_REGULATOR_EFFECT_ANALYSIS_H
#define GT2_REGULATOR_EFFECT_ANALYSIS_H

#include <boost/numeric/conversion/cast.hpp>

#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseColumnSubset.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/Exception.h>
#include <genetrail2/core/HTest.h>
#include <genetrail2/core/macros.h>
#include <genetrail2/core/MatrixIterator.h>
#include <genetrail2/core/Statistic.h>

#include <genetrail2/regulation/RegulationFile.h>
#include <genetrail2/regulation/RegulatorEffectResult.h>
#include <genetrail2/regulation/RegulatoryImpactFactors.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <tuple>

namespace GeneTrail
{
class GT2_EXPORT RegulatorEffectAnalysis
{
	using Regulation = std::tuple<size_t, size_t, double>;

  public:
	RegulatorEffectAnalysis(DenseColumnSubset* reference,
	                        DenseColumnSubset* sample,
	                        RegulationFile<double>& regulationFile)
	    : reference_(reference),
	      sample_(sample),
	      regulationFile_(regulationFile)
	{
	}

	template <typename RIF> std::vector<RegulatorEffectResult> run(RIF func, bool standardize)
	{
		std::vector<RegulatorEffectResult> results;
		for(size_t regulator : regulationFile_.regulators()) {
			auto result = run_regulator_(regulator, func);
			if(result.name != ""){
				results.emplace_back(std::move(result));	
			}
		}
		if(standardize) {
			compute_p_value(results);
		}
		return results;
	}

  private:
	template <typename RIF>
	double rif(DenseColumnSubset* reference, DenseColumnSubset* sample,
	           const std::vector<Regulation>& regulations, RIF func)
	{
		size_t regulator_idx = std::get<0>(regulations[0]);
		std::cout << "Processing: " << reference_->rowNames()[regulator_idx] << std::endl;
		RowMajorMatrixIterator<DenseColumnSubset> reference_regulator_it(reference,
		                                                      regulator_idx);
		RowMajorMatrixIterator<DenseColumnSubset> sample_regulator_it(sample,
		                                                   regulator_idx);
		double score = 0.0;
		for(const auto& reg : regulations) {
			size_t target_idx = std::get<1>(reg);
			RowMajorMatrixIterator<DenseColumnSubset> reference_target_it(reference,
			                                                   target_idx);
			RowMajorMatrixIterator<DenseColumnSubset> sample_target_it(sample,
			                                                target_idx);
			score += func.compute(reference_regulator_it, reference_target_it,
			              sample_regulator_it, sample_target_it);
		}
		return score / regulations.size();
	}

	template <typename RIF>
	RegulatorEffectResult run_regulator_(size_t regulator, RIF func)
	{
		if (!regulationFile_.checkRegulator(regulator)){
			return RegulatorEffectResult();
		}
		
		auto regulations = regulationFile_.regulator2regulations(regulator);
		RegulatorEffectResult result;
		
		if(regulations.size() > 0){
			result.name = reference_->rowNames()[regulator];
			result.score = rif(reference_, sample_, regulations, func);
			result.hits = regulations.size();
		}
		
		return std::move(result);
	}

	void compute_p_value(std::vector<RegulatorEffectResult>& results)
	{
		std::vector<double> scores;
		scores.reserve(results.size());
		for(auto& result : results) {
			scores.emplace_back(result.score);
		}
		double mean = 0.0, var = 0.0, sd = 0.0;
		size_t size;
		std::tie(mean, var, size) =
		    statistic::mean_var_size<double>(scores.begin(), scores.end());
		sd = sqrt(var);
		for(auto& result : results) {
			result.score = (result.score - mean) / sd;

			HTest::ZTest test(result.score);
			if(result.score < 0) {
				result.p_value = HTest::lowerTailedPValue(test, result.score);
			} else {
				result.p_value = HTest::upperTailedPValue(test, result.score);
			}
			result.corrected_p_value = result.p_value;
		}
	}

	DenseColumnSubset* reference_;
	DenseColumnSubset* sample_;
	RegulationFile<double>& regulationFile_;
};
}

#endif // GT2_REGULATOR_EFFECT_ANALYSIS_H
