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

#ifndef GT2_REGGAE_ANALYSIS_H
#define GT2_REGGAE_ANALYSIS_H

#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/Exception.h>
#include <genetrail2/core/macros.h>
#include <genetrail2/core/Statistic.h>

#include "RegulatorEffectResult.h"

namespace GeneTrail
{
template <typename ValueType> class GT2_EXPORT REGGAEAnalysis
{
  public:
	using value_type = ValueType;

	REGGAEAnalysis(std::vector<size_t>& sorted_targets,
	               RegulationFile<value_type>& regulationFile,
	               std::vector<std::vector<value_type>>& regulator2associations;
	               bool use_absolute_value, )
	    : sorted_targets_(sorted_targets),
		  regulationFile_(regulationFile),
	      regulator2associations_(regulator2associations),
	      use_absolute_value_(use_absolute_value),
		  max_number_of_regulators_(0),
		  biggest_regulator_idx_(0)
	{
	}

  private:
	/**
	 * Creates a sorted list of regulators REA algorithm.
	 */
	template <typename RegulatorImpactScore>
	void preprocessing_(RegulatorImpactScore impactScore)
	{
		std::vector<size_t> tf_list_tmp;
		size_t size_of_tf_list = 0;
		for(size_t targetname : sorted_targets_) {
			// Check if gene is targetted by any Regulator
			if(!regulationFile_.checkTarget(targetname)) {
				continue;
			}

			auto& targets = regulationFile_.target2regulations(targetname);

			// Get the number of regulators
			size_of_tf_list += targets.size();

			max_number_of_regulators =
			    std::max(max_number_of_regulators, targets.size());
		}

		tf_list_tmp.reserve(size_of_tf_list);

		for(size_t i = 0; i < max_number_of_regulators; ++i) {
			for(size_t targetname : sorted_targets_) {
				// Check if gene is targetted by any Regulator
				if(!regulationFile_.checkTarget(targetname)) {
					continue;
				}

				std::vector<Regulation>& regulators = regulationFile_.target2regulations(targetname);

				if(regulators.size() <= i) {
					continue;
				}

				size_t regulator_i = std::get<0>(regulators[i]);

				biggest_regulator_idx =
				    std::max(regulator_i, biggest_regulator_idx);

				tf_list_tmp.emplace_back(regulator_i);
			}
		}

		tf_list = tf_list_tmp;
	}

	std::vector<size_t> sorted_targets_;
	RegulationFile<value_type>& regulationFile_;
	std::vector<std::vector<value_type>> regulator2associations_;
	bool use_absolute_value_;
	
	
	size_t max_number_of_regulators_;
	size_t biggest_regulator_idx_;
	
	std::vector<std::vector<size_t>> regulator_inidces_;
	std::vector<RegulatorEffectResult> results_;
	std::vector<size_t> tf_list;
};
}

#endif // GT2_REGULATOR_REGGAE_ANALYSIS_H
