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

#ifndef GT2_REGULATION_FILE_H
#define GT2_REGULATION_FILE_H

#include <genetrail2/core/macros.h>

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <tuple>

namespace GeneTrail
{
template <typename ValueType> class GT2_EXPORT RegulationFile
{
  public:
	using value_type = ValueType;
	using Regulation = std::tuple<size_t, size_t, value_type>;

	RegulationFile() = delete;

	RegulationFile(size_t number_of_genes, size_t max_index)
	    : max_index_(max_index),
		  regulator_indices_(number_of_genes, max_index_),
	      target_indices_(number_of_genes, max_index_),
		  total_number_of_targets_(number_of_genes, max_index_)
	{
	}
	
	RegulationFile(size_t number_of_genes, size_t number_of_mirnas, size_t max_index)
	    : max_index_(max_index),
		  regulator_indices_(number_of_mirnas, max_index_),
	      target_indices_(number_of_genes, max_index_),
		  total_number_of_targets_(number_of_mirnas, max_index_)
	{
	}

	std::vector<Regulation>& target2regulations(size_t target)
	{
		return target2regulations_[target_indices_[target]];
	}

	bool checkTarget(size_t target)
	{
		return target < target_indices_.size() &&
		       target_indices_[target] < max_index_ &&
		       target2regulations_[target_indices_[target]].size() > 0;
	}
	
	void addRegulation(size_t regulator_idx, size_t target_idx, value_type value){
		if(regulator_indices_[regulator_idx] == max_index_) {
			regulator_indices_[regulator_idx] = regulator2regulations_.size();
			regulator2regulations_.emplace_back();
		}
		
		if(target_indices_[target_idx] == max_index_) {
			target_indices_[target_idx] = target2regulations_.size();
			target2regulations_.emplace_back();
		}
		
		const Regulation reg = std::make_tuple(regulator_idx, target_idx, value);

		regulator2regulations_[regulator_indices_[regulator_idx]].emplace_back(reg);
		target2regulations_[target_indices_[target_idx]].emplace_back(reg);
	}

	std::vector<Regulation>& regulator2regulations(size_t regulator)
	{
		return regulator2regulations_[regulator_indices_[regulator]];
	}

	bool checkRegulator(size_t regulator)
	{
		return regulator < regulator_indices_.size() &&
		       regulator_indices_[regulator] < max_index_ &&
		       regulator2regulations_[regulator_indices_[regulator]].size() > 0;
	}
	
	std::vector<size_t> regulators(){
		std::vector<size_t> regulators;
		for(size_t i=0; i<regulator_indices_.size(); ++i){ 
			auto r = regulator_indices_[i];
			if(r < max_index_ && regulator2regulations_[r].size() > 0){
				regulators.push_back(i);
			}
		}
		return std::move(regulators);
	}

	size_t maxNumberOfTargets(){
		size_t max = 0;
		for(size_t i=0; i<regulator_indices_.size(); ++i){ 
			auto r = regulator_indices_[i];
			if(r < max_index_){
				max = std::max(regulator2regulations_[r].size(), max);
			}
		}
		return max;
	}
	size_t maxNumberOfRegulators(){
		size_t max = 0;
		for(size_t i=0; i<target_indices_.size(); ++i){ 
			auto r = target_indices_[i];
			if(r < max_index_){
				max = std::max(target2regulations_[r].size(), max);
			}
		}
		return max;
	}

	void increaseNumberOfTargets(size_t regulator) {
		total_number_of_targets_[regulator]+=1;
	}

	size_t getTotalNumberOfTargets(size_t regulator) {
		return total_number_of_targets_[regulator];
	}

  private:
	size_t max_index_;
	std::vector<size_t> regulator_indices_;
	std::vector<size_t> target_indices_;

	std::vector<std::vector<Regulation>> target2regulations_;
	std::vector<std::vector<Regulation>> regulator2regulations_;
	std::vector<size_t> total_number_of_targets_;
};
}

#endif // GT2_REGULATION_FILE_PARSER_H
