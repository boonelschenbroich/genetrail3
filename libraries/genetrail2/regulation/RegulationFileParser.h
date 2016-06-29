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

#ifndef GT2_REGULATION_FILE_PARSER_H
#define GT2_REGULATION_FILE_PARSER_H

#include <genetrail2/core/macros.h>
#include <genetrail2/core/Exception.h>
#include <genetrail2/core/Matrix.h>
#include <genetrail2/core/DenseMatrix.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <unordered_set>
#include <map>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>

namespace GeneTrail
{
template <typename ValueType> class GT2_EXPORT RegulationFileParser
{
  public:
	using value_type = ValueType;
	using Regulation = std::tuple<size_t, size_t, value_type>;

	RegulationFileParser(const DenseMatrix* matrix,
	                     const std::unordered_set<size_t> test_set,
	                     const std::string& file, value_type default_value)
	    : regulator_indices_(matrix->rows(), MAX_MATRIX_INDEX),
	      target_indices_(matrix->rows(), MAX_MATRIX_INDEX)

	{
		read_(matrix, test_set, file, default_value);
	}

	std::vector<Regulation>& target2regulations(size_t target)
	{
		return target2regulations_[target_indices_[target]];
	}

	bool checkTarget(size_t target)
	{
		return target < target_indices_.size() &&
		       target_indices_[target] < MAX_MATRIX_INDEX &&
		       target2regulations_[target_indices_[target]].size() > 0;
	}

	std::vector<Regulation>& regulator2regulations(size_t regulator)
	{
		return regulator2regulations_[regulator_indices_[regulator]];
	}

	bool checkRegulator(size_t regulator)
	{
		return regulator < regulator_indices_.size() &&
		       regulator_indices_[regulator] < MAX_MATRIX_INDEX &&
		       regulator2regulations_[regulator_indices_[regulator]].size() > 0;
	}

  private:
	static constexpr Matrix::index_type MAX_MATRIX_INDEX = std::numeric_limits<Matrix::index_type>::max();
	
	void read_(const DenseMatrix* matrix, const std::unordered_set<size_t>& test_set, const std::string& file,
	           value_type default_value)
	{
		std::ifstream input(file);
		if(!input) {
			throw GeneTrail::IOError("File (" + file + ") is not open for reading");
		}

		std::vector<std::string> sline(2);
		for(std::string line; getline(input, line);) {
			boost::split(sline, line, boost::is_any_of(" \t"));
			if(sline.size() == 2) {
				addRegulation_(matrix, test_set, sline[0], sline[1], default_value);
			} else {
				throw GeneTrail::IOError("Wrong file format.");
			}
		}
	}

	void addRegulation_(const DenseMatrix* matrix, const std::unordered_set<size_t>& test_set, const std::string& regulator,
	                    const std::string& target, value_type default_value)
	{

		auto regulator_idx = matrix->rowIndex(regulator);
		auto target_idx = matrix->rowIndex(target);
		
		if(test_set.find(target_idx) == test_set.end()) {
			return;
		}

		if(regulator_idx == MAX_MATRIX_INDEX || target_idx == MAX_MATRIX_INDEX) {
			return;
		}

		if(regulator_indices_[regulator_idx] == MAX_MATRIX_INDEX) {
			regulator_indices_[regulator_idx] = regulator2regulations_.size();
			regulator2regulations_.emplace_back();
		}

		if(target_indices_[target_idx] == MAX_MATRIX_INDEX) {
			target_indices_[target_idx] = target2regulations_.size();
			target2regulations_.emplace_back();
		}

		const Regulation reg = std::make_tuple(regulator_idx, target_idx, default_value);
		regulator2regulations_[regulator_indices_[regulator_idx]].emplace_back(reg);
		target2regulations_[target_indices_[target_idx]].emplace_back(reg);
	}

	std::vector<size_t> regulator_indices_;
	std::vector<size_t> target_indices_;
	
	std::vector<std::vector<Regulation>> target2regulations_;
	std::vector<std::vector<Regulation>> regulator2regulations_;
};
}

#endif // GT2_REGULATION_FILE_PARSER_H
