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
#include "ORAGroupPreference.h"
#include "Exception.h"

#include <iostream>
#include <set>
#include <vector>
#include <algorithm>

using namespace GeneTrail;


void ORAGroupPreference::calculatePreference(const DenseMatrix& matrix, const Metadata& metadata, DenseMatrix& result) const{
	std::vector<std::string> groups = getGroups(metadata);
	std::vector<std::vector<unsigned int>> group_indices = getGroupIndices(matrix, metadata, groups);
	result = DenseMatrix(matrix.rows(), groups.size());
	result.setColNames(groups);
	result.setRowNames(matrix.rowNames());
	
	for(size_t row_index=0; row_index < matrix.rows(); ++row_index){
		std::cout << (row_index*100.0) / matrix.rows() << "% done. Current: " << matrix.rowNames()[row_index] << std::endl;
		for(size_t current_group=0; current_group < group_indices.size(); ++current_group){
			std::vector<size_t> table(4, 0);
			for(size_t group=0; group < group_indices.size(); ++group){
				for(unsigned int column_index: group_indices[group]){
					addToTable(table, matrix(row_index, column_index), current_group==group);
				}
			}
			result(row_index, current_group) = computePValue_(table);
		}
	}
}

std::vector<std::string> ORAGroupPreference::getGroups(const Metadata& metadata) const{
	std::set<std::string> groups;
	for(auto it=metadata.begin(); it != metadata.end(); ++it){
		groups.insert(get<std::string>(it->second));
	}
	return std::vector<std::string>(groups.begin(), groups.end());
}

std::vector<std::vector<unsigned int>> ORAGroupPreference::getGroupIndices(const DenseMatrix& matrix, const Metadata& metadata, const std::vector<std::string> groups) const{
	std::vector<std::vector<unsigned int>> group_indices(groups.size());
	auto col_names = matrix.colNames();
	size_t column_index = -1;
	
	for(const std::string& col_name: col_names){
		column_index++;
		if(metadata.has(col_name)){
			std::string group = get<std::string>(metadata.get(col_name));
			auto elem_it = std::find(groups.begin(), groups.end(), group);
			if(elem_it == groups.end()){
				throw IOError("Internal server error: element " + group + " was " +
				"inserted into the list of groups, but not found later");
			}
			size_t group_index = std::distance(groups.begin(), elem_it);
			group_indices[group_index].push_back(column_index);
		} else{
			throw IOError("The metadata file had no entry for the sample " + col_name + ".");
		}
	}
	return group_indices;
}

void ORAGroupPreference::addToTable(std::vector<size_t>& table, double p_value, bool is_current_group) const{
	if(is_current_group){
		table[NINDEX]++;
	}
	if(p_value < threshold_){
		table[LINDEX]++;
		if(is_current_group){
			table[KINDEX]++;
		}
	}
	table[MINDEX]++;
}

double ORAGroupPreference::computePValue_(const std::vector<size_t>& table) const{
	return computePValue_(table[MINDEX], table[LINDEX], table[NINDEX], table[KINDEX]);
}

double ORAGroupPreference::computePValue_(size_t m, size_t l, size_t n, size_t k) const {
	auto expected_k = ((double)l * n) / ((double)m);
	bool enriched = expected_k < k;
	big_float p;
	if(enriched) {
		p = hyperTest_.upperTailedPValue(m, l, n, k);
	} else {
		p = hyperTest_.lowerTailedPValue(m, l, n, k);
	}
	return enriched ? p.convert_to<double>() : -1.0 * p.convert_to<double>();
}


