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


void ORAGroupPreference::calculatePreference(
	const DenseMatrix& matrix,
	const std::vector<Metadata>& metadata,
	DenseMatrix& result
) const{
	std::map<std::string, std::vector<unsigned int>> group_indices;
	ORAGroupPreference::parseGroups(matrix, metadata, group_indices);
	std::vector<std::string> groups;
	for(const auto& entry: group_indices) groups.push_back(entry.first);
	
	result = DenseMatrix(matrix.rows(), group_indices.size());
	result.setColNames(groups);
	result.setRowNames(matrix.rowNames());
	
	for(size_t row_index=0; row_index < matrix.rows(); ++row_index){
		size_t idx_non_empty = -1;
		for(const auto& current_entry: group_indices){
			if(current_entry.second.empty()) continue;
			idx_non_empty++;
			std::vector<size_t> table(4, 0);
			for(const auto& other_entry: group_indices){
				for(const auto column_index: other_entry.second){
					addToTable(table, matrix(row_index, column_index),
							   current_entry.first == other_entry.first);
				}
			}
			result(row_index, idx_non_empty) = computePValue_(table);
		}
	}
}

void ORAGroupPreference::parseGroups(
	const DenseMatrix& matrix,
	const std::vector<Metadata>& metadata,
	std::map<std::string, std::vector<unsigned int>>& group_indices)
{
	auto col_names = matrix.colNames();
	size_t column_index = -1;
	
	for(const std::string& col_name: col_names){
		column_index++;
		std::string groupName = "";
		for(const auto& meta: metadata){
			if(meta.has(col_name)){
				std::string group = get<std::string>(meta.get(col_name));
				groupName = groupName == "" ? group : groupName + "_" + group;
			} else{
				throw IOError("The metadata file had no entry for the sample " + col_name + ".");
			}
		}
		auto entry = group_indices.find(groupName);
		if(entry == group_indices.end()){
			group_indices.emplace(groupName, std::vector<unsigned int>(1, column_index));
		} else{
			entry->second.push_back(column_index);
		}
	}
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
	// Ensures that E != 0
	if(l == 0) return -1.0;
	if(l == m) return 1.0;
	if(n == 0) return 1.0;
	if(n == m) return 1.0;
	
	// Calculate E and chi^2
	auto p_11 = (l*n)/((double)m*m), p_12 = ((m-l)*n)/((double)m*m);
	auto p_21 = (l*(m-n))/((double)m*m), p_22 = ((m-l)*(m-n))/((double)m*m);
	auto E_11 = m*p_11, E_12 = m*p_12;
	auto E_21 = m*p_21, E_22 = m*p_22;
	auto o_11 = (double)k,   o_12 = (double)n-k;
	auto o_21 = (double)l-k, o_22 = (m-l)-o_12;
	auto chi = (o_11-E_11)*(o_11-E_11)/E_11 +
			   (o_12-E_12)*(o_12-E_12)/E_12 +
			   (o_21-E_21)*(o_21-E_21)/E_21 +
			   (o_22-E_22)*(o_22-E_22)/E_22;
	
	// Calculate p-value from chi^2
	double res = boost::math::cdf(complement(test_, chi));
	if(o_11 < E_11){
		res = -1*res;
	}
	return res;
}


