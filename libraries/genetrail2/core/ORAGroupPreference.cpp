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
	std::vector<std::string> groups = ORAGroupPreference::getGroups(metadata);
	std::vector<std::vector<unsigned int>> group_indices = ORAGroupPreference::getGroupIndices(matrix, metadata, groups);
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

std::vector<std::string> ORAGroupPreference::getGroups(const Metadata& metadata){
	std::set<std::string> groups;
	for(auto it=metadata.begin(); it != metadata.end(); ++it){
		groups.insert(get<std::string>(it->second));
	}
	return std::vector<std::string>(groups.begin(), groups.end());
}

std::vector<std::vector<unsigned int>> ORAGroupPreference::getGroupIndices(const DenseMatrix& matrix, const Metadata& metadata, const std::vector<std::string> groups){
	return ORAGroupPreference::getGroupIndices_<unsigned int>(matrix, metadata, groups, 0);
}

std::vector<std::vector<std::string>> ORAGroupPreference::getGroupedSamples(const DenseMatrix& matrix, const Metadata& metadata, const std::vector<std::string> groups){
	return ORAGroupPreference::getGroupIndices_<std::string>(matrix, metadata, groups, "");
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


