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

#include "GroupedScores.h"
#include "MatrixTools.h"

using namespace GeneTrail;


void GroupedScores::calculateGroupedScores(
	DenseMatrix& matrix,
	const std::vector<Metadata>& meta,
	const std::string method,
	DenseMatrix& result
) const{
	std::map<std::string, std::vector<unsigned int>> group_indices;
	ORAGroupPreference::parseGroups(matrix, meta, group_indices);
	std::vector<std::string> groups;
	for(const auto& entry: group_indices) groups.push_back(entry.first);
	
	result = DenseMatrix(matrix.rows(), group_indices.size());
	result.setColNames(groups);
	result.setRowNames(matrix.rowNames());
	MatrixHTest htest;
	MatrixTools mtools;
	
	size_t group_idx=-1;
	for(const auto& current_group: groups){
		group_idx++;
		Samples g1, g2;
		createSamples(matrix, group_indices, current_group, g1, g2);
		auto subset = mtools.splitMatrix(matrix, g1, g2);
		auto gene_set = htest.test(method, std::get<0>(subset), std::get<1>(subset));
		addToResult(result, gene_set, group_idx);
	}
	return;
}

void GroupedScores::createSamples(
	const DenseMatrix& matrix,
	const std::map<std::string, std::vector<unsigned int>>& group_indices,
	const std::string& current_group,
	Samples& g1, Samples& g2
) const{
	auto col_names = matrix.colNames();
	auto entry = group_indices.find(current_group);
	
	// Here, we use that the vector of indices is sorted
	unsigned int current_lookup = 0;
	for(unsigned int i=0; i < col_names.size(); i++){
		const auto& col_name = col_names[i];
		if(current_lookup < entry->second.size() && entry->second[current_lookup] == i){
			g1.push_back(col_name);
			current_lookup++;
		} else{
			g2.push_back(col_name);
		}
	}
	return;
}

void GroupedScores::addToResult(DenseMatrix& result, const Scores& gene_set, size_t column_idx) const{
	size_t row_idx = -1;
	for(const auto& score : gene_set.scores()) {
		row_idx++;
		result(row_idx, column_idx) = score;
	}
}


