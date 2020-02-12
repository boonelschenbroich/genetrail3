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
	
	result = DenseMatrix(matrix.rows()+1, group_indices.size());
	result.setColNames(groups);
	std::vector<std::string> row_names(matrix.rowNames());
	row_names.insert(row_names.begin(), "GroupSizes");
	result.setRowNames(row_names);
	addGroupSizeToResult(result, group_indices);
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

void GroupedScores::addGroupSizeToResult(
	DenseMatrix& result,
	const std::map<std::string, std::vector<unsigned int>>& group_indices
) const{
	size_t col_idx = -1;
	for(const auto& indices: group_indices) {
		col_idx++;
		result(0, col_idx) = indices.second.size();
	}
}

void GroupedScores::addToResult(DenseMatrix& result, const Scores& gene_set, size_t column_idx) const{
	size_t row_idx = 0;  // First row is already filled
	for(const auto& score : gene_set.scores()) {
		row_idx++;
		result(row_idx, column_idx) = score;
	}
}

void GroupedScores::calculateGroupedScores(
	DenseMatrix& matrix,
	const std::vector<std::string>& sampleGroups,
	const std::vector<std::string>& referenceGroups,
	const std::string&,
	DenseMatrix& result
) const{
	result = DenseMatrix(matrix.rows()-1, 1);
	std::vector<std::string> row_names(matrix.rowNames());
	row_names.erase(row_names.begin());
	result.setRowNames(row_names);
	
	MatrixHTest htest;
	MatrixTools mtools;
	
	auto subset = mtools.splitMatrix(matrix, sampleGroups, referenceGroups);
	auto sample = std::get<0>(subset);
	auto reference = std::get<1>(subset);
	std::vector<size_t> sampleIndices, referenceIndices;
	getIndices(sample, sampleGroups, sampleIndices);
	getIndices(reference, referenceGroups, referenceIndices);
	std::vector<size_t> sampleSizes, referenceSizes;
	getGroupSizes(sample, sampleIndices, sampleSizes);
	getGroupSizes(reference, referenceIndices, referenceSizes);
	size_t n_sample = getTotalGroupSizes(sampleSizes);
	size_t n_reference = getTotalGroupSizes(referenceSizes);
	
	for(size_t idx_row=1; idx_row < matrix.rows(); idx_row++){
		double sampleScore = getMeanScore(idx_row, sample, sampleIndices, sampleSizes, n_sample);
		double referenceScore = getMeanScore(idx_row, reference, referenceIndices, referenceSizes, n_reference);
		result(idx_row-1, 0) = sampleScore - referenceScore;
	}
	return;
}

void GroupedScores::getIndices(
	const DenseColumnSubset& sub,
	const std::vector<std::string>& groups,
	std::vector<size_t>& indices
) const{
	for(const auto& group: groups){
		indices.push_back(sub.colIndex(group));
	}
}

void GroupedScores::getGroupSizes(
	const DenseColumnSubset& sub,
	const std::vector<size_t>& indices,
	std::vector<size_t>& sizes
) const{
	for(auto idx: indices){
		sizes.push_back(sub(0, idx));
	}
}

size_t GroupedScores::getTotalGroupSizes(
	const std::vector<size_t>& sizes
) const{
	size_t n = 0;
	for(auto size: sizes){
		n += size;
	}
	return n;
}

double GroupedScores::getMeanScore(
	size_t idx_row,
	const DenseColumnSubset& sub,
	const std::vector<size_t>& indices,
	const std::vector<size_t>& sizes,
	size_t n
) const{
	double mean = 0.0;
	for(size_t idx=0; idx < indices.size(); idx++){
		mean += sizes[idx] * sub(idx_row, indices[idx]);
	}
	return mean/n;
}



