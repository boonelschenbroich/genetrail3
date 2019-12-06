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


void GroupedScores::calculateGroupedScores(DenseMatrix& matrix, const Metadata& metadata,
							const std::string method, DenseMatrix& result) const{
	auto groups = ORAGroupPreference::getGroups(metadata);
	auto group_samples = ORAGroupPreference::getGroupedSamples(matrix, metadata, groups);
	result = DenseMatrix(matrix.rows(), groups.size());
	result.setColNames(groups);
	result.setRowNames(matrix.rowNames());
	MatrixHTest htest;
	MatrixTools mtools;
	
	for(size_t group_idx=0; group_idx < groups.size(); group_idx++){
		auto g1 = getNames(group_samples, group_idx);
		auto g2 = getOtherNames(group_samples, group_idx);
		auto subset = mtools.splitMatrix(matrix, g1, g2);
		auto gene_set = htest.test(method, std::get<0>(subset), std::get<1>(subset));
		addToResult(result, gene_set, group_idx);
	}
}

auto GroupedScores::getNames(const std::vector<Samples>& group_samples, size_t group_idx) const -> Samples{
	return Samples(group_samples[group_idx].begin(), group_samples[group_idx].end());
}

auto GroupedScores::getOtherNames(const std::vector<Samples>& group_samples, size_t group_idx) const -> Samples{
	Samples result;
	for(size_t i=0; i < group_samples.size(); i++){
		if(i == group_idx) continue;
		result.insert(result.begin(), group_samples[i].begin(), group_samples[i].end());
	}
	return result;
}

void GroupedScores::addToResult(DenseMatrix& result, const Scores& gene_set, size_t column_idx) const{
	size_t row_idx = -1;
	for(const auto& score : gene_set.scores()) {
		row_idx++;
		result(row_idx, column_idx) = score;
	}
}


