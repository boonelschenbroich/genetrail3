/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2013 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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

#include "MatrixTools.h"
#include "Exception.h"

#include <ostream>
#include <vector>

using namespace GeneTrail;

std::vector<unsigned int> MatrixTools::getIndices(const DenseMatrix& matrix, const std::vector<std::string>& colnames, const std::string& groupname)
{
	std::vector<unsigned int> indices;
	indices.reserve(colnames.size());
	for(const auto& s : colnames) {
		if(matrix.hasCol(s)) {
			indices.emplace_back(matrix.colIndex(s));
		} else {
			std::cerr << "WARNING: Could not find column \"" + s + "\".\n";
		}
	}

	if(indices.empty() && !colnames.empty()) {
		throw EmptyGroup(groupname);
	}

	return indices;
}

std::tuple<DenseColumnSubset, DenseColumnSubset> MatrixTools::splitMatrix(DenseMatrix& matrix, const std::vector<std::string>& reference, const std::vector<std::string>& test)
{
	return std::make_tuple(
		DenseColumnSubset(
			&matrix, getIndices(matrix, reference, "reference")),
		DenseColumnSubset(
			&matrix, getIndices(matrix, test, "test")));
}
