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
#ifndef GT2_GROUPED_SCORES_H
#define GT2_GROUPED_SCORES_H

#include "macros.h"
#include "ORAGroupPreference.h"
#include "DenseMatrix.h"
#include "DenseColumnSubset.h"
#include "Metadata.h"
#include "MatrixHTest.h"

#include <utility>
#include <tuple>
#include <vector>

namespace GeneTrail {

    class GT2_EXPORT GroupedScores {
		public:
			using Samples = std::vector<std::string>;
			
			GroupedScores() = default;
			
			/**
			 * This method accepts a (gene x sample) expression matrix and a Metadata
			 * object with a key for each sample providing a std::string representing
			 * the group of the corresponding sample. This method calculates #group
			 * many tests (specified by method) testing a group versus all other groups.
			 * The resulting (gene x group) matrix is stored in 'result'.
			 */
			void calculateGroupedScores(
				DenseMatrix& matrix,
				const std::vector<Metadata>& metadata,
				const std::string method, DenseMatrix& result
			) const;
			
			/**
			 * Use this method if mean-fold-quotient is to be calculated and the group
			 * means are already precomputed. The first row in the matrix is interpreted
			 * as number of elements in the groups.
			 */
			void calculateGroupedScores(
				DenseMatrix& matrix,
				const std::vector<std::string>& sampleGroups,
				const std::vector<std::string>& referenceGroups,
				const std::string& method,
				DenseMatrix& result
			) const;
			
		private:
			void createSamples(
				const DenseMatrix& matrix,
				const std::map<std::string, std::vector<unsigned int>>& group_indices,
				const std::string& current_group,
				Samples& g1, Samples& g2
			) const;
			
			void addGroupSizeToResult(
				DenseMatrix& result,
				const std::map<std::string, std::vector<unsigned int>>& group_indices
			) const;
			
			void addToResult(DenseMatrix& result, const Scores& gene_set, size_t column_idx) const;
			
			void getIndices(
				const DenseColumnSubset& sub,
				const std::vector<std::string>& groups,
				std::vector<size_t>& indices
			) const;
			
			void getGroupSizes(
				const DenseColumnSubset& sub,
				const std::vector<size_t>& indices,
				std::vector<size_t>& sizes
			) const;
			
			size_t getTotalGroupSizes(
				const std::vector<size_t>& sizes
			) const;
			
			double getMeanScore(
				size_t idx_row,
				const DenseColumnSubset& sub,
				const std::vector<size_t>& indices,
				const std::vector<size_t>& sizes,
				size_t n
			) const;
	};
}

#endif // GT2_CORE_OVER_REPRESENTATION_ANALYSIS_H

