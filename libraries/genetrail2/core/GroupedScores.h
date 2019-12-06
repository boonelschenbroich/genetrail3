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
			void calculateGroupedScores(DenseMatrix& matrix, const Metadata& metadata,
										const std::string method, DenseMatrix& result) const;
			
		private:
			Samples getNames(const std::vector<Samples>& group_samples, size_t group_idx) const;
			
			Samples getOtherNames(const std::vector<Samples>& group_samples, size_t group_idx) const;
			
			void addToResult(DenseMatrix& result, const Scores& gene_set, size_t column_idx) const;
	};
}

#endif // GT2_CORE_OVER_REPRESENTATION_ANALYSIS_H

