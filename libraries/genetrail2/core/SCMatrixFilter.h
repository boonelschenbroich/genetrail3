/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2013-2014 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *               2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_SC_MATRIX_FILTER_H
#define GT2_SC_MATRIX_FILTER_H

#include "macros.h"
#include "DenseMatrix.h"
#include "DenseColumnSubset.h"
#include "DenseMatrixWriter.h"

#include <istream>
#include <fstream>
#include <vector>
#include <limits>

namespace GeneTrail{
	struct GT2_EXPORT FilterParams{
		double min_total_count = 0;
		double max_total_count = std::numeric_limits<double>::max();
		double min_features = 0;
		double max_features = std::numeric_limits<double>::max();
		double max_mito = 1;
		double nonzero_threshold = 0.001;
		std::string mito_genes = "";
		
		std::string out_matrix = "";
		std::string out_statistics = "";
	};
	
	class GT2_EXPORT SCMatrixFilter{
	public:
		SCMatrixFilter() = default;
		
		void filterMatrix(DenseMatrix& matrix, const FilterParams& params);
		
	private:
		void fillColumnStatistics(const DenseMatrix& matrix, std::vector<double>& total_count, std::vector<double>& mito_count,
		                      std::vector<double>& nonzero_features, const FilterParams& params);
		bool passFilter(double total_count, double nonzero_features, const FilterParams& params);
		void writeFilteredMatrix(const DenseColumnSubset& result, const FilterParams& params);
		void writeFiles(const std::string& sampleOutDirs, const std::string& matrixOutFiles);
		void writeStatisticsFile(const std::vector<double>& total_count, const std::vector<double>& mito_count,
											 const std::vector<double>& nonzero_features, const std::vector<std::string>& keep,
											 const std::vector<size_t>& keep_idx, const FilterParams& params);
		void writeStatisticsLine(std::ofstream& writer, const std::string& row_name,
											 const std::vector<double>& features,
											 const std::vector<size_t>& keep_idx);

	private:
		
	};
}

#endif //GT2_DENSE_MATRIX_READER_H

