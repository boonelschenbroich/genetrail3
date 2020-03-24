/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2020 Nico Gerstner <dstoeckel@bioinf.uni-sb.de>
 *               2020 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#include "SCMatrixFilter.h"

#include <vector>
#include <fstream>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

#include "Exception.h"
#include <iostream>

namespace GeneTrail
{
	 void SCMatrixFilter::filterMatrix(DenseMatrix& matrix, std::set<std::string> mito_genes, const FilterParams& params){
		std::vector<double> total_count(matrix.cols(), 0.0);
		std::vector<double> mito_count(matrix.cols(), 0.0);
		std::vector<double> nonzero_features(matrix.cols(), 0.0);

		fillColumnStatistics(matrix, mito_genes, total_count, mito_count, nonzero_features, params);
		
		std::cout << "Filtering cells..." << std::endl;
		std::vector<std::string> keep;
		std::vector<size_t> keep_idx;
		for(size_t idx_cell=0; idx_cell < matrix.cols(); idx_cell++){
			if(passFilter(total_count[idx_cell], nonzero_features[idx_cell], mito_count[idx_cell], params)){
				keep.push_back(matrix.colNames()[idx_cell]);
				keep_idx.push_back(idx_cell);
			}
		}
		std::cout << "Writing matrix..." << std::endl;
		DenseColumnSubset result = DenseColumnSubset(&matrix, keep_idx.begin(), keep_idx.end());
		
		std::cout << "Writing statistics file..." << std::endl;
		writeFilteredMatrix(result, params);
		writeStatisticsFile(total_count, mito_count, nonzero_features, keep, keep_idx, params);
	}
	
	void SCMatrixFilter::fillColumnStatistics(const DenseMatrix& matrix, std::set<std::string> mito_genes, std::vector<double>& total_count, std::vector<double>& mito_count,
		                      std::vector<double>& nonzero_features, const FilterParams& params
	){
		std::vector<std::string> rownames = matrix.rowNames();
		for(size_t idx_row=0; idx_row < matrix.rows(); idx_row++){
			if(idx_row % 1000 == 0){
				std::cout << "Inspecting row " << (idx_row+1) << "/" << matrix.rows() << std::endl;
			}
			auto it = mito_genes.find(rownames[idx_row]);
			for(size_t idx_col=0; idx_col < matrix.cols(); idx_col++){
				double v = matrix(idx_row, idx_col);
				total_count[idx_col] += v;
				if(v > params.nonzero_threshold){
					nonzero_features[idx_col]++;
				}
				if(it != mito_genes.end()) {
					mito_count[idx_col] += v;
				}
			}
		}
	}
	
	bool SCMatrixFilter::passFilter(double total_count, double nonzero_features, double mito_count, const FilterParams& params){
		if(total_count > params.max_total_count) return false;
		if(total_count < params.min_total_count) return false;
		if(nonzero_features > params.max_features) return false;
		if(nonzero_features < params.min_features) return false;
		if((mito_count / total_count) > params.max_mito) return false;
		return true;
	}
	
	void SCMatrixFilter::writeFilteredMatrix(const DenseColumnSubset& result, const FilterParams& params){
		DenseMatrixWriter writer;
		std::ofstream out(params.out_matrix);
		writer.writeText(out, result);
	}
	
	void SCMatrixFilter::writeStatisticsFile(const std::vector<double>& total_count, const std::vector<double>& mito_count,
											 const std::vector<double>& nonzero_features, const std::vector<std::string>& keep,
											 const std::vector<size_t>& keep_idx, const FilterParams& params
	){
		std::ofstream writer(params.out_statistics);
		bool first = true;
		for(const auto& header: keep){
			writer << (first ? "" : "\t") << header;
			first = false;
		}
		writer << std::endl;

		std::vector<double> mito_percentage;
		mito_percentage.reserve(mito_count.size());
		for(size_t i=0; i<total_count.size(); ++i){
			mito_percentage.emplace_back(mito_count[i] / total_count[i]);
		}
		
		writeStatisticsLine(writer, "total_count", total_count, keep_idx);
		//writeStatisticsLine(writer, "mito_count", mito_count, keep_idx);
		writeStatisticsLine(writer, "mito_percentage", mito_percentage, keep_idx);
		writeStatisticsLine(writer, "nonzero_features", nonzero_features, keep_idx);
		writer.close();
	}
	
	void SCMatrixFilter::writeStatisticsLine(std::ofstream& writer, const std::string& row_name,
											 const std::vector<double>& features,
											 const std::vector<size_t>& keep_idx
	){
		writer << row_name;
		for(size_t keep: keep_idx){
			writer << '\t' << features[keep];
		}
		writer << std::endl;
	}
}




