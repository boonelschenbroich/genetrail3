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
	 void SCMatrixFilter::filterMatrix(const std::string& matrix, std::set<std::string> mito_genes, const FilterParams& params){
		std::vector<double> total_count;
		std::vector<double> mito_count;
		std::vector<double> nonzero_features;

		fillColumnStatistics(matrix, mito_genes, total_count, mito_count, nonzero_features, params);
		
		std::cout << "Filtering cells..." << std::endl;
		std::vector<std::string> keep;
		std::vector<size_t> keep_idx;
		for(size_t idx_cell=0; idx_cell < cols; idx_cell++){
			if(passFilter(total_count[idx_cell], nonzero_features[idx_cell], mito_count[idx_cell], params)){
				keep.push_back(col_names[idx_cell]);
				keep_idx.push_back(idx_cell);
			}	
		}
		
		std::cout << "Writing matrix..." << std::endl;
		writeFilteredMatrix(matrix, keep_idx, params);
		
		std::cout << "Writing statistics file..." << std::endl;
		writeStatisticsFile(total_count, mito_count, nonzero_features, keep, keep_idx, params);
	}
	
	void SCMatrixFilter::fillColumnStatistics(const std::string& matrix, std::set<std::string> mito_genes, std::vector<double>& total_count, std::vector<double>& mito_count,
		                      std::vector<double>& nonzero_features, const FilterParams& params
	){
		std::ifstream reader(matrix);
		std::vector<std::string> fields;
		int row_count = 0;
		
		// parse header
		std::string line;
		std::getline(reader, line);
		boost::split(col_names, line, boost::is_any_of("\t,"), boost::token_compress_on);
		cols = col_names.size();
		
		total_count = std::vector<double>(cols, 0.0);
		mito_count = std::vector<double>(cols, 0.0);
		nonzero_features = std::vector<double>(cols, 0.0);
		
		while(std::getline(reader, line)){
			boost::split(fields, line, boost::is_any_of("\t,"), boost::token_compress_on);
			row_count++;
			if(fields.empty()) continue;
			
			if(row_count % 500 == 0){
				std::cout << "Inspecting row " << (row_count+1) << std::endl;
			}
			row_names.push_back(fields[0]);
			auto it = mito_genes.find(fields[0]);
			
			size_t idx_col;
			std::string name;
			double v;
			std::istringstream iss(line);
			iss >> name;
			for(size_t i=1; i < fields.size(); i++){
				iss >> v;
				idx_col = i-1;
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
	
	void SCMatrixFilter::writeFilteredMatrix(const std::string& matrix, const std::vector<size_t>& keep_idx, const FilterParams& params){
		std::ifstream reader(matrix);
		std::ofstream writer(params.out_matrix);
		std::vector<std::string> fields;
		std::string line;
		int row_count = 0;
		int col_offset = 0;
		
		if(keep_idx.empty()) return;
		
		while(std::getline(reader, line)){
			boost::split(fields, line, boost::is_any_of("\t,"), boost::token_compress_on);
			row_count++;
			
			if(row_count % 1000 == 0){
				std::cout << "Writing filtered row " << (row_count+1) << std::endl;
			}
			
			// write row name if exists
			if(row_count > 1){
				writer << fields[0] << "\t";
			}
			// write remaining of the columns
			writer << fields[keep_idx[0] + col_offset];
			for(size_t i=1; i < keep_idx.size(); i++){
				writer << "\t" << fields[keep_idx[i] + col_offset];
			}
			writer << std::endl;
			
			if(col_offset == 0) col_offset = 1;
		}
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




