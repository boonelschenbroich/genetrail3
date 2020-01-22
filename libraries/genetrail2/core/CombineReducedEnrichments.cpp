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

#include "CombineReducedEnrichments.h"

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
	void CombineReducedEnrichments::writeFiles(const std::string& sampleOutDirs, const std::string& matrixOutFiles){
		parseMatrixOutFiles(matrixOutFiles);
		parseSampleOutDirs(sampleOutDirs);
		for(size_t idx=0; idx < categoryDBs.size(); idx++){
			write(idx);
		}
	}
	
	void CombineReducedEnrichments::parseMatrixOutFiles(const std::string& matrixOutFiles){
		categoryDBs.clear();
		allWriters.clear();
		std::ifstream reader(matrixOutFiles);
		const std::string split_string = "\t";
		std::vector<std::string> fields;
		
		if(!reader) throw IOError("Invalid input stream!");
		for(std::string line; std::getline(reader, line);){
			boost::split(fields, line, boost::is_any_of(split_string), boost::token_compress_on);
			if(fields.size() < 2) continue;
			categoryDBs.emplace_back(fields[0]);
			allWriters.emplace_back(fields[1]);
			if(!allWriters.back().is_open()){
				throw IOError("Could not open output file: " + fields[1]);
			}
		}
	}
	
	void CombineReducedEnrichments::parseSampleOutDirs(const std::string& sampleOutDirs){
		samples.clear();
		dirs.clear();
		const std::string split_string = "\t";
		std::vector<std::string> fields;
		
		std::ifstream reader(sampleOutDirs);
		if(!reader) throw IOError("Invalid input stream!");
		
		for(std::string line; std::getline(reader, line);){
			boost::split(fields, line, boost::is_any_of(split_string), boost::token_compress_on);
			if(fields.size() < 2) continue;
			samples.emplace_back(fields[0]);
			dirs.emplace_back(fields[1]);
		}
	}
	
	void CombineReducedEnrichments::write(size_t idx){
		const auto& category = categoryDBs[idx];
		auto& writer = allWriters[idx];
		std::map<std::string, std::vector<std::string>> result;
		
		for(size_t idx_sample=0; idx_sample < samples.size(); idx_sample++){
			std::cout << "Processing category " << (idx+1) << "/";
			std::cout << categoryDBs.size() << " - sample " << (idx_sample+1);
			std::cout << "/" << samples.size() << std::endl;
			readSample(result, idx_sample, category);
		}
		
		std::cout << "Writing matrix for category " << (idx+1) << "/" << categoryDBs.size() << std::endl;
		writeMatrix(result, writer);
	}
	
	void CombineReducedEnrichments::readSample(std::map<std::string, std::vector<std::string>>& result, size_t idx_sample, const std::string& categoryDB){
		const std::string file = dirs[idx_sample] + "/" + categoryDB + ".txt";
		
		const std::string split_string = "\t";
		std::vector<std::string> fields;
		std::string category;
		std::ifstream reader(file);
		if(!reader) throw IOError("Could not open input file " + file);
		
		for(std::string line; std::getline(reader, line);){
			if(boost::starts_with(line, "#")) continue;
			boost::split(fields, line, boost::is_any_of(split_string), boost::token_compress_on);
			if(fields.size() < 2) continue;
			category = fields[0];
			auto search = result.find(category);
			if(search == result.end()){
				result.emplace(category, std::vector<std::string>(samples.size()));
				search = result.find(category);
			}
			search->second[idx_sample] = fields[1];
		}
	}
	
	void CombineReducedEnrichments::writeHeader(unsigned int idx){
		auto& writer = allWriters[idx];
		bool first = true;
		for(const auto& sample: samples){
			writer << (first ? "" : "\t") << sample;
			first = false;
		}
		writer << "\n";
	}
	
	void CombineReducedEnrichments::writeMatrix(std::map<std::string, std::vector<std::string>>& result, std::ofstream& writer){
		bool first = true;
		for(const auto& sample: samples){
			writer << (first ? "" : "\t") << sample;
			first = false;
		}
		writer << std::endl;
		
		for(auto& categoryResults: result){
			std::string category = categoryResults.first;
			boost::replace_all(category, "_", " ");
			writer << category;
			for(const auto& res: categoryResults.second){
				writer << "\t" << res;
			}
			writer << std::endl;
		}
	}
}




