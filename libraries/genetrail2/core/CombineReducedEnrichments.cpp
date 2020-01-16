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
	void CombineReducedEnrichments::readConfig(const std::string& sampleOutDirs, const std::string& matrixOutFiles){
		parseMatrixOutFiles(matrixOutFiles);
		parseSampleOutDirs(sampleOutDirs);
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
		allReaders.clear();
		allReaders.resize(categoryDBs.size());
		samples.clear();
		std::ifstream reader(sampleOutDirs);
		const std::string split_string = "\t";
		std::vector<std::string> fields;
		if(!reader) throw IOError("Invalid input stream!");
		
		for(std::string line; std::getline(reader, line);){
			boost::split(fields, line, boost::is_any_of(split_string), boost::token_compress_on);
			if(fields.size() < 2) continue;
			samples.emplace_back(fields[0]);
			const auto& outDir = fields[1];
			
			int idx = -1;
			for(const auto& categoryDB: categoryDBs){
				idx++;
				std::string outfile = outDir + "/" + categoryDB + ".txt";
				allReaders[idx].emplace_back(outfile);
				if(!allReaders[idx].back().is_open()){
					throw new IOError(
						"Could not open input file: " + outfile
					);
				}
			}
		}
	}
	
	void CombineReducedEnrichments::write(){
		for(unsigned int idx=0; idx < categoryDBs.size(); idx++){
			writeHeader(idx);
			std::cout << "Writing category " << idx+1 << " of " <<categoryDBs.size() << std::endl;
			writeBody(idx);
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
	
	void CombineReducedEnrichments::writeBody(unsigned int idx){
		const std::string split_string = "\t";
		std::vector<std::string> fields;
		std::string line;
		bool good = !allReaders[idx].empty();
		auto& writer = allWriters[idx];
		
		while(good){
			bool first = true;
			bool comment = false;
			for(auto& reader: allReaders[idx]){
				if(!std::getline(reader, line)){
					good = false;
					break;
				}
				if(boost::starts_with(line, "#")){
					comment = true;
					continue;
				}
				boost::split(fields, line, boost::is_any_of(split_string), boost::token_compress_on);
				boost::replace_all(fields[0], "_", " ");
				writer << (first ? fields[0] + "\t" : "\t") << fields[1];
				first = false;
			}
			if(comment){
				comment = false;
			} else {
				writer << "\n";
			}
		}
	}
}




