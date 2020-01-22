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

#ifndef GT2_DENSE_MATRIX_READER_H
#define GT2_DENSE_MATRIX_READER_H

#include "macros.h"

#include <istream>
#include <vector>
#include <map>
#include <set>
#include <initializer_list>

namespace GeneTrail{
	class GT2_EXPORT CombineReducedEnrichments{
	public:
		CombineReducedEnrichments() = default;
		
		void writeFiles(const std::string& sampleOutDirs, const std::string& matrixOutFiles);

	private:
		using WriterPerCategoryDB = std::vector<std::ofstream>;
		using CategoryDBs = std::vector<std::string>;
		using Samples = std::vector<std::string>;
		
		WriterPerCategoryDB allWriters;
		CategoryDBs categoryDBs;
		std::vector<std::string> samples;
		std::vector<std::string> dirs;
		
		void parseMatrixOutFiles(const std::string& matrixOutFiles);
		void parseSampleOutDirs(const std::string& sampleOutDirs);
		void write(size_t idx);
		void readSample(std::map<std::string, std::vector<std::string>>& result, size_t idx_sample, const std::string& category);
		void writeMatrix(std::map<std::string, std::vector<std::string>>& result, std::ofstream& writer);
		void writeHeader(unsigned int idx);
	};
}

#endif //GT2_DENSE_MATRIX_READER_H

