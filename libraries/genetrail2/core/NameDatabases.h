/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2016 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_CORE_NAME_DATABASES_H
#define GT2_CORE_NAME_DATABASES_H

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "macros.h"

#include "DenseMatrix.h"

#include <string>
#include <iostream>
#include <fstream>

namespace GeneTrail
{
	
struct GT2_EXPORT MatrixNameDatabase
{
	MatrixNameDatabase(DenseMatrix* matrix) : matrix_(matrix) {}

	std::string operator()(size_t index) { return matrix_->rowNames()[index]; }

	size_t operator()(std::string name) { return matrix_->rowIndex(name); }

	size_t size() { return matrix_->rows(); }

	DenseMatrix* matrix_;
};
	
struct GT2_EXPORT MapNameDatabase
{
	MapNameDatabase() {}

	MapNameDatabase(std::string file) { initialize(file); }

	std::string operator()(size_t index) { return index2name[index]; }

	size_t operator()(std::string name)
	{
		std::map<std::string, size_t>::iterator iter = name2index.find(name);
		if(iter == name2index.end()) {
			index2name.emplace_back(name);
			name2index[name] = index2name.size() - 1;
		}
		return name2index[name];
	}

	size_t size() { return index2name.size(); }

	void initialize(std::string file)
	{
		std::cout << "INFO: Initializing name database" << std::endl;
		std::ifstream input(file);
		if(!input) {
			throw GeneTrail::IOError("File (" + file +
			                         ") is not open for reading");
		}

		std::vector<std::string> sline(2);
		for(std::string line; getline(input, line);) {
			boost::split(sline, line, boost::is_any_of(" \t"), boost::token_compress_on);
			(*this)(sline[0]);
			(*this)(sline[1]);
		}
	}

	std::map<std::string, size_t> name2index;
	std::vector<std::string> index2name;
};
}

#endif // GT2_CORE_NAME_DATABASES_H
