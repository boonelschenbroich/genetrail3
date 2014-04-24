/*
 * GeneTrail2 - An efficent library for interpreting genetic data
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

#include "RMAExpressMatrixReader.h"

#include <cstdint>

#include "DenseMatrix.h"

namespace GeneTrail
{
	DenseMatrix RMAExpressMatrixReader::read(std::istream& input, unsigned int opts) const
	{
		//TODO: interpret the passed options

		if(readString_(input) != "RMAExpressionValues") {
			return DenseMatrix(0,0);
		}

		uint32_t version;
		input.read((char*)&version, 4);

		if(version != 1) {
			return DenseMatrix(0,0);
		}

		skipString_(input); // Skip the RMAExpress version number
		skipString_(input); // Skip the chiptype

		uint32_t nrows, ncols;
		input.read((char*)&ncols, 4);
		input.read((char*)&nrows, 4);

		std::vector<std::string> colnames(ncols);
		std::vector<std::string> rownames(nrows);

		readStringVector_(input, colnames);
		readStringVector_(input, rownames);

		DenseMatrix result(std::move(rownames), std::move(colnames));

		input.read((char*)result.matrix().data(), nrows * ncols * sizeof(double));

		return result;
	}

	void RMAExpressMatrixReader::readStringVector_(std::istream& input, std::vector<std::string>& vec) const
	{
		for(auto& s : vec) {
			s = readString_(input);
		}
	}

	void RMAExpressMatrixReader::skipString_(std::istream& input) const
	{
		uint32_t size;
		input.read((char*)&size, 4);

		input.seekg(size, std::ios::cur);
	}

	std::string RMAExpressMatrixReader::readString_(std::istream& input) const
	{
		uint32_t size;
		input.read((char*)&size, 4);

		std::string result(size, '\0');
		input.read(&result[0], size);

		return result;
	}


}

