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

#ifndef GT2_DENSE_MATRIX_WRITER_H
#define GT2_DENSE_MATRIX_WRITER_H

#include "config.h"

#include <ostream>
#include <vector>

namespace GeneTrail
{
	class DenseMatrix;

	class GT2_EXPORT DenseMatrixWriter
	{
		public:
			void writeText(std::ostream& output, const DenseMatrix& matrix) const;

			/**
			 * Reads a matrix from a binary file.
			 * 
			 * \see DenseMatrixReader::binaryRead_
			 */
			void writeBinary(std::ostream& output, const DenseMatrix& matrix) const;

	private:
			void writeChunkHeader_(std::ostream& output, uint8_t type, uint64_t size) const;
			void writeHeader_     (std::ostream& output, const DenseMatrix& matrix) const;
			void writeData_       (std::ostream& output, const DenseMatrix& matrix) const;
			void writeRowNames_   (std::ostream& output, const DenseMatrix& matrix) const;
			void writeColNames_   (std::ostream& output, const DenseMatrix& matrix) const;
			void writeNames_      (std::ostream& output, const std::vector<std::string>& names) const;
	};
}

#endif //GT2_DENSE_MATRIX_WRITER_H
