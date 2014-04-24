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

#include "SparseMatrixWriter.h"
#include "SparseMatrix.h"

#include <iterator>
#include <iostream>

namespace GeneTrail
{
	uint64_t SparseMatrixWriter::writeBinary(std::ostream& output, const SparseMatrix& matrix) const
	{
		return writeBinary_(output, matrix) + writeData_(output, matrix);
	}

	void SparseMatrixWriter::writeText(std::ostream& output, const SparseMatrix& matrix) const
	{
		writeText_(output, matrix);

		// Write the column names
		output << matrix.rowName(0);
		if(matrix.rowName(0) == "") {
			std::cerr << "Warning: empty row name supplied in column 0.\n";
		}

		for(SparseMatrix::index_type i = 1; i < matrix.rows(); ++i) {
			if(matrix.rowName(i) == "") {
				std::cerr << "Warning: empty row name supplied in column " << i << ".\n";
			}

			output << "\t" << matrix.rowName(i);
		}

		output << "\n";

		// Write the matrix as triplets
		for(int k = 0; k < matrix.matrix().outerSize(); ++k) {
			for (SparseMatrix::SMatrix::InnerIterator it(matrix.matrix(), k); it; ++it) {
				output << it.row() << "\t" << it.col() << "\t" << it.value() << "\n";
			}
		}
	}

	uint64_t SparseMatrixWriter::writeData_(std::ostream& output, const SparseMatrix& matrix) const
	{
		size_t total = 0;

		// Write data
		uint64_t n = matrix.matrix().nonZeros() * sizeof(SparseMatrix::value_type);
		total += writeChunkHeader_(output, 0x5, n);
		output.write((const char*)matrix.matrix().valuePtr(), n);
		total += n;

		// Write outer indices
		n = (matrix.matrix().outerSize() + 1) * sizeof(SparseMatrix::SMatrix::Index);
		total += writeChunkHeader_(output, 0x3, n);
		output.write((const char*)matrix.matrix().outerIndexPtr(), n);
		total += n;

		// Write inner indices
		n = matrix.matrix().nonZeros() * sizeof(SparseMatrix::SMatrix::Index);
		writeChunkHeader_(output, 0x4, n);
		output.write((const char*)matrix.matrix().innerIndexPtr(), n);
		total += n;

		return total;
	}

}
