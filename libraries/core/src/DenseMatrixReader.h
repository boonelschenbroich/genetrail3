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

#ifndef GT2_DENSE_MATRIX_READER_H
#define GT2_DENSE_MATRIX_READER_H

#include "config.h"

#include <istream>
#include <vector>

namespace GeneTrail
{
	class DenseMatrix;

	class GT2_EXPORT DenseMatrixReader
	{
		public:
			enum ReaderOptions
			{
			NO_OPTIONS          = 0,
			READ_ROW_NAMES      = 1 << 0,
			READ_COL_NAMES      = 1 << 1,
			TRANSPOSE           = 1 << 2,
			ADDITIONAL_COL_NAME = 1 << 3
			};

			static unsigned int defaultOptions();

			virtual DenseMatrix read(std::istream& input, unsigned int opts = defaultOptions()) const;

		private:
			DenseMatrix textRead_ (std::istream& input, unsigned int opts) const;
			
			enum ChunkType {
				HEADER   = 0x00,
				ROWNAMES = 0x01,
				COLNAMES = 0x02,
				DATA     = 0x03
			};

			/**
			 * If isBinary_ returns true, binaryRead_ will attempt to read a binary
			 * matrix from the specified file.
			 * 
			 * The file format is based on chunks, that, with the exception of the
			 * header can be stored in an arbitrary order.
			 * 
			 * The basic chunk format is as follows:
			 * 
			 * CHUNK-ID  : uint8_t   --- The identifier of the chunk
			 * CHUNK-SIZE: uint64_t  --- The size of the CHUNK-DATA field in bytes
			 *                           forwarding the file cursor CHUNK-SIZE bytes
			 *                           will either reach the next CHUNK or EOF
			 * CHUNK-DATA: This field depends on the current chunk
			 * 
			 * Chunk Types (ID):
			 *  * HEADER (0x00):
			 *   - ROW-COUNT:     uint32_t --- The number of rows of the matrix
			 *   - COL-COUNT:     uint32_t --- The number of columns of the matrix
			 *   - STORAGE-ORDER: uint8_t  --- Zero if the matrix is stored in row major
			 *                                 non-zero if it is stored in column major
			 *  * ROWNAMES (0x01):
			 *   Contains ROW-COUNT zero terminated names.
			 *
			 *  * COLNAMES (0x02):
			 *   Contains COL-COUNT zero terminated names.
			 * 
			 *  * DATA (0x03):
			 *   Contains ROW-COUNT * COL-COUNT double values of entries
			 *   in the storage order specified in the header.
			 */
			DenseMatrix binaryRead_(std::istream& input, unsigned int opts = NO_OPTIONS) const;

			void skipEmptyLines_(std::istream& input, std::string& line) const;

			/**
			 * This method checks the magic number of a stream in order to decide
			 * whether the matrix is stored in binary or text format.
			 * 
			 * The magic number searched for is the signed 12 byte sequence
			 * "BINARYMATRIX"
			 * 
			 */
			bool isBinary_(std::istream& input) const;

			void readRowNames_(std::istream& input, DenseMatrix& result, uint64_t chunk_size) const;
			void readColNames_(std::istream& input, DenseMatrix& result, uint64_t chunk_size) const;
			void readNames_   (std::istream& input, DenseMatrix& result, uint64_t chunk_size, std::vector<std::string>& names) const;
			void readData_    (std::istream& input, DenseMatrix& result, uint8_t storage_order) const;
			void readChunkHeader_(std::istream& input, uint8_t& chunk_type, uint64_t& chunk_size) const;
			DenseMatrix readHeader_(std::istream& input, uint8_t& storage_order) const;
			
	};
}

#endif //GT2_DENSE_MATRIX_READER_H
