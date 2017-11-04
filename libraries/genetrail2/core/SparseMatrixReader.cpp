/*
 * GeneTrail2 - An efficient library for interpreting genetic data
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

#include "SparseMatrixReader.h"

#include <vector>
#include <deque>
#include <iostream>
#include <cstring>
#include <iostream>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>

#include "SparseMatrix.h"
#include "Exception.h"

namespace GeneTrail
{
	unsigned int SparseMatrixReader::defaultOptions()
	{
		return READ_COL_NAMES | READ_ROW_NAMES;
	}

	void SparseMatrixReader::skipEmptyLines_(std::istream& input, std::string& line) const
	{
		while(std::getline(input, line))
		{
			// TODO this is dangerous if there is an empty row name
			boost::trim(line);

			if(!line.empty())
			{
				break;
			}
		}
	}

	bool SparseMatrixReader::isBinary_(std::istream& input) const
	{
		char magic[12];
		input.read(magic, 12);
		return strncmp(magic, "BINARYMATRIX", 12) == 0;
	}

	SparseMatrix SparseMatrixReader::read(std::istream& input, unsigned int opts) const
	{
		if(!input) {
			throw IOError("Invalid input stream!");
		}

		if(isBinary_(input)) {
			return binaryRead_(input);
		}

		input.seekg(0, std::ios::beg);
		return textRead_(input, opts);
	}

	void SparseMatrixReader::readChunkHeader_(std::istream& input, uint8_t& chunk_type, uint64_t& size) const
	{
		input.read(reinterpret_cast<char*>(&chunk_type), 1);
		input.read(reinterpret_cast<char*>(&size), 8);
	}
	
	SparseMatrix SparseMatrixReader::readHeader_(std::istream& input, uint8_t& storage_order) const
	{
		uint32_t row_count;
		uint32_t col_count;

		input.read(reinterpret_cast<char*>(&row_count), 4);
		input.read(reinterpret_cast<char*>(&col_count), 4);
		input.read(reinterpret_cast<char*>(&storage_order), 1);

		return SparseMatrix(row_count, col_count);
	}

	void SparseMatrixReader::readRowNames_(std::istream& input, SparseMatrix& result, uint64_t chunk_size) const
	{
		// Reserve some storage
		std::vector<std::string> names(result.rows());

		readNames_(input, result, chunk_size, names);

		result.setRowNames(names);
	}

	void SparseMatrixReader::readColNames_(std::istream& input, SparseMatrix& result, uint64_t chunk_size) const
	{
		// Reserve some storage
		std::vector<std::string> names(result.cols());

		readNames_(input, result, chunk_size, names);

		result.setColNames(names);
	}

	void SparseMatrixReader::readNames_(std::istream& input, SparseMatrix&, uint64_t chunk_size, std::vector<std::string>& names) const
	{
		unsigned int i = 0;

		uint64_t bytes_read = 0;

		const int BUFFER_SIZE = 1024;

		char buffer[BUFFER_SIZE];
		while(input.good() && i < names.size() && bytes_read < chunk_size) {
			input.getline(buffer, BUFFER_SIZE, '\0');

			bytes_read += input.gcount(); // Keep track of the number of bytes read

			names[i] += buffer;

			if((input.rdstate() & std::istream::failbit) != 0) {
				input.clear();
			} else {
				++i;
			}
		}

		checkByteMismatch(chunk_size, bytes_read);
	}

	void SparseMatrixReader::readInnerData_(std::istream& input, SparseMatrix& result, uint64_t chunk_size) const
	{
		uint64_t bytes_read = 0;
		result.matrix().resizeNonZeros(chunk_size / sizeof(SparseMatrix::SMatrix::Index));

		// As the internal storage format of matrix is column major this is quite efficient...
		input.read(reinterpret_cast<char*>(result.matrix().innerIndexPtr()), chunk_size);
		bytes_read += input.gcount();

		checkByteMismatch(chunk_size, bytes_read);
	}

	void SparseMatrixReader::readOuterData_(std::istream& input, SparseMatrix& result, uint64_t chunk_size) const
	{
		uint64_t bytes_read = 0;

		// As the internal storage format of matrix is column major this is quite efficient...
		input.read(reinterpret_cast<char*>(result.matrix().outerIndexPtr()), chunk_size);
		bytes_read += input.gcount();

		checkByteMismatch(chunk_size, bytes_read);
	}

	void SparseMatrixReader::readData_(std::istream& input, SparseMatrix& result, uint64_t chunk_size) const
	{
		uint64_t bytes_read = 0;
		result.matrix().resizeNonZeros(chunk_size / sizeof(SparseMatrix::value_type));

		// As the internal storage format of matrix is column major this is quite efficient...
		input.read(reinterpret_cast<char*>(result.matrix().valuePtr()), chunk_size);
		bytes_read += input.gcount();

		checkByteMismatch(chunk_size, bytes_read);
	}

	void SparseMatrixReader::checkByteMismatch(uint64_t expected, uint64_t actual) const
	{
		if(actual != expected) {
			throw IOError("Could not read the specified amount of bytes: " + boost::lexical_cast<std::string>(expected)
			            + ". Got:" + boost::lexical_cast<std::string>(actual));
		}
	}

	SparseMatrix SparseMatrixReader::binaryRead_(std::istream& input, unsigned int) const
	{
		uint8_t chunk_type = 0;
		uint64_t chunk_size = 0;
		
		// Read the header
		readChunkHeader_(input, chunk_type, chunk_size);

		if(chunk_type != 0x0)
		{
			throw IOError("Unexpected chunk: expected 0 (matrix header), got " + boost::lexical_cast<std::string>(chunk_type));
		}

		if(chunk_size != 9) {
			throw IOError("Inconsistent header size: expected 9 got " + boost::lexical_cast<std::string>(chunk_size));
		}

		uint8_t  storage_order;
		SparseMatrix result = readHeader_(input, storage_order);
		result.matrix().makeCompressed();
		
		while(input.good()) {
			readChunkHeader_(input, chunk_type, chunk_size);

			if(input.eof()) {
				break;
			}

			switch(chunk_type) {
				case SparseMatrixReader::HEADER:
					throw IOError("Unexpected chunk: did not expect header chunk!");
				case SparseMatrixReader::ROWNAMES:
					readRowNames_(input, result, chunk_size);
					break;
				case SparseMatrixReader::COLNAMES:
					readColNames_(input, result, chunk_size);
					break;
				case SparseMatrixReader::OUTER_DATA:
					readOuterData_(input, result, chunk_size);
					break;
				case SparseMatrixReader::INNER_DATA:
					readInnerData_(input, result, chunk_size);
					break;
				case SparseMatrixReader::DATA:
					readData_(input, result, chunk_size);
					break;
				default:
					std::cout << "Unknown chunk " << chunk_type << " Skipping!" << std::endl;
					input.seekg(chunk_size, std::ios::cur);
			}
		}

		return result;
	}

	SparseMatrix SparseMatrixReader::textRead_(std::istream& input, unsigned int) const
	{
		std::string line;
		// Get the first interesting line
		skipEmptyLines_(input, line);


		std::vector<std::string> row_names;
		std::vector<std::string> col_names;

		boost::split(row_names, line, boost::is_any_of(" \t"), boost::token_compress_on);
		skipEmptyLines_(input, line);
		boost::split(col_names, line, boost::is_any_of(" \t"), boost::token_compress_on);

		typedef Eigen::Triplet<SparseMatrix::value_type> T;
		std::deque<T> triplets;

		while(std::getline(input, line))
		{
			boost::trim(line);

			if(line.empty())
			{
				continue;
			}

			std::vector<std::string> fields;
			boost::split(fields, line, boost::is_any_of(" \t"), boost::token_compress_on);

			if(fields.size() != 3) {
				throw IOError("Expected 3 identifiers, got " + boost::lexical_cast<std::string>(fields.size()));
			}

			try {
				triplets.push_back(T(
					boost::lexical_cast<int>(fields[0]),
					boost::lexical_cast<int>(fields[1]),
					boost::lexical_cast<SparseMatrix::value_type>(fields[2])
				));
			} catch(boost::bad_lexical_cast& e) {
				throw IOError(std::string("Error while parsing: ") + e.what());
			}
		}

		SparseMatrix res(std::move(row_names), std::move(col_names));
		res.matrix().setFromTriplets(triplets.begin(), triplets.end());

		return res;
	}
}
