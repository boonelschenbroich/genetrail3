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

#ifndef GT2_METADATA_READER_H
#define GT2_METADATA_READER_H

#include "macros.h"
#include "Metadata.h"

#include <istream>
#include <string>

namespace GeneTrail
{
	class GT2_EXPORT MetadataReader
	{
		public:
			/**
			 * Reads a Metadata from the provided input stream.
			 *
			 * @param input a (seekable) stream
			 * @param column the column in the metadata file that should be parsed along
			 *               with the row names
			 * 
			 * @throws IOError an IOError is thrown if the provided stream is invalid or the
			 *                 metadata description is invalid.
			 */
			virtual Metadata readMetadataFile(std::istream& input, const std::string& column) const;

		private:

			Metadata textRead_ (std::istream& input, const std::string& column) const;
			void skipEmptyLines_(std::istream& input, std::string& line) const;
			Metadata readHeader_(std::istream& input, uint8_t& storage_order) const;
	};
}

#endif //GT2_DENSE_MATRIX_READER_H

