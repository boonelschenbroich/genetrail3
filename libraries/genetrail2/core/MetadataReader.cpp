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

#include "MetadataReader.h"

#include <vector>
#include <iostream>
#include <string>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>

#include "Exception.h"

namespace GeneTrail{
	void MetadataReader::skipEmptyLines_(std::istream& input, std::string& line) const{
		while(std::getline(input, line)){
			// TODO this is dangerous if there is an empty row name
			boost::trim(line);
			if(!line.empty()) break;
		}
	}

	Metadata MetadataReader::readMetadataFile(std::istream& input, const std::string& column) const
	{
		if(!input) throw IOError("Invalid input stream!");
		return textRead_(input, column);
	}

	Metadata MetadataReader::textRead_(std::istream& input, const std::string& column) const{
		std::string line;
		// Get the first interesting line
		skipEmptyLines_(input, line);

		std::vector<std::string> fields;
		boost::split(fields, line, boost::is_any_of("\t"), boost::token_compress_on);
		
		unsigned int column_index = 0;
		for(const std::string& field: fields){
			if(field == column) break;
			++column_index;
		}
		if(column_index == fields.size()){
			throw IOError("The metadata file doesn't contain the requested column " + column + ".");
		}
		column_index++;

		Metadata metadata;
		unsigned int cur_line = 1;
		unsigned int num_fields = fields.size() + 1;

		while(std::getline(input, line)){
			boost::trim(line);
			if(line.empty()) continue;

			boost::split(fields, line, boost::is_any_of("\t"), boost::token_compress_on);

			if(fields.size() != num_fields){
				throw IOError(
					"Expected " + boost::lexical_cast<std::string>(num_fields) +
					" columns in line " + boost::lexical_cast<std::string>(cur_line) +
					", got " + boost::lexical_cast<std::string>(fields.size())
				);
			}
			metadata.add(fields[0], fields[column_index]);
		}
		return metadata;
	}
}

