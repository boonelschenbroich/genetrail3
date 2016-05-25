/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_GEO_GPL_PARSER_H
#define GT2_GEO_GPL_PARSER_H

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <limits>

#include "GEO.h"

#include "macros.h"

namespace GeneTrail
{
	class GT2_EXPORT GPL_Parser : public GEO
	{
		public:

		GPL_Parser();

		~GPL_Parser();

		void downloadGPLFile(const std::string& filename,
		                     const std::string& geo_dir);

		std::map<std::string, std::string>
		annotate(const std::string& geo_dir,
				 const std::string& platform,
		         const std::string& mappingColumn = "Gene ID");

		std::map<std::string, std::string>
		readGPLFile(const std::string& filename,
		            const std::string& mappingColumn = "Gene ID");
	};
}

#endif //GT2_GEO_GPL_PARSER_H
