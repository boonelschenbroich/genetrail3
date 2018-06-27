/**
* GeneTrail2 - An efficent library for interpreting genetic data
* Copyright (C) 2018 Tim Kehl tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_CONTINGENCY_MATRIX_PARSER_H
#define GT2_CONTINGENCY_MATRIX_PARSER_H

#include <genetrail2/core/macros.h>
#include <genetrail2/core/Exception.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <vector>
#include <functional>
#include <tuple>

namespace GeneTrail
{
template <typename IntType> class GT2_EXPORT ContingencyMatrixParser
{
  public:
	using int_type = IntType;
	using contingency_table = std::tuple<std::string, int_type, int_type, int_type, int_type>;
	using contingency_tables = std::vector<contingency_table>;

	ContingencyMatrixParser(const std::string& file_name)
	:ctables_()
	{
		read_(file_name);
	}

	const contingency_tables& getContingencyTables() const {
		return ctables_;
	}

  private:
	
	void read_(const std::string& file_name)
	{
		std::ifstream input(file_name);
		if(!input) {
			throw GeneTrail::IOError("File (" + file_name + ") is not open for reading");
		}

		std::vector<std::string> sline;
		for(std::string line; getline(input, line);) {
			boost::trim_if(line, boost::is_any_of("\t "));
			boost::split(sline, line, boost::is_any_of(" \t"), boost::token_compress_on);
			if(sline.size() == 5) {
				ctables_.emplace_back(
					std::make_tuple(
						sline[0],
						boost::numeric_cast<int_type>(boost::lexical_cast<int>(sline[1])),
						boost::numeric_cast<int_type>(boost::lexical_cast<int>(sline[2])),
						boost::numeric_cast<int_type>(boost::lexical_cast<int>(sline[3])),
						boost::numeric_cast<int_type>(boost::lexical_cast<int>(sline[4]))
					)
				);
			} else {
				throw GeneTrail::IOError("Wrong file format.");
			}
		}
	}

	contingency_tables ctables_;
};
}

#endif // GT2_CONTINGENCY_MATRIX_PARSER_H
