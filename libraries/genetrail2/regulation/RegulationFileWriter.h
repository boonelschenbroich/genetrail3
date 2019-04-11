/**
* GeneTrail2 - An efficent library for interpreting genetic data
* Copyright (C) 2016 Tim Kehl tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_REGULATION_FILE_WRITER_H
#define GT2_REGULATION_FILE_WRITER_H

#include <genetrail2/core/macros.h>
#include <genetrail2/core/Exception.h>
#include <genetrail2/core/Matrix.h>
#include <genetrail2/core/DenseMatrix.h>

#include "RegulationFile.h"

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <unordered_set>
#include <map>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>

namespace GeneTrail
{
template <typename NameDatabase, typename ValueType> class GT2_EXPORT RegulationFileWriter
{
  public:
	using value_type = ValueType;
	using Regulation = std::tuple<size_t, size_t, value_type>;

	RegulationFileWriter(NameDatabase& name_database,
	                     const std::unordered_set<size_t> test_set,
	                     const std::string& file, value_type default_value)
	    : regulation_file_(name_database.size(), MAX_MATRIX_INDEX)

	{
		read_(name_database, test_set, file, default_value);
	}

	RegulationFile<value_type>& getRegulationFile() { return regulation_file_; }

  private:
	static constexpr Matrix::index_type MAX_MATRIX_INDEX =
	    std::numeric_limits<Matrix::index_type>::max();

	void read_(NameDatabase& name_database,
	           const std::unordered_set<size_t>& test_set,
	           const std::string& file, value_type default_value)
	{
		std::ifstream input(file);
		if(!input) {
			throw GeneTrail::IOError("File (" + file +
			                         ") is not open for reading");
		}

		std::vector<std::string> sline;
		for(std::string line; getline(input, line);) {
			boost::trim_if(line, boost::is_any_of("\t "));
			boost::split(sline, line, boost::is_any_of(" \t"), boost::token_compress_on);
			if(sline.size() == 2) {
				addRegulation_(name_database, test_set, sline[0], sline[1],
				               default_value);
			} else if(sline.size() == 3) {
				addRegulation_(name_database, test_set, sline[0], sline[1],
				               boost::lexical_cast<value_type>(sline[2]));
			} else {
				throw GeneTrail::IOError("Wrong file format.");
			}
		}
	}

	void addRegulation_(NameDatabase& name_database,
	                    const std::unordered_set<size_t>& test_set,
	                    const std::string& regulator, const std::string& target,
	                    value_type value)
	{
		auto regulator_idx = name_database(regulator);
		auto target_idx = name_database(target);

		if(regulator_idx == MAX_MATRIX_INDEX ||
		   target_idx == MAX_MATRIX_INDEX) {
			return;
		}

		regulation_file_.increaseNumberOfTargets(regulator_idx);

		if(test_set.find(target_idx) == test_set.end()) {
			return;
		}

		regulation_file_.addRegulation(regulator_idx, target_idx, value);
	}

	RegulationFile<value_type> regulation_file_;
};
}

#endif // GT2_REGULATION_FILE_PARSER_H
 
