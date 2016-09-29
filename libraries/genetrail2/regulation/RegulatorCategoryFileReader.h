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

#ifndef GT2_REGULATOR_CATEGORY_FILE_READER_H
#define GT2_REGULATOR_CATEGORY_FILE_READER_H

#include <genetrail2/core/Category.h>
#include <genetrail2/core/EntityDatabase.h>
#include <genetrail2/core/Exception.h>
#include <genetrail2/core/macros.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <map>
#include <algorithm>
#include <fstream>
#include <memory>

namespace GeneTrail
{
class GT2_EXPORT RegulatorCategoryFileReader
{
  public:
	RegulatorCategoryFileReader(const std::string& file, EntityDatabase* db)
	    : file_(file), db_(db), categories_(), reference_(db_, "reference")
	{
	}

	RegulatorCategoryFileReader(RegulatorCategoryFileReader&&) = default;
	RegulatorCategoryFileReader&
	operator=(RegulatorCategoryFileReader&&) = default;
	RegulatorCategoryFileReader(const RegulatorCategoryFileReader&) = default;

	void parse() { read_(); }

	Category getReference() { return reference_; }

	std::map<std::string, Category>& getCategories() { return categories_; }

  private:
	void read_();

	std::string file_;
	EntityDatabase* db_;
	std::map<std::string, Category> categories_;
	Category reference_;
};
}

#endif // GT2_REGULATOR_CATEGORY_FILE_READER_H
