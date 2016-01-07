/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2014-2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#include "GMTFile.h"

#include "Exception.h"
#include "EntityDatabase.h"

#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/find_iterator.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <functional>

using namespace boost;

// A UnaryFunction object as required by the boost::make_transform_iterator
// function.
using Range = iterator_range<std::string::const_iterator>;
struct copy_range_f : public std::unary_function<Range, std::string>
{
	std::string operator()(const Range& r) const
	{
		return trim_copy(copy_range<std::string>(r));
	}
};

namespace GeneTrail
{
	GMTFile::GMTFile(const std::shared_ptr<EntityDatabase>& db,
	                 const std::string& path, FileOpenMode mode)
	    : CategoryDatabaseFile(path, mode), entity_database_(db)
	{
		if(mode == FileOpenMode::READ) {
			advanceLine_();
		}
	}

	CategoryDatabase GMTFile::read()
	{
		if(!isValid_() || !isReading()) {
			throw IOError("File is not open for reading");
		}

		CategoryDatabase result(entity_database_);

		while(isValid_()) {
			readCategory_(result);
			advanceLine_();
		}

		return result;
	}

	void GMTFile::readCategory_(CategoryDatabase& db)
	{
		auto split_it =
		    make_split_iterator(next_line_, first_finder("\t", is_equal()));
		decltype(split_it) end_it;

		if(split_it == end_it) {
			//TODO: Better exception
			throw IOError("To few arguments in line");
		}

		auto name = copy_range<std::string>(*split_it++);

		if(split_it == end_it) {
			//TODO: Better exception
			throw IOError("To few arguments in line");
		}

		auto url = copy_range<std::string>(*split_it++);

		auto cit = make_transform_iterator(split_it, copy_range_f());
		auto cend = make_transform_iterator(end_it, copy_range_f());

		// Create the category using our newly created transform iterators
		// we move name and reference, as we do not need them any longer
		auto& c = db.addCategory(cit, cend);
		c.setName(std::move(name));
		c.setReference(std::move(url));
	}

	bool GMTFile::write(const CategoryDatabase& db)
	{
		if(!isValid_() || !isWriting()) {
			throw IOError("File is not open for writing");
		}

		for(const auto& cat : db) {
			(*out_strm_) << cat.name() << '\t' << cat.reference();

			for(const auto& s : cat.names()) {
				(*out_strm_) << '\t' << s;
			}

			(*out_strm_) << '\n';
		}

		return true;
	}

	void GMTFile::advanceLine_()
	{
		do {
			std::getline(*in_strm_, next_line_);
			trim(next_line_);
		} while(isValid_() && next_line_ == "");
	}
}

