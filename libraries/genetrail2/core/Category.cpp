/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#include "Category.h"

#include <algorithm>
#include <exception>
#include <iostream>

namespace GeneTrail
{
	Category::Category(EntityDatabase* database, const std::string& name)
	    : name_(name), database_(database)
	{
	}

	Category::Category(EntityDatabase* database)
	    : database_(database)
	{
	}

	const std::string& Category::name() const { return name_; }

	void Category::setName(const std::string& n) { name_ = n; }
	void Category::setName(std::string&& n) { name_ = std::move(n); }

	const std::string& Category::reference() const { return reference_; }

	void Category::setReference(std::string r) { reference_ = std::move(r); }

	const Metadata& Category::metadata() const { return metadata_; }

	Metadata& Category::metadata() { return metadata_; }

	bool Category::contains(const std::string& id) const
	{
		return contains(database_->index(id));
	}

	bool Category::contains(size_t i) const
	{
		return container_.find(i) != container_.end();
	}

	bool Category::insert(const std::string& id)
	{
		return insert(database_->index(id));
	}

	bool Category::insert(size_t i) { return container_.emplace(i).second; }

	size_t Category::size() const { return container_.size(); }

	bool Category::empty() const { return container_.empty(); }

	bool Category::operator<(const Category& o) const
	{
		return std::lexicographical_compare(
		    container_.begin(), container_.end(), o.container_.begin(),
		    o.container_.end());
	}

	bool Category::operator==(const Category& o) const
	{
		if(size() != o.size()) {
			return false;
		}

		auto it = begin();
		auto jt = o.begin();

		for(; it != end(); ++it, ++jt) {
			if(*it != *jt) {
				return false;
			}
		}

		return true;
	}

	Category Category::intersect(const std::string& name, const Category& a,
	                             const Category& b)
	{
		Category result(a.database_);

		result.setName(name);

		if(a.database_ != b.database_) {
			throw std::invalid_argument("EntityDatabases not compatible.");
		}

		std::set_intersection(
		    a.container_.begin(), a.container_.end(), b.container_.begin(),
		    b.container_.end(),
		    std::inserter(result.container_, std::end(result.container_)));

		return result;
	}

	Category Category::combine(const std::string& name, const Category& a,
	                           const Category& b)
	{
		Category result(a.database_);
		result.setName(name);

		if(a.database_ != b.database_) {
			throw std::invalid_argument("EntityDatabases not compatible.");
		}

		result.container_.reserve(std::max(a.size(), b.size()));

		std::set_union(
		    a.container_.begin(), a.container_.end(), b.container_.begin(),
		    b.container_.end(),
		    std::inserter(result.container_, std::end(result.container_)));

		return result;
	}

	std::ostream& operator<<(std::ostream& strm, const Category& cat)
	{
		strm << &cat << std::endl << std::endl;
		strm << cat.name() << "\t" << cat.reference();

		std::copy(cat.names().begin(), cat.names().end(),
		          std::ostream_iterator<std::string>(strm, "\t"));

		return strm;
	}
}
