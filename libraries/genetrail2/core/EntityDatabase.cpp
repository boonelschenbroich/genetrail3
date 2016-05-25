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
#include "EntityDatabase.h"

#include "Exception.h"

namespace GeneTrail
{
	GT2_EXPORT std::shared_ptr<EntityDatabase> EntityDatabase::global =
	    std::make_shared<EntityDatabase>();

	void EntityDatabase::clear()
	{
		name_to_index_.clear();
		db_.clear();
	}

	size_t EntityDatabase::index(const std::string& name)
	{
		auto res = name_to_index_.find(name);

		if(res == name_to_index_.end()) {
			res = name_to_index_.emplace(name, db_.size()).first;
			db_.push_back(name);
		}

		return res->second;
	}

	size_t EntityDatabase::index(const std::string& name) const
	{
		auto res = name_to_index_.find(name);

		if(res == name_to_index_.end()) {
			throw UnknownEntry(name);
		}

		return res->second;
	}
}
