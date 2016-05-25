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
#include "CategoryDatabase.h"

namespace GeneTrail
{

CategoryDatabase::CategoryDatabase(const std::shared_ptr<EntityDatabase>& db)
    : entity_database_(db)
{
}

void CategoryDatabase::reserve(size_t capacity)
{
	categories_.reserve(capacity);
}

CategoryDatabase::iterator CategoryDatabase::begin()
{
	return categories_.begin();
}

CategoryDatabase::iterator CategoryDatabase::end() { return categories_.end(); }

CategoryDatabase::const_iterator CategoryDatabase::begin() const
{
	return categories_.begin();
}

CategoryDatabase::const_iterator CategoryDatabase::end() const
{
	return categories_.end();
}

CategoryDatabase::const_iterator CategoryDatabase::cbegin() const
{
	return categories_.cbegin();
}

CategoryDatabase::const_iterator CategoryDatabase::cend() const
{
	return categories_.cend();
}

const std::string& CategoryDatabase::name() const { return name_; }

void CategoryDatabase::setName(const std::string& name) { name_ = name; }

Editor& CategoryDatabase::editor() { return editor_; }
const Editor& CategoryDatabase::editor() const { return editor_; }

void CategoryDatabase::setEditor(const Editor& editor) { editor_ = editor; }

const std::string& CategoryDatabase::creationDate() const
{
	return creation_date_;
}

void CategoryDatabase::setCreationDate(const std::string& date)
{
	creation_date_ = date;
}

const std::string& CategoryDatabase::sourceUrl() const { return source_url_; }

void CategoryDatabase::setSourceUrl(const std::string& url)
{
	source_url_ = url;
}

const std::string& CategoryDatabase::identifier() const { return identifier_; }

void CategoryDatabase::setIdentifier(const std::string& identifer)
{
	identifier_ = identifer;
}

const Category& CategoryDatabase::operator[](size_t i) const
{
	return categories_[i];
}

Category& CategoryDatabase::operator[](size_t i) { return categories_[i]; }

size_t CategoryDatabase::size() const { return categories_.size(); }

const Metadata& CategoryDatabase::metadata() const { return metadata_; }

Metadata& CategoryDatabase::metadata() { return metadata_; }
}
