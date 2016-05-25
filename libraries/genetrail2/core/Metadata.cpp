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
#include "Metadata.h"
#include "Exception.h"

#include "compat.h"

namespace GeneTrail
{
Metadata::Metadata(const Metadata& metadata)
{
	if(metadata.data_ != nullptr) {
		data_ = std::make_unique<Map>(*metadata.data_);
	}
}

Metadata& Metadata::operator=(const Metadata& metadata)
{
	if(this == &metadata) {
		return *this;
	}

	if(metadata.data_ != nullptr) {
		data_ = std::make_unique<Map>(*metadata.data_);
	} else {
		data_ = nullptr;
	}

	return *this;
}

bool Metadata::has(const std::string& key) const
{
	if(data_ == nullptr) {
		return false;
	}

	return data_->find(key) != data_->end();
}

const Metadata::Value& Metadata::operator[](const std::string& key) const
{
	return get(key);
}

const Metadata::Value& Metadata::get(const std::string& key) const
{
	if(data_ == nullptr) {
		throw InvalidKey(key);
	}

	auto res = data_->find(key);
	if(res != data_->end()) {
		return res->second;
	}

	throw InvalidKey(key);
}

Metadata::Value& Metadata::operator[](const std::string& key)
{
	ensureData_();

	return (*data_)[key];
}

void Metadata::ensureData_()
{
	if(data_ == nullptr) {
		data_ = std::make_unique<Map>();
	}
}

bool Metadata::remove(const std::string& key)
{
	if(data_ == nullptr) {
		return false;
	}

	return data_->erase(key) > 0;
}
}
