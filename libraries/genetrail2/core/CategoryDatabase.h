/*
 * GeneTrail2 - An efficent library for interpreting genetic data
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
#ifndef GT2_CORE_CATEGORYDATABASE_H
#define GT2_CORE_CATEGORYDATABASE_H

#include "Category.h"
#include "Editor.h"
#include "Metadata.h"

#include "macros.h"

#include <memory>
#include <string>
#include <vector>

namespace GeneTrail
{

/**
 * Class representing a collection of Categories.
 */
class GT2_EXPORT CategoryDatabase
{
  public:
	using iterator = std::vector<Category>::iterator;
	using const_iterator = std::vector<Category>::const_iterator;

	/**
	 * Constructor
	 *
	 * @param db An EntityDatabase instance used for all contained categories.
	 */
	explicit CategoryDatabase(const std::shared_ptr<EntityDatabase>& db);

	/**
	 * Creates and returns a new Category.
	 * This method ensures, that the Category uses the correct EntityDatabase.
	 */
	template <typename... Args> Category& addCategory(Args... args)
	{
		categories_.emplace_back(entity_database_.get(),
		                         std::forward<Args>(args)...);

		return categories_.back();
	}

	/**
	 * Reserve storage for new Categories.
	 *
	 * @param capacity The number of elements that should at least
	 *                 fit into the CategoryDatabase.
	 */
	void reserve(size_t capacity);

	/**
	 * @returns the number of stored categories.
	 */
	size_t size() const;

	iterator begin();
	iterator end();

	const_iterator begin() const;
	const_iterator end() const;

	const_iterator cbegin() const;
	const_iterator cend() const;

	const std::string& name() const;
	void setName(const std::string& name);

	Editor& editor();
	const Editor& editor() const;
	void setEditor(const Editor& editor);

	const std::string& creationDate() const;
	void setCreationDate(const std::string& date);

	const std::string& sourceUrl() const;
	void setSourceUrl(const std::string& url);

	const std::string& identifier() const;
	void setIdentifier(const std::string& identifer);

	/**
	 * Returns the metadata attached to the database.
	 */
	const Metadata& metadata() const;

	/**
	 * Returns the metadata attached to the database.
	 */
	Metadata& metadata();

	/// Returns the used EntityDatabase instance
	const std::shared_ptr<EntityDatabase>& entityDatabase() const;

	/// Get the i-th stored category
	const Category& operator[](size_t i) const;
	/// Get the i-th stored category
	Category& operator[](size_t i);

  private:
	std::string name_;
	Editor editor_;

	std::string creation_date_;
	std::string source_url_;
	std::string identifier_;

	Metadata metadata_;
	std::shared_ptr<EntityDatabase> entity_database_;

	std::vector<Category> categories_;
};
}

#endif // GT2_CORE_CATEGORYDATABASE_H
