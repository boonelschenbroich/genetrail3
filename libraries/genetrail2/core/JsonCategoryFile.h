/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2016 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_JSONCATEGORYFILE_H
#define GT2_CORE_JSONCATEGORYFILE_H

#include "CategoryDatabaseFile.h"
#include "Category.h"

#include "macros.h"

#include <memory>
#include <string>

namespace GeneTrail
{

class GT2_EXPORT JsonCategoryFile : public CategoryDatabaseFile
{
  public:
	/**
	 * Constructor.
	 * @param db The EntityDatabase that should be used by the
	 *           read CategoryDatabase.
	 * @param path The path to the file that should be opened
	 * @param mode Should the file be used for reading or writing?
	 */
	JsonCategoryFile(const std::shared_ptr<EntityDatabase>& db,
	                 const std::string& path,
	                 FileOpenMode mode = FileOpenMode::READ);

	// Allow moving a JsonCategoryFile
	JsonCategoryFile(JsonCategoryFile&&) = default;
	JsonCategoryFile& operator=(JsonCategoryFile&&) = default;

	// Disallow copying a JsonCategoryFile
	JsonCategoryFile(const JsonCategoryFile&) = delete;
	JsonCategoryFile& operator=(const JsonCategoryFile&) = delete;

	/// Replace the used EntityDatabase
	void setEntityDatabase(const std::shared_ptr<EntityDatabase>& db);
	/// Get the used EntityDatabase
	const std::shared_ptr<EntityDatabase>& getEntityDatabase() const;

	/// Read a CategoryDatabase
	CategoryDatabase read() override;
	/// Write a CategoryDatabase to file
	bool write(const CategoryDatabase& cat) override;

  private:
	std::shared_ptr<EntityDatabase> entity_database_;
};


}

#endif // GT2_CORE_JSONCATEGORYFILE_H
