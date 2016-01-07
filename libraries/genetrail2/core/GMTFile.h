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
#ifndef GT2_CORE_GMTFILE_H
#define GT2_CORE_GMTFILE_H

#include "CategoryDatabaseFile.h"
#include "CategoryDatabase.h"

#include "macros.h"

namespace GeneTrail
{
	class GT2_EXPORT GMTFile : public CategoryDatabaseFile
	{
	  public:
	    GMTFile(const std::shared_ptr<EntityDatabase>& db,
	            const std::string& path,
	            FileOpenMode mode = FileOpenMode::READ);

	    GMTFile(GMTFile&&) = default;
		GMTFile& operator=(GMTFile&&) = default;

		/// GMTFile is not copy constructible
		/// as file streams are not copy constructible
		GMTFile(const GMTFile&) = delete;

		/// GMTFile cannot be copied
		/// as file streams cannotbe copied
		GMTFile& operator=(const GMTFile&) = delete;

		CategoryDatabase read() override;
		bool write(const CategoryDatabase& db) override;

	  private:
		void advanceLine_();
		void readCategory_(CategoryDatabase& db);

		std::string next_line_;
		std::shared_ptr<EntityDatabase> entity_database_;
	};
}

#endif //GT2_CORE_GMTFILE_H

