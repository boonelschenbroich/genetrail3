/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
 *               2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_TEXT_FILE_H
#define GT2_CORE_TEXT_FILE_H

#include "macros.h"
#include "Exception.h"
#include "File.h"

#include <functional>
#include <vector>
#include <set>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace GeneTrail
{
	class GT2_EXPORT TextFile : public File<std::vector<std::string>>
	{
		public:

		TextFile(const std::string& path, const std::string& delimiter, const std::set<std::string>& skipSymbols = std::set<std::string>(), FileOpenMode mode = FileOpenMode::READ);

		TextFile(TextFile&&) = default;
		TextFile& operator=(TextFile&&) = default;

		/// TextFile is not copy constructible
		/// as file streams are not copy constructible
		TextFile(const TextFile&) = delete;

		std::vector<std::string> read();

		bool write(const std::vector<std::string>& r);

		private:
		std::string next_line_;
		std::string delimiter_;
		std::set<std::string> skipSymbols_;

		void advanceLine_();
	};
}

#endif //GT2_CORE_TEXT_FILE_READER_H

