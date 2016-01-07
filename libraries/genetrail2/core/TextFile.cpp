/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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
#include "TextFile.h"

using namespace boost;

namespace GeneTrail
{
	TextFile::TextFile(const std::string& path, const std::string& delimiter, const std::set<std::string>& skipSymbols, FileOpenMode mode)
		: File<std::vector<std::string>>(path, mode),delimiter_(delimiter), skipSymbols_(skipSymbols)
	{
		if(mode == FileOpenMode::READ) {
			advanceLine_();
		}
	}

	std::vector<std::string> TextFile::read()
	{
		if(!isValid_() || !isReading()) {
			throw IOError("File is not open for reading");
		}

		for(const auto& symbol : skipSymbols_)
		{
			if(boost::starts_with(next_line_, symbol))
			{
				advanceLine_();
				return read();
			}
		}

		std::vector<std::string> sline;
		boost::split(sline, next_line_, boost::is_any_of(delimiter_));
		for(auto& s : sline)
		{
			boost::trim(s);
		}

		// Prefetch the next line from the file.
		advanceLine_();

		return sline;
	}

	void TextFile::advanceLine_()
	{
		do {
			std::getline(*in_strm_, next_line_);
			trim(next_line_);
		} while(isValid_() && next_line_ == "");
	}

	bool TextFile::write(const std::vector<std::string>& r)
	{
		//TODO
		return false;
	}
}

