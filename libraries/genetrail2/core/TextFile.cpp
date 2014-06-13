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

