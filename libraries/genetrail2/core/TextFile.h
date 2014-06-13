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

		TextFile(const std::string& path, const std::string& delimiter, const std::set<std::string>& skipSymbols, FileOpenMode mode = FileOpenMode::READ);

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

