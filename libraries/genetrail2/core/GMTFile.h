#ifndef GT2_CORE_GMTFILE_H
#define GT2_CORE_GMTFILE_H

#include "CategoryFile.h"
#include "Category.h"

#include "macros.h"

namespace GeneTrail
{
	class GT2_EXPORT GMTFile : public CategoryFile
	{
		public:
		GMTFile(const std::string& path, FileOpenMode mode = FileOpenMode::READ);

		GMTFile(GMTFile&&) = default;
		GMTFile& operator=(GMTFile&&) = default;

		/// GMTFile is not copy constructible
		/// as file streams are not copy constructible
		GMTFile(const GMTFile&) = delete;

		/// GMTFile cannot be copied
		/// as file streams cannotbe copied
		Category& operator=(const GMTFile&) = delete;

		Category read();
		bool write(const Category& cat);

		private:
		std::string next_line_;

		void advanceLine_();
	};
}

#endif //GT2_CORE_GMTFILE_H

