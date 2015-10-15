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

