#include "EntityDatabase.h"

#include "Exception.h"

namespace GeneTrail
{
	GT2_EXPORT std::shared_ptr<EntityDatabase> EntityDatabase::global =
	    std::make_shared<EntityDatabase>();

	void EntityDatabase::clear()
	{
		name_to_index_.clear();
		db_.clear();
	}

	size_t EntityDatabase::index(const std::string& name)
	{
		auto res = name_to_index_.find(name);

		if(res == name_to_index_.end()) {
			res = name_to_index_.emplace(name, db_.size()).first;
			db_.push_back(name);
		}

		return res->second;
	}

	size_t EntityDatabase::index(const std::string& name) const
	{
		auto res = name_to_index_.find(name);

		if(res == name_to_index_.end()) {
			throw UnknownEntry(name);
		}

		return res->second;
	}
}
