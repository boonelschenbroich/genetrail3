#include "Category.h"

#include <algorithm>

namespace GeneTrail
{
	Category::Category(EntityDatabase* database, const std::string& name)
	    : name_(name), database_(database)
	{
	}

	Category::Category(EntityDatabase* database)
	    : database_(database)
	{
	}

	const std::string& Category::name() const { return name_; }

	void Category::setName(const std::string& n) { name_ = n; }
	void Category::setName(std::string&& n) { name_ = std::move(n); }

	const std::string& Category::reference() const { return reference_; }

	void Category::setReference(std::string r) { reference_ = std::move(r); }

	bool Category::contains(const std::string& id) const
	{
		return contains(database_->index(id));
	}

	bool Category::contains(size_t i) const
	{
		return container_.find(i) != container_.end();
	}

	bool Category::insert(const std::string& id)
	{
		return insert(database_->index(id));
	}

	bool Category::insert(size_t i) { return container_.emplace(i).second; }

	size_t Category::size() const { return container_.size(); }

	bool Category::empty() const { return container_.empty(); }

	bool Category::operator<(const Category& o) const
	{
		return std::lexicographical_compare(
		    container_.begin(), container_.end(), o.container_.begin(),
		    o.container_.end());
	}

	bool Category::operator==(const Category& o) const
	{
		if(size() != o.size()) {
			return false;
		}

		auto it = begin();
		auto jt = o.begin();

		for(; it != end(); ++it, ++jt) {
			if(*it != *jt) {
				return false;
			}
		}

		return true;
	}

	Category Category::intersect(const std::string& name, const Category& a,
	                             const Category& b)
	{
		Category result(a.database_);

		result.setName(name);

		if(a.database_ != b.database_) {
			throw "";
		}

		std::set_intersection(
		    a.container_.begin(), a.container_.end(), b.container_.begin(),
		    b.container_.end(),
		    std::inserter(result.container_, std::end(result.container_)));

		return result;
	}

	Category Category::combine(const std::string& name, const Category& a,
	                           const Category& b)
	{
		Category result(a.database_);
		result.setName(name);

		if(a.database_ != b.database_) {
			throw "";
		}

		result.container_.reserve(std::max(a.size(), b.size()));

		std::set_union(
		    a.container_.begin(), a.container_.end(), b.container_.begin(),
		    b.container_.end(),
		    std::inserter(result.container_, std::end(result.container_)));

		return result;
	}

	std::ostream& operator<<(std::ostream& strm, const Category& cat)
	{
		strm << &cat << std::endl << std::endl;
		strm << cat.name() << '\t' << cat.reference();

		std::copy(cat.names().begin(), cat.names().end(),
		          std::ostream_iterator<std::string>(strm, "\t"));

		return strm;
	}
}
