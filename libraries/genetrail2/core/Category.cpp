#include "Category.h"

#include <algorithm>

namespace GeneTrail
{
	Category::Category(std::string name) : name_(std::move(name))
	{
	}

	Category::Category(std::string name,
	                   const std::shared_ptr<Category>& parent)
	    : name_(std::move(name)), parent_(parent)
	{
	}

	const std::string& Category::name() const
	{
		return name_;
	}

	void Category::setName(std::string n)
	{
		name_ = std::move(n);
	}

	const std::string& Category::reference() const
	{
		return reference_;
	}

	void Category::setReference(std::string r)
	{
		reference_ = std::move(r);
	}

	bool Category::contains(const std::string& id) const
	{
		return container_.find(id) != container_.end();
	}

	bool Category::insert(std::string id)
	{
		return container_.emplace(std::move(id)).second;
	}

	const std::shared_ptr<Category>& Category::getParent()
	{
		return parent_;
	}

	size_t Category::size() const
	{
		return container_.size();
	}

	bool Category::empty() const
	{
		return container_.empty();
	}

	Category Category::intersect(std::string name, const Category& a, const Category& b)
	{
		Category result(std::move(name));

		std::set_intersection(
		    a.container_.begin(), a.container_.end(), b.container_.begin(),
		    b.container_.end(),
		    std::inserter(result.container_, std::end(result.container_)));

		return result;
	}

	Category Category::combine(std::string name, const Category& a, const Category& b)
	{
		Category result(std::move(name));
		result.container_.reserve(std::max(a.size(),b.size()));

		std::set_union(
		    a.container_.begin(), a.container_.end(), b.container_.begin(),
		    b.container_.end(),
		    std::inserter(result.container_, std::end(result.container_)));

		return result;
	}
}

