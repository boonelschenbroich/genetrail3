#ifndef GT2_CORE_CATEGORY_H
#define GT2_CORE_CATEGORY_H

#include <forward_list>
#include <memory>
#include <string>

#include <boost/container/flat_set.hpp>

#include "macros.h"

namespace GeneTrail
{
	class GT2_EXPORT Category
	{
		private:
		typedef boost::container::flat_set<std::string> Container;

		public:
		typedef Container::iterator iterator;
		typedef Container::const_iterator const_iterator;

		static Category intersect(std::string name, const Category& a,
		                          const Category& b);
		static Category combine(std::string name, const Category& a,
		                        const Category& b);

		template <typename InputIterator>
		Category(std::string name, InputIterator first, InputIterator last)
		    : container_(first, last), name_(std::move(name))
		{
		}

		explicit Category(std::string name);
		explicit Category(std::string name,
		                  const std::shared_ptr<Category>& parent);

		Category() = default;
		Category(const Category&) = default;
		Category(Category&&) = default;

		Category& operator=(const Category&) = default;
		Category& operator=(Category&&) = default;

		//
		// Setters and Getters
		//

		const std::string& name() const;
		void setName(std::string n);

		const std::string& reference() const;
		void setReference(std::string r);

		//
		// Element access
		//
		bool contains(const std::string& id) const;
		bool insert(std::string id);

		const std::shared_ptr<Category>& getParent();

		size_t size() const;
		bool empty() const;

		//
		// Iterators
		//
		iterator begin() { return container_.begin(); }
		const_iterator begin() const { return container_.begin(); }

		iterator end() { return container_.end(); }
		const_iterator end() const { return container_.end(); }

		//
		// Operations
		//
		friend Category intersect(std::string name, const Category& a, const Category& b);
		friend Category combine(std::string name, const Category& a, const Category& b);

		private:
		Container container_;

		std::string name_;
		std::string reference_;

		std::shared_ptr<Category> parent_;
		std::forward_list<std::shared_ptr<Category>> children_;
	};
}

#endif // GT2_CORE_CATEGORY_H
