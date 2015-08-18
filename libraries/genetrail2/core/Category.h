#ifndef GT2_CORE_CATEGORY_H
#define GT2_CORE_CATEGORY_H

#include <forward_list>
#include <memory>
#include <string>

#include <boost/container/flat_set.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "macros.h"

#include "EntityDatabase.h"

namespace GeneTrail
{
	class GT2_EXPORT Category
	{
		private:
		typedef boost::container::flat_set<size_t> Container;

		public:
		typedef Container::iterator iterator;
		typedef Container::const_iterator const_iterator;
		class NamesProxy
		{
			public:
			using const_iterator = boost::transform_iterator<
			    std::reference_wrapper<const EntityDatabase>, Category::const_iterator>;

			NamesProxy(const Category& c, const EntityDatabase* db)
			    : cat_(c), db_(db)
			{
			}

			const_iterator begin() const
			{
				return boost::make_transform_iterator(cat_.begin(), std::ref(*db_));
			}

			const_iterator end() const
			{
				return boost::make_transform_iterator(cat_.end(), std::ref(*db_));
			}

			private:
			const Category& cat_;
			const EntityDatabase* db_;
		};

		static Category intersect(std::string name, const Category& a,
		                          const Category& b);
		static Category combine(std::string name, const Category& a,
		                        const Category& b);

		template <typename InputIterator>
		Category(InputIterator first, InputIterator last)
			: container_(first, last), database_(EntityDatabase::global)
		{
		}

		template <typename InputIterator>
		Category(std::string name, InputIterator first, InputIterator last)
		    : name_(std::move(name)), database_(EntityDatabase::global)
		{
			container_.reserve(std::distance(first, last));
			database_->transform(
			    first, last, std::inserter(container_, std::end(container_)));
		}

		explicit Category(std::string name);
		explicit Category(std::string name,
		                  const std::shared_ptr<Category>& parent);

		Category() : database_(EntityDatabase::global){}
		Category(const Category&) = default;
		Category(Category&&) = default;

		Category& operator=(const Category&) = default;
		Category& operator=(Category&&) = default;

		NamesProxy names() const { return NamesProxy(*this, database_.get()); }

		//
		// Setters and Getters
		//

		const std::string& name() const;
		void setName(std::string n);

		const std::string& reference() const;
		void setReference(std::string r);

		template<typename InputIterator>
		void replaceAll(InputIterator begin, InputIterator end) {
			container_.clear();
			container_.reserve(std::distance(begin, end));
			container_.insert(begin, end);
		}

		//
		// Element access
		//
		bool contains(const std::string& id) const;
		bool contains(size_t i) const;
		bool insert(const std::string& id);
		bool insert(size_t i);

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

		bool operator<(const Category& o) const;

		bool operator==(const Category& o) const;

		//
		// Operations
		//
		friend Category intersect(std::string name, const Category& a,
		                          const Category& b);
		friend Category combine(std::string name, const Category& a,
		                        const Category& b);

		friend std::ostream& operator<<(std::ostream& strm,
		                                const Category& cat);

		private:
		Container container_;

		std::string name_;
		std::string reference_;

		std::shared_ptr<EntityDatabase> database_;

		std::shared_ptr<Category> parent_;
		std::forward_list<std::shared_ptr<Category>> children_;
	};

	GT2_EXPORT std::ostream& operator<<(std::ostream& strm,
	                                    const Category& cat);
}

#endif // GT2_CORE_CATEGORY_H
