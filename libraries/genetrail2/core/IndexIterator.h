#ifndef INDEX_ITERATOR_H
#define INDEX_ITERATOR_H

#include <iterator>

#include "macros.h"

namespace GeneTrail
{
	template <typename Value_t, typename Container_t>
	class GT2_EXPORT IndexIterator : public std::iterator<std::random_access_iterator_tag, Value_t>
	{
		protected:

		Container_t* container_;
		int index_;

		public:

		IndexIterator() : container_(0), index_(0)
		{
		}

		IndexIterator(Container_t& container, int index)
		    : container_(&container), index_(index)
		{
		}

		bool operator==(const IndexIterator& other)
		{
			return container_ == other.container_ && index_ == other.index_;
		}

		bool operator!=(const IndexIterator& other)
		{
			return !(*this == other);
		}

		Value_t& operator*()
		{
			return (*container_)[index_];
		}

		Value_t const& operator*() const
		{
			return (*container_)[index_];
		}

		Value_t* operator->()
		{
			return &((*container_)[index_]);
		}

		Value_t const* operator->() const
		{
			return &((*container_)[index_]);
		}

		IndexIterator& operator++()
		{
			++index_;
			return *this;
		}

		IndexIterator operator++(int)
		{
			IndexIterator prev(*this);
			operator++();
			return prev;
		}

		IndexIterator& operator--()
		{
			--index_;
			return *this;
		}

		IndexIterator operator--(int)
		{
			IndexIterator prev(*this);
			operator--();
			return prev;
		}

		friend IndexIterator operator+(const IndexIterator& a, int b)
		{
			IndexIterator ret(a);
			ret += b;
			return ret;
		}

		friend IndexIterator operator-(const IndexIterator& a, int b)
		{
			IndexIterator ret(a);
			ret -= b;
			return ret;
		}

		friend IndexIterator operator+(int a, const IndexIterator& b)
		{
			IndexIterator ret(b);
			ret += a;
			return ret;
		}

		friend IndexIterator operator-(int a, const IndexIterator& b)
		{
			IndexIterator ret(b);
			ret -= a;
			return ret;
		}

		int operator-(const IndexIterator& other) const
		{
			return index_ - other.index_;
		}

		bool operator<(const IndexIterator& other)
		{
			return container_ == other.container_ && index_ < other.index_;
		}

		bool operator<=(const IndexIterator& other)
		{
			return container_ == other.container_ && index_ <= other.index_;
		}

		bool operator>(const IndexIterator& other)
		{
			return container_ == other.container_ && index_ > other.index_;
		}

		bool operator>=(const IndexIterator& other)
		{
			return container_ == other.container_ && index_ >= other.index_;
		}

		IndexIterator& operator+=(int b)
		{
			index_ += b;
		}

		IndexIterator& operator-=(int b)
		{
			index_ -= b;
		}

		Value_t& operator[](int i)
		{
			return (*container_)[i];
		}

		Value_t const& operator[](int i) const
		{
			return (*container_)[i];
		}
	};

	template <typename Value_t, typename Container_t>
	inline IndexIterator<Value_t, Container_t>
	begin(Container_t& container)
	{
		return IndexIterator<Value_t, Container_t>(container, 0);
	}

	template <typename Value_t, typename Container_t>
	inline IndexIterator<Value_t, Container_t>
	end(Container_t& container)
	{
		return IndexIterator<Value_t, Container_t>(container,
		                                            container.size());
	}

	template <typename Value_t, typename Container_t>
	inline IndexIterator<const Value_t, const Container_t>
	begin(const Container_t& container)
	{
		return IndexIterator<const Value_t, const Container_t>(container, 0);
	}

	template <typename Value_t, typename Container_t>
	inline IndexIterator<const Value_t, const Container_t>
	end(const Container_t& container)
	{
		return IndexIterator<const Value_t, const Container_t>(
		    container, container.size());
	}
}

#endif //INDEX_ITERATOR_H

