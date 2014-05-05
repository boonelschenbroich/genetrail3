#include "TestSet.h"

namespace GeneTrail
{

	TestSet TestSet::readSet(std::istream& input)
	{
	}

	TestSet TestSet::readSetWithoutScores(std::istream& input)
	{
	}

	TestSet::iterator TestSet::begin()
	{
		return container_.begin();
	}

	TestSet::const_iterator TestSet::begin() const
	{
		return container_.begin();
	}

	TestSet::iterator TestSet::end()
	{
		return container_.end();
	}

	TestSet::const_iterator TestSet::end() const
	{
		return container_.end();
	}

	TestSet::Element& TestSet::operator[](size_t i)
	{
		return container_[i];
	}

	const TestSet::Element& TestSet::operator[](size_t i) const
	{
		return container_[i];
	}

	Category TestSet::toCategory(const std::string& name) const
	{
		Category res(name);

		for(const auto& it : container_) {
			res.insert(it.first);
		}

		return res;
	}
}
