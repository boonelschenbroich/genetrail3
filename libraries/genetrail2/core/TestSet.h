#ifndef GT2_CORE_TESTSET_H
#define GT2_CORE_TESTSET_H

#include "macros.h"
#include "Category.h"

#include <utility>
#include <vector>

namespace GeneTrail
{

	class GT2_EXPORT TestSet
	{
		private:
		typedef std::vector<std::pair<std::string, double>> Container;

		public:
		typedef Container::iterator iterator;
		typedef Container::const_iterator const_iterator;
		typedef std::pair<std::string, double> Element;

		static TestSet readSet(std::istream& input);
		static TestSet readSetWithoutScores(std::istream& input);

		TestSet() = default;
		TestSet(const TestSet&) = default;
		TestSet(TestSet&&) = default;

		TestSet& operator=(const TestSet&) = default;
		TestSet& operator=(TestSet&&) = default;

		iterator begin();
		iterator end();

		const_iterator begin() const;
		const_iterator end() const;

		const Element& operator[](size_t i) const;
		Element& operator[](size_t i);

		Category toCategory(const std::string& name = "") const;

		private:
		std::vector<Element> container_;
	};
}

#endif // GT2_CORE_TESTSET_H

