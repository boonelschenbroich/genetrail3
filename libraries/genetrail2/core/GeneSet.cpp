#include "GeneSet.h"

#include "Category.h"
#include "GeneSetFilters.h"

#include <algorithm>
#include <cmath>

namespace GeneTrail
{
	/**
	 * This struct is used as comparator for sorting
	 */
	struct increasing_compare
	{
		bool operator()(const GeneSet::Element& a, const GeneSet::Element& b)
		{
			return a.second < b.second;
		}
	};

	/**
	 * This struct is used as comparator for sorting
	 */
	struct decreasing_compare
	{
		bool operator()(const GeneSet::Element& a, const GeneSet::Element& b)
		{
			return a.second > b.second;
		}
	};

	/**
	 * This struct is used as comparator for sorting
	 */
	struct absolute_compare
	{
		bool operator()(const GeneSet::Element& a, const GeneSet::Element& b)
		{
			return std::abs(a.second) > std::abs(b.second);
		}
	};

	GeneSet::Container GeneSet::getSortedScores(bool decreasing) const
	{
		if(decreasing) {
			return getDecreasinglySortedScores();
		}

		return getIncreasinglySortedScores();
	}

	GeneSet::Container GeneSet::getIncreasinglySortedScores() const
	{
		Container sorted_scores(container_);
		std::sort(sorted_scores.begin(), sorted_scores.end(),
		          increasing_compare());
		return sorted_scores;
	}

	GeneSet::Container GeneSet::getDecreasinglySortedScores() const
	{
		Container sorted_scores(container_);
		std::sort(sorted_scores.begin(), sorted_scores.end(),
		          decreasing_compare());
		return sorted_scores;
	}

	GeneSet::Container GeneSet::getAbsoluteSortedScores() const
	{
		Container sorted_scores(container_);
		std::sort(sorted_scores.begin(), sorted_scores.end(),
		          absolute_compare());
		return sorted_scores;
	}

	GeneSet::Container
	GeneSet::intersect(const std::vector<Element>& scores,
	                   const std::set<std::string>& myset) const
	{
		Container inter;

		for(auto& elem : scores) {
			std::set<std::string>::iterator setIt;
			setIt = myset.find(elem.first);

			if(setIt != myset.end()) {
				inter.push_back(elem);
			}
		}

		return inter;
	}

	GeneSet::Container
	GeneSet::intersect(const std::set<std::string>& set) const
	{
		return intersect(container_, set);
	}

	GeneSet::Container
	GeneSet::sortAndIntersect(const std::set<std::string>& set,
	                          bool decreasing) const
	{
		return intersect(getSortedScores(decreasing), set);
	}

	GeneSet::Container GeneSet::getFirstK(const std::vector<Element>& scores,
	                                      int k) const
	{
		return GeneSet::Container{scores.begin(), scores.begin() + k};
	}

	std::vector<std::string>
	GeneSet::getIdentifier(const std::vector<GeneSet::Element>& scores) const
	{
		std::vector<std::string> s(scores.size());

		for(size_t i = 0; i < scores.size(); ++i) {
			s[i] = scores[i].first;
		}

		return s;
	}

	std::vector<std::string> GeneSet::getIdentifier() const
	{
		return getIdentifier(container_);
	}

	std::vector<std::string> GeneSet::getSortedIdentifier(bool decreasing) const
	{
		return getIdentifier(getSortedScores(decreasing));
	}

	std::vector<std::string> GeneSet::getDecreasinglySortedIdentifier() const
	{
		return getIdentifier(getSortedScores(true));
	}

	std::vector<std::string> GeneSet::getIncreasinglySortedIdentifier() const
	{
		return getIdentifier(getSortedScores(false));
	}

	std::vector<std::string> GeneSet::getAbsoluteSortedIdentifier() const
	{
		return getIdentifier(getAbsoluteSortedScores());
	}

	Category GeneSet::toCategory(const std::string& name) const
	{
		Category res(name);
		for(const auto& it : container_) {
			res.insert(it.first);
		}
		return res;
	}

	GeneSet& GeneSet::filter(GeneSetFilter::GeneSetFilter* gf)
	{
		if(gf == nullptr) {
			return *this;
		}

		// Call the filters setup routine. This is needed for
		// filters that for example need to look at the score distribution.
		gf->setup(*this);

		auto new_end =
		    std::remove_if(container_.begin(), container_.end(),
		                   [gf](const Element& e) { return gf->filter(e); });

		container_.erase(new_end, container_.end());

		return *this;
	}

	void abs(GeneSet& gene_set)
	{
		transform(gene_set, [](double d) { return std::abs(d); });
	}

	void sqrt(GeneSet& gene_set)
	{
		transform(gene_set, [](double d) { return std::sqrt(d); });
	}

	void log(GeneSet& gene_set)
	{
		transform(gene_set, [](double d) { return std::log(d); });
	}

	void log2(GeneSet& gene_set)
	{
		transform(gene_set, [](double d) { return std::log2(d); });
	}

	void log10(GeneSet& gene_set)
	{
		transform(gene_set, [](double d) { return std::log10(d); });
	}

	void pow(GeneSet& gene_set, int n)
	{
		std::transform(gene_set.begin(), gene_set.end(), gene_set.begin(),
		               [n](GeneSet::Element& e) {
			e.second = std::pow(e.second, n);
			return e;
		});
	}

	void pow2(GeneSet& gene_set) { pow(gene_set, 2); }
}

