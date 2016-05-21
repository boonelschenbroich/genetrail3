/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public
 * License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 */
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

	Category GeneSet::toCategory(const std::string& name, EntityDatabase* db) const
	{
		Category res(db);
		res.setName(name);
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

	void Overloads::abs(GeneSet& gene_set)
	{
		transform(gene_set, [](double d) { return std::abs(d); });
	}

	void Overloads::sqrt(GeneSet& gene_set)
	{
		transform(gene_set, [](double d) { return std::sqrt(d); });
	}

	void Overloads::log(GeneSet& gene_set)
	{
		transform(gene_set, [](double d) { return std::log(d); });
	}

	void Overloads::log2(GeneSet& gene_set)
	{
		transform(gene_set, [](double d) { return std::log2(d); });
	}

	void Overloads::log10(GeneSet& gene_set)
	{
		transform(gene_set, [](double d) { return std::log10(d); });
	}

	void Overloads::pow(GeneSet& gene_set, int n)
	{
		std::transform(gene_set.begin(), gene_set.end(), gene_set.begin(),
		               [n](GeneSet::Element& e) {
			e.second = std::pow(e.second, n);
			return e;
		});
	}

	void Overloads::pow2(GeneSet& gene_set) { pow(gene_set, 2); }
}

