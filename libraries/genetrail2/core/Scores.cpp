#include "Scores.h"
#include "Category.h"
#include "GeneSet.h"

#include <algorithm>

namespace GeneTrail
{
	Scores::ScoresProxy::ScoresProxy(const std::vector<Score>* data)
	    : data_(data)
	{
	}

	Scores::NamesProxy::NamesProxy(const std::vector<Score>* data) : data_(data)
	{
	}

	Scores::Scores(const std::vector<Score>& data)
	    : data_(data), isSortedByName_(false)
	{
	}

	Scores::Scores(std::vector<Score>&& data)
	    : data_(std::move(data)), isSortedByName_(false)
	{
	}

	Scores::Scores(const GeneTrail::GeneSet& gene_set) : isSortedByName_(true)
	{
		data_.reserve(gene_set.size());

		// Insert the entries of the gene set. emplace_back
		// ensures, that isSortedByName_ is updated properly.
		for(const auto& entry : gene_set) {
			emplace_back(entry.first, entry.second);
		}
	}

	Scores::Scores(GeneSet&& gene_set) : isSortedByName_(true)
	{
		data_.reserve(gene_set.size());

		// Insert the entries of the gene set. emplace_back
		// ensures, that isSortedByName_ is updated properly.
		for(auto& entry : gene_set) {
			emplace_back(std::move(entry.first), entry.second);
		}
	}

	Scores::Scores(size_t size) : isSortedByName_(true) { data_.reserve(size); }

	Scores Scores::subset(const Category& c) const
	{
		if(isSortedByName_) {
			return subsetMerge_(c);
		}

		return subsetFind_(c);
	}

	Scores Scores::subsetMerge_(const Category& c) const
	{
		Scores result(std::min(size(), c.size()));

		auto scoresIt = begin();
		auto categoryIt = c.begin();

		auto predicate = [](const Score& a, const Score& b) {
			return a.name() < b.name();
		};

		while(scoresIt != end() && categoryIt != c.end()) {
			//TODO: Temporary can be optimized away using C++14 heterogenous
			//      lookup.
			Score dummy(*categoryIt, 0.0);
			scoresIt = std::lower_bound(scoresIt, end(), dummy, predicate);

			if(scoresIt->name() == *categoryIt) {
				result.emplace_back(*scoresIt);
				++scoresIt;
			}
			++categoryIt;
		}

		return result;
	}

	Scores Scores::subsetFind_(const Category& c) const
	{
		Scores result(std::min(size(), c.size()));

		for(const auto& entry : *this) {
			if(c.contains(entry.name())) {
				result.emplace_back(entry);
			}
		}

		return result;
	}

	void Scores::sortByName()
	{
		// Nothing to do here
		if(isSortedByName_) {
			return;
		}

		isSortedByName_ = true;

		std::sort(
		    data_.begin(), data_.end(),
		    [](const Score& a, const Score& b) { return a.name() < b.name(); });
	}

	void Scores::sortByScore()
	{
		// The data is only sorted by name if there is at most one item present.
		isSortedByName_ = size() <= 1;

		std::sort(data_.begin(), data_.end(),
		          [](const Score& a,
		             const Score& b) { return a.score() < b.score(); });
	}

	bool Scores::contains(const std::string& name) const
	{
		// TODO: This might be replaced with heterogenous lookup from
		//      C++14
		Score dummy(name, 0.0);
		if(isSortedByName_) {
			return std::binary_search(data_.begin(), data_.end(), dummy,
			                          [](const Score& a, const Score& s) {
				return a.name() < s.name();
			});
		} else {
			return std::find_if(data_.begin(), data_.end(), [&name](const Score& a) {return a.name() == name; }) != data_.end();
		}
	}
}
