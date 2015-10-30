#include "Scores.h"
#include "Category.h"
#include "GeneSet.h"

#include <algorithm>

namespace GeneTrail
{
	Scores::IndexProxy::IndexProxy(const std::vector<Score>* data) : data_(data)
	{
	}

	Scores::ConstScoresProxy::ConstScoresProxy(const std::vector<Score>* data)
	    : data_(data)
	{
	}

	Scores::ScoresProxy::ScoresProxy(std::vector<Score>* data)
	    : data_(data)
	{
	}

	Scores::NamesProxy::NamesProxy(const std::vector<Score>* data,
	                               const EntityDatabase* db)
	    : data_(data), db_(db)
	{
	}

	Scores::Scores(const std::vector<Score>& data,
	               const std::shared_ptr<EntityDatabase>& db)
	    : data_(data), isSortedByIndex_(false), db_(db)
	{
	}

	Scores::Scores(std::vector<Score>&& data,
	               const std::shared_ptr<EntityDatabase>& db)
	    : data_(std::move(data)), isSortedByIndex_(false), db_(db)
	{
	}

	Scores::Scores(const GeneTrail::GeneSet& gene_set,
	               const std::shared_ptr<EntityDatabase>& db)
	    : isSortedByIndex_(true), db_(db)
	{
		data_.reserve(gene_set.size());

		// Insert the entries of the gene set. emplace_back
		// ensures, that isSortedByName_ is updated properly.
		for(const auto& entry : gene_set) {
			emplace_back(entry.first, entry.second);
		}
	}

	Scores::Scores(GeneSet&& gene_set,
	               const std::shared_ptr<EntityDatabase>& db)
	    : isSortedByIndex_(true), db_(db)
	{
		data_.reserve(gene_set.size());

		// Insert the entries of the gene set. emplace_back
		// ensures, that isSortedByName_ is updated properly.
		for(auto& entry : gene_set) {
			emplace_back(std::move(entry.first), entry.second);
		}
	}

	Scores::Scores(const std::shared_ptr<EntityDatabase>& db)
	    : isSortedByIndex_(true), db_(db)
	{
	}
	Scores::Scores(size_t size, const std::shared_ptr<EntityDatabase>& db)
	    : isSortedByIndex_(true), db_(db)
	{
		data_.reserve(size);
	}

	Scores Scores::subset(const Category& c) const
	{
		if(isSortedByIndex_) {
			return subsetMerge_(c);
		}

		return subsetFind_(c);
	}

	Scores Scores::subsetMerge_(const Category& c) const
	{
		auto n = size();
		Scores result(std::min(n, c.size()), db_);

		auto scoresIt = begin();
		auto categoryIt = c.begin();

		auto predicate = [](const Score& a,
		                    const Score& b) { return a.index() < b.index(); };

		while(scoresIt != end() && categoryIt != c.end()) {
			auto search_end = end();

			// This is a heuristic that tries to guess the position
			// of the current category element in the scores vector
			// For this purpose we assume that the scores vector contains
			// almost all identifiers in the entity database and we can simply
			// perform a lookup in the scores vector to get an approximate
			// location
			if(*categoryIt < n) {
				auto cat_lookup = begin() + *categoryIt;
				auto scores_index = cat_lookup->index();

				if(scores_index == *categoryIt) {
					// If we hit what we were looking for, we
					// are done and can continue.
					scoresIt = cat_lookup + 1;
					++categoryIt;
					result.emplace_back(*cat_lookup);
					continue;
				} else if(scores_index > *categoryIt) {
					// If the found index is larger than what we were looking
					// for move the end of the search range here.
					search_end = cat_lookup;
				} else {
					// Otherwise the searched element must be somewhere behind
					// the examined position.
					scoresIt = cat_lookup;
				}
			}

			// TODO: Temporary can be optimized away using C++14 heterogenous
			//      lookup.
			Score dummy(*categoryIt, 0.0);
			scoresIt = std::lower_bound(scoresIt, search_end, dummy, predicate);

			if(scoresIt->index() == *categoryIt) {
				result.emplace_back(*scoresIt);
				++scoresIt;
			}
			++categoryIt;
		}

		return result;
	}

	Scores Scores::subsetFind_(const Category& c) const
	{
		assert(c.db_.get() == c.entityDatabase());

		Scores result(std::min(size(), c.size()), db_);

		for(const auto& entry : *this) {
			if(c.contains(entry.index())) {
				result.emplace_back(entry);
			}
		}

		return result;
	}

	std::vector<size_t> Scores::subsetIndices(const Category& c) const
	{
		std::vector<size_t> result;
		result.reserve(std::min(size(), c.size()));

		for(size_t i = 0; i < size(); ++i) {
			if(c.contains(data_[i].index())) {
				result.emplace_back(i);
			}
		}

		return result;
	}

	void Scores::sortByIndex()
	{
		// Nothing to do here
		if(isSortedByIndex_) {
			return;
		}

		isSortedByIndex_ = true;

		std::sort(data_.begin(), data_.end(),
		          [](const Score& a,
		             const Score& b) { return a.index() < b.index(); });
	}

	void Scores::sortByName()
	{
		isSortedByIndex_ = size() <= 1;

		std::sort(data_.begin(), data_.end(),
		          [this](const Score& a, const Score& b) {
			return a.name(*db_) < b.name(*db_);
		});
	}

	void Scores::sortByScore(Order order)
	{
		// The data is only sorted by name if there is at most one item present.
		isSortedByIndex_ = size() <= 1;

		auto inc = [](const Score& a,
		              const Score& b) { return a.score() < b.score(); };

		auto dec = [](const Score& a,
		              const Score& b) { return a.score() > b.score(); };

		switch(order) {
			case Order::Increasing:
				std::sort(data_.begin(), data_.end(), inc);
				break;
			case Order::Decreasing:
				std::sort(data_.begin(), data_.end(), dec);
				break;
		}
	}

	bool Scores::contains(const std::string& name) const
	{
		// TODO: This might be replaced with heterogenous lookup from
		//      C++14
		Score dummy(*db_, name, 0.0);
		if(isSortedByIndex_) {
			return std::binary_search(data_.begin(), data_.end(), dummy,
			                          [](const Score& a, const Score& s) {
				return a.index() < s.index();
			});
		} else {
			return std::find_if(data_.begin(), data_.end(),
			                    [&dummy](const Score& a) {
				       return a.index() == dummy.index();
				   }) != data_.end();
		}
	}
}
