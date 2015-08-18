#include "Scores.h"
#include "Category.h"
#include "GeneSet.h"

#include <algorithm>

namespace GeneTrail
{
	Scores::IndexProxy::IndexProxy(const std::vector<Score>* data) : data_(data)
	{
	}

	Scores::ScoresProxy::ScoresProxy(const std::vector<Score>* data)
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
		Scores result(std::min(size(), c.size()), db_);

		auto scoresIt = begin();
		auto categoryIt = c.begin();

		auto predicate = [](const Score& a,
		                    const Score& b) { return a.index() < b.index(); };

		while(scoresIt != end() && categoryIt != c.end()) {
			// TODO: Temporary can be optimized away using C++14 heterogenous
			//      lookup.
			Score dummy(*categoryIt, 0.0);
			scoresIt = std::lower_bound(scoresIt, end(), dummy, predicate);

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
		Scores result(std::min(size(), c.size()));

		for(const auto& entry : *this) {
			if(c.contains(entry.index())) {
				result.emplace_back(entry);
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

	void Scores::sortByScore()
	{
		// The data is only sorted by name if there is at most one item present.
		isSortedByIndex_ = size() <= 1;

		std::sort(data_.begin(), data_.end(),
		          [](const Score& a,
		             const Score& b) { return a.score() < b.score(); });
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
