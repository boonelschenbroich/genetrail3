#ifndef GT2_SCORES_H
#define GT2_SCORES_H

#include "macros.h"

#include <boost/iterator/transform_iterator.hpp>

#include <string>
#include <vector>

namespace GeneTrail
{
	class Category;
	class GeneSet;

	class GT2_EXPORT Score
	{
		public:
		Score(std::string&& n, double s) : name_(std::move(n)), score_(s) {}
		Score(const std::string& n, double s) : name_(n), score_(s) {}
		const std::string& name() const { return name_; }
		double score() const { return score_; }

		private:
		std::string name_;
		double score_;
	};

	class GT2_EXPORT Scores
	{
		public:
		using iterator = std::vector<Score>::iterator;
		using const_iterator = std::vector<Score>::const_iterator;

		class NamesProxy
		{
			private:
				struct ExtractName {
					const std::string& operator()(const Score& s) const {
						return s.name();
					};
				};
			const std::vector<Score>* data_;

			public:
			using const_iterator = boost::transform_iterator<ExtractName, Scores::const_iterator>;

			NamesProxy(const std::vector<Score>* data);

			const_iterator begin() const
			{
				return boost::make_transform_iterator(data_->begin(), ExtractName());
			}
			const_iterator end() const
			{
				return boost::make_transform_iterator(data_->end(), ExtractName());
			}
		};

		class ScoresProxy
		{
			private:
				struct ExtractScore{
					double operator()(const Score& s) const {
						return s.score();
					}
				};
			const std::vector<Score>* data_;

			public:
			using const_iterator = boost::transform_iterator<ExtractScore, Scores::const_iterator>;

			ScoresProxy(const std::vector<Score>* data);

			const_iterator begin() const
			{
				return boost::make_transform_iterator(data_->begin(), ExtractScore());
			}

			const_iterator end() const
			{
				return boost::make_transform_iterator(data_->end(), ExtractScore());
			}
		};

		using ConstScoreIterator = ScoresProxy::const_iterator;
		using ConstNameIterator = NamesProxy::const_iterator;

		explicit Scores(size_t size = 0);
		explicit Scores(std::vector<Score>&& data);
		explicit Scores(const std::vector<Score>& data);
		explicit Scores(const GeneSet& gene_set);
		explicit Scores(GeneSet&& gene_set);

		Scores(const Scores& o) = default;
		Scores(Scores&& o) = default;
		~Scores() = default;

		Scores& operator=(const Scores&) = default;
		Scores& operator=(Scores&&) = default;

		template <typename... Ts> void emplace_back(Ts&&... ts)
		{
			data_.emplace_back(std::forward<Ts>(ts)...);

			isSortedByName_ = size() <= 1 || (isSortedByName_ && data_[size() - 1].name() >= data_[size() - 2].name());
		}

		iterator begin() { return data_.begin(); }
		iterator end() { return data_.end(); }

		const_iterator begin() const
		{
			return data_.begin();
		}

		const_iterator end() const { return data_.end(); }

		NamesProxy names() const { return NamesProxy(&data_); };
		ScoresProxy scores() const { return ScoresProxy(&data_); };
		size_t size() const { return data_.size(); }

		Scores subset(const Category& c) const;

		const Score& operator[](size_t i) const { return data_[i]; }

		void sortByName();
		void sortByScore();
		bool isSortedByName() { return isSortedByName_; }
		bool contains(const std::string& name) const;

		private:
		Scores subsetMerge_(const Category& c) const;
		Scores subsetFind_(const Category& c) const;

		std::vector<Score> data_;
		bool isSortedByName_;
	};
}

#endif // GT2_SCORES_H
