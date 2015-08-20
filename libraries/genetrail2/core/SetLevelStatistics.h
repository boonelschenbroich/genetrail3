#ifndef GT2_SET_LEVEL_STATISTICS_H
#define GT2_SET_LEVEL_STATISTICS_H

#include "HTest.h"
#include "Scores.h"
#include "GeneSetEnrichmentAnalysis.h"
#include "WeightedGeneSetEnrichmentAnalysis.h"
#include "OverRepresentationAnalysis.h"
#include "OneSampleTTest.h"

#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

namespace GeneTrail
{
	class Category;
	class DenseMatrix;
	class DenseMatrixSubset;

	namespace SetLevelStatistics
	{
		struct Direct
		{
		};
		struct Indirect
		{
		};
		struct Scores
		{
		};
		struct Identifiers
		{
		};
		struct Matrix
		{
		};
	}

	class StatisticsEnrichment
	{
		public:
		using RowWiseMode = SetLevelStatistics::Indirect;
		using InputType = SetLevelStatistics::Scores;

		using _viter = Scores::ConstScoreIterator;
		using Statistics = std::function<double(_viter, _viter)>;

		StatisticsEnrichment(const Statistics& test, const Scores& scores)
		    : test_(test), scores_(scores)
		{
			scores_.sortByIndex();
			expected_value_ =
			    test_(scores_.scores().begin(), scores_.scores().end());
		}

		void setInputScores(const Scores& scores)
		{
			scores_ = scores;
			scores_.sortByIndex();
			expected_value_ =
			    test_(scores_.scores().begin(), scores_.scores().end());
		}

		bool canUseCategory(const Category&, size_t) const { return true; }

		std::tuple<double, double> computeScore(const Category& c) const
		{
			auto intersection = scores_.subset(c);
			auto score = test_(intersection.scores().begin(),
			                   intersection.scores().end());

			return std::make_tuple(score, expected_value_);
		}

		private:
		Statistics test_;
		Scores scores_;
		double expected_value_;
	};

	template <typename Test> class HTestEnrichmentBase
	{
		public:
		using RowWiseMode = SetLevelStatistics::Direct;
		using InputType = SetLevelStatistics::Scores;

		HTestEnrichmentBase(const Scores& scores) : scores_(scores) {}

		double computeRowWisePValue(Test& test, EnrichmentResult* result)
		{
			if(result->enriched) {
				return HTest::upperTailedPValue(test, result->score);
			} else {
				return HTest::lowerTailedPValue(test, result->score);
			}
		}

		protected:
		Scores scores_;
	};

	template <typename Test>
	class HTestEnrichment : public HTestEnrichmentBase<Test>
	{
		public:
		using HTestEnrichmentBase<Test>::HTestEnrichmentBase;

		void setInputScores(const Scores& scores) { this->scores_ = scores; }

		bool canUseCategory(const Category&, size_t hits) const
		{
			return hits > 2;
		}

		double computeRowWisePValue(EnrichmentResult* result) {
			return HTestEnrichmentBase<Test>::computeRowWisePValue(test_, result);
		}

		std::tuple<double, double> computeScore(const Category& c)
		{
			using namespace boost;

			auto filter_pos =
			    [&c](const Score& s) { return c.contains(s.index()); };
			auto filter_neg =
			    [&c](const Score& s) { return !c.contains(s.index()); };

			auto it_tmp_beg = make_filter_iterator(
			    filter_pos, this->scores_.begin(), this->scores_.end());
			auto it_tmp_end = make_filter_iterator(
			    filter_pos, this->scores_.end(), this->scores_.end());
			auto jt_tmp_beg = make_filter_iterator(
			    filter_neg, this->scores_.begin(), this->scores_.end());
			auto jt_tmp_end = make_filter_iterator(
			    filter_neg, this->scores_.end(), this->scores_.end());

			auto project_to_score = [](const Score& s) { return s.score(); };
			auto it_beg = make_transform_iterator(it_tmp_beg, project_to_score);
			auto it_end = make_transform_iterator(it_tmp_end, project_to_score);
			auto jt_beg = make_transform_iterator(jt_tmp_beg, project_to_score);
			auto jt_end = make_transform_iterator(jt_tmp_end, project_to_score);

			auto score = HTest::test(test_, it_beg, it_end, jt_beg, jt_end);

			return std::make_tuple(score, 0.0);
		}

		private:
		Test test_;
	};

	template <typename T>
	class HTestEnrichment<OneSampleTTest<T>>
	    : public HTestEnrichmentBase<OneSampleTTest<T>>
	{
		public:
		HTestEnrichment(const Scores& scores)
		    : HTestEnrichmentBase<OneSampleTTest<T>>(scores),
		      test_(1e-5, getExpectedValue_())
		{
		}

		void setInputScores(const Scores& scores)
		{
			this->scores_ = scores;
			//TODO: update test here
		}

		bool canUseCategory(const Category&, size_t hits) const
		{
			return hits > 1;
		}

		double computeRowWisePValue(EnrichmentResult* result) {
			return HTestEnrichmentBase<OneSampleTTest<T>>::computeRowWisePValue(test_, result);
		}

		std::tuple<double, double> computeScore(const Category& c)
		{
			using namespace boost;

			auto filter_pos =
			    [&c](const Score& s) { return c.contains(s.index()); };

			auto it_tmp_beg = make_filter_iterator(
			    filter_pos, this->scores_.begin(), this->scores_.end());
			auto it_tmp_end = make_filter_iterator(
			    filter_pos, this->scores_.end(), this->scores_.end());

			auto project_to_score = [](const Score& s) { return s.score(); };
			auto it_beg = make_transform_iterator(it_tmp_beg, project_to_score);
			auto it_end = make_transform_iterator(it_tmp_end, project_to_score);

			auto score = HTest::test(this->test_, it_beg, it_end);

			return std::make_tuple(score, 0.0);
		}

		private:
		T getExpectedValue_()
		{
			return statistic::mean<T>(this->scores_.scores().begin(),
			                          this->scores_.scores().end());
		}

		private:
		OneSampleTTest<T> test_;
	};

	class Ora
	{
		public:
		using RowWiseMode = SetLevelStatistics::Direct;
		using InputType = SetLevelStatistics::Identifiers;

		Ora(const Category& reference_set, const Category& test_set)
		    : test_(reference_set, test_set){};

		bool canUseCategory(const Category&, size_t) const { return true; }

		std::tuple<double, double> computeScore(const Category& c) const
		{
			auto score = test_.computeScore(c);
			return std::make_tuple(score, test_.expectedNumberOfHits(c));
		}

		double computeRowWisePValue(EnrichmentResult* result) const
		{
			return test_.computePValue(*result->category);
		}

		private:
		OverRepresentationAnalysis test_;
	};

	class KolmogorovSmirnov
	{
		public:
		using RowWiseMode = SetLevelStatistics::Direct;
		using InputType = SetLevelStatistics::Scores;

		template <typename Iterator>
		KolmogorovSmirnov(const Iterator& beginIds, const Iterator& endIds)
		    : ids_(beginIds, endIds)
		{
		}

		void setInputScores(Scores scores)
		{
			scores.sortByScore();
			ids_.assign(scores.names().begin(), scores.names().end());
		}

		bool canUseCategory(const Category&, size_t) const { return true; }

		std::tuple<double, double> computeScore(const Category& category)
		{
			auto score =
			    test_.computeRunningSum(category, ids_.begin(), ids_.end());
			return std::make_tuple(score, 0.0);
		}

		double computeRowWisePValue(EnrichmentResult* result)
		{
			auto intersection_size = test_.intersectionSize(
			    *result->category, ids_.begin(), ids_.end());

			if(result->score > 0.0) {
				return test_.computeRightPValue(ids_.size(), intersection_size,
				                                result->score)
				    .convert_to<double>();
			} else {
				return test_.computeLeftPValue(ids_.size(), intersection_size,
				                               result->score)
				    .convert_to<double>();
			}
		}

		private:
		std::vector<std::string> ids_;
		GeneSetEnrichmentAnalysis<big_float, int64_t> test_;
	};

	class WeightedKolmogorovSmirnov
	{
		public:
		using RowWiseMode = SetLevelStatistics::Indirect;
		using InputType = SetLevelStatistics::Scores;

		WeightedKolmogorovSmirnov(const Scores& scores) : test_(scores) {}

		void setInputScores(const Scores& scores) { test_.setScores(scores); }

		bool canUseCategory(const Category&, size_t) const { return true; }

		std::tuple<double, double> computeScore(const Category& category)
		{
			auto score = test_.computeRunningSum(category);

			return std::make_tuple(score, 0.0);
		}

		private:
		WeightedGeneSetEnrichmentAnalysis<double> test_;
	};
}

#endif // GT2_SET_LEVEL_STATISTICS_H
