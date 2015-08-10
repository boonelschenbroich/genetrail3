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
	}

	class StatisticsEnrichment
	{
		public:
		using RowWiseMode = SetLevelStatistics::Indirect;

		using _viter = Scores::ConstScoreIterator;
		using Statistics = std::function<double(_viter, _viter)>;

		StatisticsEnrichment(const Statistics& test, const Scores& scores)
		    : test_(test), scores_(scores)
		{
		}

		bool canUseCategory(const Category&, size_t) const { return true; }

		double computeScore(const Category& c) const
		{
			auto intersection = scores_.subset(c);
			return test_(intersection.scores().begin(),
			             intersection.scores().end());
		}

		private:
		Statistics test_;
		const Scores& scores_;
	};

	template <typename Test> class HTestEnrichmentBase
	{
		public:
		using RowWiseMode = SetLevelStatistics::Direct;
		HTestEnrichmentBase(const Scores& scores) : scores_(scores) {}

		double computeRowWisePValue(EnrichmentResult* result)
		{
			if(result->score > 0.0) {
				return HTest::upperTailedPValue(test_, result->score);
			} else {
				return HTest::lowerTailedPValue(test_, result->score);
			}
		}

		protected:
		Scores scores_;
		Test test_;
	};

	template <typename Test>
	class HTestEnrichment : public HTestEnrichmentBase<Test>
	{
		public:
		using HTestEnrichmentBase<Test>::HTestEnrichmentBase;

		bool canUseCategory(const Category&, size_t hits) const
		{
			return hits > 2;
		}

		double computeScore(const Category& c)
		{
			using namespace boost;

			auto filter_pos =
			    [&c](const Score& s) { return c.contains(s.name()); };
			auto filter_neg =
			    [&c](const Score& s) { return !c.contains(s.name()); };

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

			return HTest::test(this->test_, it_beg, it_end, jt_beg, jt_end);
		}
	};

	template <typename T>
	class HTestEnrichment<OneSampleTTest<T>>
	    : public HTestEnrichmentBase<OneSampleTTest<T>>
	{
		public:
		using HTestEnrichmentBase<OneSampleTTest<T>>::HTestEnrichmentBase;

		bool canUseCategory(const Category&, size_t hits) const
		{
			return hits > 1;
		}

		double computeScore(const Category& c)
		{
			using namespace boost;

			auto filter_pos =
			    [&c](const Score& s) { return c.contains(s.name()); };

			auto it_tmp_beg = make_filter_iterator(
			    filter_pos, this->scores_.begin(), this->scores_.end());
			auto it_tmp_end = make_filter_iterator(
			    filter_pos, this->scores_.end(), this->scores_.end());

			auto project_to_score = [](const Score& s) { return s.score(); };
			auto it_beg = make_transform_iterator(it_tmp_beg, project_to_score);
			auto it_end = make_transform_iterator(it_tmp_end, project_to_score);

			return HTest::test(this->test_, it_beg, it_end);
		}
	};

	class Ora
	{
		public:
		using RowWiseMode = SetLevelStatistics::Direct;

		Ora(const Category& reference_set, const Category& test_set)
		    : test_(reference_set, test_set){};

		bool canUseCategory(const Category&, size_t) const { return true; }

		double computeScore(const Category& c) const
		{
			return test_.computeScore(c);
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

		template <typename Iterator>
		KolmogorovSmirnov(const Iterator& beginIds, const Iterator& endIds)
		    : intersectionSize_(0), ids_(beginIds, endIds)
		{
		}

		bool canUseCategory(const Category&, size_t) const { return true; }

		double computeScore(const Category& category)
		{
			intersectionSize_ =
			    test_.intersectionSize(category, ids_.begin(), ids_.end());
			return test_.computeRunningSum(category, ids_.begin(), ids_.end());
		}

		double computeRowWisePValue(EnrichmentResult* result)
		{
			if(result->score > 0.0) {
				return test_.computeRightPValue(intersectionSize_, ids_.size(),
				                                result->score)
				    .convert_to<double>();
			} else {
				return test_.computeLeftPValue(intersectionSize_, ids_.size(),
				                               result->score)
				    .convert_to<double>();
			}
		}

		private:
		size_t intersectionSize_;
		std::vector<std::string> ids_;
		GeneSetEnrichmentAnalysis<big_float, int64_t> test_;
	};

	class WeightedKolmogorovSmirnov
	{
		public:
		using RowWiseMode = SetLevelStatistics::Indirect;

		WeightedKolmogorovSmirnov(const Scores& scores)
		    : test_(scores)
		{
		}

		bool canUseCategory(const Category&, size_t) const { return true; }

		double computeScore(const Category& category)
		{
			return test_.computeRunningSum(category);
		}

		private:
		WeightedGeneSetEnrichmentAnalysis<double> test_;
	};
}

#endif // GT2_SET_LEVEL_STATISTICS_H
