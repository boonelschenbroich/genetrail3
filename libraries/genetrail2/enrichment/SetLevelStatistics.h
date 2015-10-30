/*
 * GeneTrail2 - An efficent library for interpreting genetic data
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

#ifndef GT2_ENRICHMENT_SET_LEVEL_STATISTICS_H
#define GT2_ENRICHMENT_SET_LEVEL_STATISTICS_H

#include "EnrichmentResult.h"

#include <genetrail2/core/HTest.h>
#include <genetrail2/core/Scores.h>
#include <genetrail2/core/GeneSetEnrichmentAnalysis.h>
#include <genetrail2/core/WeightedGeneSetEnrichmentAnalysis.h>
#include <genetrail2/core/OverRepresentationAnalysis.h>
#include <genetrail2/core/OneSampleTTest.h>

#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

namespace GeneTrail
{
	class Category;
	class DenseMatrix;
	class DenseMatrixSubset;

	using Indices = std::vector<size_t>;
	using IndexIterator = Indices::const_iterator;

	namespace StatTags
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
		struct SupportsIndices
		{
		};
		struct DoesNotSupportIndices
		{
		};
	}

	template <typename Mode = StatTags::Indirect,
	          typename Indices = StatTags::DoesNotSupportIndices,
	          typename Input = StatTags::Scores>
	class SetLevelStatistics
	{
		public:
		using RowWiseMode = Mode;
		using InputType = Input;
		using SupportsIndices = Indices;
	};

	class GT2_EXPORT StatisticsEnrichment : public SetLevelStatistics<>
	{
		public:
		using _viter = Scores::ConstScoreIterator;
		using Statistics = std::function<double(_viter, _viter)>;

		StatisticsEnrichment(const Statistics& test, const Scores& scores);
		StatisticsEnrichment(const Statistics& test, const Scores& scores,
		                     double cached);

		void setInputScores(const Scores& scores);

		bool canUseCategory(const Category&, size_t) const { return true; }

		std::tuple<double, double> computeScore(const Category& c) const;

		protected:
		virtual double cacheStatistic_() const;
		virtual double getExpectedValue_(const Category&) const;

		Statistics test_;
		Scores scores_;
		double cached_;
	};

	class GT2_EXPORT SumEnrichment final : public StatisticsEnrichment
	{
		public:
		SumEnrichment(const Scores& scores)
		    : StatisticsEnrichment(statistic::sum<double, _viter>, scores, 0.0)
		{
			cached_ = cacheStatistic_();
		}

		protected:
		double cacheStatistic_() const override;

		double getExpectedValue_(const Category& c) const override
		{
			return c.size() * cached_;
		}
	};

	class MaxMeanEnrichment final : public StatisticsEnrichment
	{
		public:
		MaxMeanEnrichment(const Scores& scores)
		    : StatisticsEnrichment(statistic::max_mean<double, _viter>, scores,
		                           0.0)
		{
		}

		protected:
		double cacheStatistic_() const override { return 0.0; }
	};

	class MedianEnrichment final : public StatisticsEnrichment
	{
		public:
		MedianEnrichment(const Scores& scores)
		    : StatisticsEnrichment(statistic::median<double, _viter>, scores)
		{
		}
	};

	class MeanEnrichment final : public StatisticsEnrichment
	{
		public:
		MeanEnrichment(const Scores& scores)
		    : StatisticsEnrichment(statistic::mean<double, _viter>, scores)
		{
		}
	};

	template <typename Test>
	class HTestEnrichmentBase : public SetLevelStatistics<StatTags::Direct>
	{
		public:
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
			return hits > 1;
		}

		double computeRowWisePValue(EnrichmentResult* result)
		{
			return HTestEnrichmentBase<Test>::computeRowWisePValue(test_,
			                                                       result);
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
			// TODO: update test here
		}

		bool canUseCategory(const Category&, size_t hits) const
		{
			return hits > 1;
		}

		double computeRowWisePValue(EnrichmentResult* result)
		{
			return HTestEnrichmentBase<OneSampleTTest<T>>::computeRowWisePValue(
			    test_, result);
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

	class Ora : public SetLevelStatistics<StatTags::Direct,
	                                      StatTags::DoesNotSupportIndices,
	                                      StatTags::Identifiers>
	{
		public:
		Ora(const Category& reference_set, const Category& test_set)
		    : test_(reference_set, test_set){};

		bool canUseCategory(const Category&, size_t) const { return true; }

		std::tuple<double, double> computeScore(const Category& c) const
		{
			//!! IMPORTANT:
			//!! Here we need to compare the number of hits and the expected number of hits.
			return std::make_tuple(test_.numberOfHits(c), test_.expectedNumberOfHits(c));
		}

		double computeRowWisePValue(EnrichmentResult* result) const
		{
			return test_.computePValue(*result->category);
		}

		private:
		OverRepresentationAnalysis test_;
	};

	class KolmogorovSmirnov : public SetLevelStatistics<StatTags::Direct, StatTags::SupportsIndices>
	{
		public:
		template <typename Iterator>
		KolmogorovSmirnov(const Iterator& beginIds, const Iterator& endIds, Order order)
		    : order_(order), ids_(beginIds, endIds)
		{
		}

		Order getOrder() const {
			return order_;
		}

		void setInputScores(Scores scores)
		{
			scores.sortByScore(order_);
			ids_.assign(scores.indices().begin(), scores.indices().end());
		}

		bool canUseCategory(const Category&, size_t) const { return true; }

		std::tuple<double, double> computeScore(const Category& category)
		{
			auto score =
			    test_.computeRunningSum(category, ids_.begin(), ids_.end());
			return std::make_tuple(score, 0.0);
		}

		std::tuple<double, double> computeScore(IndexIterator begin, IndexIterator end)
		{
			auto score =
			    test_.computeRunningSum(ids_.size(), begin, end);
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
		Order order_;
		std::vector<size_t> ids_;
		GeneSetEnrichmentAnalysis<big_float, int64_t> test_;
	};

	class WeightedKolmogorovSmirnov
	    : public SetLevelStatistics<StatTags::Indirect, StatTags::SupportsIndices>
	{
		public:
		WeightedKolmogorovSmirnov(const Scores& scores, Order order)
		    : scores_(scores.db()), order_(order)
		{
			setInputScores(scores);
		}

		Order getOrder() const {
			return order_;
		}

		void setInputScores(Scores scores)
		{
			scores.sortByScore(order_);
			test_.setScores(std::move(scores));
		}

		bool canUseCategory(const Category&, size_t) const { return true; }

		std::tuple<double, double> computeScore(const Category& category)
		{
			auto score = test_.computeRunningSum(category);
			return std::make_tuple(score, 0.0);
		}

		std::tuple<double, double> computeScore(IndexIterator begin,
		                                        IndexIterator end)
		{
			auto score = test_.computeRunningSum(begin, end);
			return std::make_tuple(score, 0.0);
		}

		private:
		Order order_;
		WeightedGeneSetEnrichmentAnalysis<double> test_;
	};
}

#endif // GT2_SET_LEVEL_STATISTICS_H
