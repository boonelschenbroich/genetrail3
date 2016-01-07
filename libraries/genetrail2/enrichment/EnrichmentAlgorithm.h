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
#ifndef GT2_ENRICHMENT_ENRICHMENT_ALGORITHM_H
#define GT2_ENRICHMENT_ENRICHMENT_ALGORITHM_H

#include "SetLevelStatistics.h"

#include <genetrail2/core/compat.h>
#include <genetrail2/core/Exception.h>
#include <genetrail2/core/macros.h>

#include <memory>

namespace GeneTrail
{
	enum class PValueMode { RowWise, ColumnWise, Restandardize };

	class GT2_EXPORT EnrichmentAlgorithm
	{
		public:
		EnrichmentAlgorithm(PValueMode mode);

		PValueMode pValueMode() const { return mode_; }
		bool pValuesComputed() const;

		virtual bool canUseCategory(const Category& c, size_t hits) const = 0;
		virtual bool rowWisePValueIsDirect() const = 0;
		virtual bool supportsIndices() const = 0;
		virtual Order getOrder() const = 0;

		virtual void setScores(const Scores& scores) = 0;

		virtual std::unique_ptr<EnrichmentResult>
		computeEnrichment(const std::shared_ptr<Category>& c) = 0;

		virtual std::tuple<double, double>
		computeEnrichmentScore(const Category& c) = 0;

		virtual std::tuple<double, double>
		computeEnrichmentScore(IndexIterator begin, IndexIterator end) = 0;

		private:
		PValueMode mode_;
	};

	using EnrichmentAlgorithmPtr = std::unique_ptr<EnrichmentAlgorithm>;

	namespace internal
	{
		template <typename Statistics>
		class EnrichmentWorker final : public EnrichmentAlgorithm
		{
			public:
			template <typename... Ts>
			EnrichmentWorker(PValueMode mode, Ts&&... ts)
			    : EnrichmentAlgorithm(mode),
			      statistics_(std::forward<Ts>(ts)...)
			{
			}

			void setScores(const Scores& scores) override
			{
				setScoresDispatch_(scores, typename Statistics::InputType());
			}

			bool canUseCategory(const Category& c, size_t hits) const override
			{
				return statistics_.canUseCategory(c, hits);
			}

			std::tuple<double, double>
			computeEnrichmentScore(const Category& c) override
			{
				return statistics_.computeScore(c);
			}

			std::tuple<double, double>
			computeEnrichmentScore(IndexIterator begin,
			                       IndexIterator end) override
			{
				return computeEnrichmentScoreDispatch_(
				    typename Statistics::SupportsIndices(), begin, end);
			}

			std::unique_ptr<EnrichmentResult>
			computeEnrichment(const std::shared_ptr<Category>& c) override
			{
				auto result = std::make_unique<EnrichmentResult>(c);

				std::tie(result->score, result->expected_score) =
				    computeEnrichmentScore(*c);
				result->enriched = result->score > result->expected_score;

				computePValueDispatch_(result.get(),
				                       typename Statistics::RowWiseMode());

				return result;
			}

			bool rowWisePValueIsDirect() const override
			{
				return rowWisePValueIsDirectDispatch_(
				    typename Statistics::RowWiseMode());
			}

			bool supportsIndices() const override
			{
				return supportsIndicesDispatch_(
				    typename Statistics::SupportsIndices());
			}

			Order getOrder() const override
			{
				return getOrderDispatch_(typename Statistics::SupportsIndices());
			}

			private:
			void setScoresDispatch_(const Scores& scores, StatTags::Scores)
			{
				statistics_.setInputScores(scores);
			}

			void setScoresDispatch_(const Scores&, StatTags::Identifiers) {}

			void computePValueDispatch_(EnrichmentResult* result,
			                            StatTags::Direct)
			{
				if(pValuesComputed()) {
					result->pvalue = statistics_.computeRowWisePValue(result);
				}
			}

			void computePValueDispatch_(EnrichmentResult*, StatTags::Indirect)
			{
			}

			bool rowWisePValueIsDirectDispatch_(StatTags::Direct) const
			{
				return true;
			}

			bool rowWisePValueIsDirectDispatch_(StatTags::Indirect) const
			{
				return false;
			}

			bool supportsIndicesDispatch_(StatTags::SupportsIndices) const
			{
				return true;
			}

			bool supportsIndicesDispatch_(StatTags::DoesNotSupportIndices) const
			{
				return false;
			}

			Order getOrderDispatch_(StatTags::SupportsIndices) const
			{
				return statistics_.getOrder();
			}

			Order getOrderDispatch_(StatTags::DoesNotSupportIndices) const
			{
				throw NotImplemented(__FILE__, __LINE__,
				                    "This type does not support the concept of "
				                    "sort orders. If you encounter this message "
				                    "this means there is a bug.");
			}

			std::tuple<double, double>
			computeEnrichmentScoreDispatch_(StatTags::SupportsIndices,
			                                IndexIterator begin,
			                                IndexIterator end)
			{
				return statistics_.computeScore(begin, end);
			}

			std::tuple<double, double>
			    computeEnrichmentScoreDispatch_(StatTags::DoesNotSupportIndices,
			                                    IndexIterator, IndexIterator)
			{
				throw NotImplemented(__FILE__, __LINE__,
				                     "This type does not implement the "
				                     "computation of enrichment scores using "
				                     "indices.");
			}

			Statistics statistics_;
		};
	}

	template <typename Statistics, typename... Ts>
	EnrichmentAlgorithmPtr GT2_EXPORT
	    createEnrichmentAlgorithm(PValueMode mode, Ts&&... ts)
	{
		return std::make_unique<internal::EnrichmentWorker<Statistics>>(
		    mode, std::forward<Ts>(ts)...);
	}
}

#endif // GT2_ENRICHMENT_ENRICHMENT_ALGORITHM_H
