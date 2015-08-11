#ifndef GT2_ENRICHMENT_ALGORITHM_H
#define GT2_ENRICHMENT_ALGORITHM_H

#include "compat.h"
#include "SetLevelStatistics.h"
#include "macros.h"

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

		virtual void setScores(const Scores& scores) = 0;
		virtual bool canUseCategory(const Category& c, size_t hits) const = 0;
		virtual std::unique_ptr<EnrichmentResult>
		computeEnrichment(const std::shared_ptr<Category>& c) = 0;
		virtual double computeEnrichmentScore(const Category& c) = 0;
		virtual bool rowWisePValueIsDirect() const = 0;

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

			double computeEnrichmentScore(const Category& c) override
			{
				return statistics_.computeScore(c);
			}

			std::unique_ptr<EnrichmentResult>
			computeEnrichment(const std::shared_ptr<Category>& c) override
			{
				auto result = std::make_unique<EnrichmentResult>(c);

				result->score = computeEnrichmentScore(*c);

				computePValueDispatch_(result.get(),
				                       typename Statistics::RowWiseMode());

				return result;
			}

			bool rowWisePValueIsDirect() const override
			{
				return rowWisePValueIsDirectDispatch_(
				    typename Statistics::RowWiseMode());
			}

			private:
			void setScoresDispatch_(const Scores& scores, SetLevelStatistics::Scores) {
				statistics_.setInputScores(scores);
			}

			void setScoresDispatch_(const Scores&, SetLevelStatistics::Identifiers) {
			}

			void computePValueDispatch_(EnrichmentResult* result,
			                            SetLevelStatistics::Direct)
			{
				if(pValuesComputed()) {
					result->pvalue = statistics_.computeRowWisePValue(result);
				}
			}

			void computePValueDispatch_(EnrichmentResult*,
			                            SetLevelStatistics::Indirect)
			{
			}

			bool
			    rowWisePValueIsDirectDispatch_(SetLevelStatistics::Direct) const
			{
				return true;
			}

			bool rowWisePValueIsDirectDispatch_(
			    SetLevelStatistics::Indirect) const
			{
				return false;
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

#endif // GT2_ENRICHMENT_ALGORITHM_H
