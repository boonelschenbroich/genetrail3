#include "GeneSetFilters.h"

#include "GeneSet.h"

#include <stdexcept>

namespace GeneTrail
{
	namespace GeneSetFilter
	{
		UpperQuantileFilter::UpperQuantileFilter(double quantile)
		    : quantile_(quantile), threshold_(0.0)
		{
			if(quantile_ < 0.0 || quantile_ > 1.0) {
				throw std::out_of_range("Quantile must lie in [0, 1]");
			}
		}

		void UpperQuantileFilter::setup(const GeneSet& gene_set)
		{
			if(gene_set.empty()) {
				return;
			}

			auto sorted = gene_set.getDecreasinglySortedScores();

			int i = std::floor(sorted.size() * quantile_);
			threshold_ = sorted[i].second;
		}

		LowerQuantileFilter::LowerQuantileFilter(double quantile)
		    : quantile_(quantile), threshold_(0.0)
		{
			if(quantile_ < 0.0 || quantile_ > 1.0) {
				throw std::out_of_range("Quantile must lie in [0, 1]");
			}
		}

		void LowerQuantileFilter::setup(const GeneSet& gene_set)
		{
			if(gene_set.empty()) {
				return;
			}

			auto sorted = gene_set.getIncreasinglySortedScores();

			int i = std::floor(sorted.size() * quantile_);
			threshold_ = sorted[i].second;
		}

		QuantileFilter::QuantileFilter(double quantile)
		    : quantile_(quantile), lower_threshold_(0.0), upper_threshold_(0.0)
		{
			if(quantile_ < 0.0 || quantile_ > 1.0) {
				throw std::out_of_range("Quantile must lie in [0, 1]");
			}
		}

		void QuantileFilter::setup(const GeneSet& gene_set)
		{
			if(gene_set.empty()) {
				return;
			}

			auto sorted = gene_set.getIncreasinglySortedScores();

			int i = std::floor(sorted.size() * quantile_ * 0.5);
			lower_threshold_ = sorted[i].second;
			upper_threshold_ = sorted[sorted.size() - i - 1].second;
		}
	}
}
