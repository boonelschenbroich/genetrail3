/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2014-2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *               2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_PVALUE_H
#define GT2_CORE_PVALUE_H

#include "Exception.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <functional>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/optional/optional.hpp>

namespace GeneTrail
{
	enum class MultipleTestingCorrection {
		Bonferroni,
		Sidak,
		Holm,
		HolmSidak,
		Finner,
		BenjaminiHochberg,
		BenjaminiYekutieli,
		Hochberg,
		Simes,
		GSEA
	};

	template <typename float_type> using pValue = typename std::pair<std::string, float_type>;
	template <typename float_type = double>
	class pvalue
	{
		public:
		/**
		 * Cummulative min function.
		 *
		 * @param pvalues P-values
		 * @return Adjusted p-values
		 */
		static std::vector<pValue<float_type>>
		stepUp(const std::vector<pValue<float_type>>& pvalues)
		{
			std::vector<pValue<float_type>> adjustedPValues(pvalues);

			if(adjustedPValues.empty()) {
				return adjustedPValues;
			}

			adjustedPValues.back().second = std::min(adjustedPValues.back().second, 1.0);
			for(int i = adjustedPValues.size() - 2; i >= 0; --i) {
				adjustedPValues[i].second =
				    std::min(std::min(adjustedPValues[i + 1].second,
				                      adjustedPValues[i].second),
				             1.0);
			}
			return adjustedPValues;
		}

		/**
		 * Cummulative max function.
		 *
		 * @param pvalues P-values
		 * @return Adjusted p-values
		 */
		static std::vector<pValue<float_type>>
		stepDown(const std::vector<pValue<float_type>>& pvalues)
		{
			std::vector<pValue<float_type>> adjustedPValues(pvalues);
			adjustedPValues[0].second = std::min(adjustedPValues[0].second, 1.0);
			for(size_t i = 1; i < adjustedPValues.size(); ++i) {
				adjustedPValues[i].second =
				    std::min(std::max(adjustedPValues[i - 1].second,
				                      adjustedPValues[i].second),
				             1.0);
			}
			return adjustedPValues;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		// One step adjustments
		////////////////////////////////////////////////////////////////////////////////////////////////////

		template <typename func>
		static std::vector<pValue<float_type>> adjustBySize(const std::vector<pValue<float_type>>& pvalues,
		                                   func f)
		{
			std::vector<pValue<float_type>> adjustedPValues(pvalues);
			size_t n = pvalues.size();
			for(unsigned int i = 0; i < adjustedPValues.size(); ++i) {
				adjustedPValues[i].second =
				    std::min(f(adjustedPValues[i], n), 1.0);
			}
			return adjustedPValues;
		}

		static float_type bonferroni_func(pValue<float_type> p, int n)
		{
			return p.second * n;
		}

		/**
		 * Bonferroni p-value adjustment.
		 *
		 * Reference:
		 *
		 * @param pvalues P-values
		 * @return Adjusted p-values
		 */
		static std::vector<pValue<float_type>>
		bonferroni(const std::vector<pValue<float_type>>& pvalues)
		{
			return adjustBySize(pvalues, bonferroni_func);
		}

		static float_type sidak_func(pValue<float_type> p, int n)
		{
			return 1.0 - std::pow((1.0 - p.second), n);
		}

		/**
		 * Sidak p-value adjustment.
		 *
		 * Reference: Zbyněk Šidák (1967). Rectangular confidence regions for
		 *the means of multivariate normal distributions
		 *
		 * @param pvalues P-values
		 * @return Adjusted p-values
		 */
		static std::vector<pValue<float_type>> sidak(const std::vector<pValue<float_type>>& pvalues)
		{
			return adjustBySize(pvalues, sidak_func);
		}

		template <typename func>
		static std::vector<pValue<float_type>>
		adjustByOrder(const std::vector<pValue<float_type>>& pvalues, func f)
		{
			std::vector<pValue<float_type>> adjustedPValues(pvalues);
			std::sort(
			    adjustedPValues.begin(), adjustedPValues.end(),
			    [](const pValue<float_type>& a, const pValue<float_type>& b) {
				    return a.second < b.second;
				});
			size_t n = adjustedPValues.size();
			for(unsigned int i = 0; i < n; ++i) {
				adjustedPValues[i].second = f(adjustedPValues[i], n, i + 1);
			}
			return adjustedPValues;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		// Step-down adjustments
		////////////////////////////////////////////////////////////////////////////////////////////////////

		static float_type hochberg_func(pValue<float_type> p, int n, int i)
		{
			return p.second * (n - i + 1);
		}

		/**
		 * Holm p-value adjustment.
		 *
		 * Reference: Holm, S. (1979).  A simple sequentially rejective multiple
		 *test procedure.
		 *
		 * @param pvalues P-values
		 * @return Adjusted p-values
		 */
		static std::vector<pValue<float_type>> holm(const std::vector<pValue<float_type>>& pvalues)
		{
			return stepDown(adjustByOrder(pvalues, hochberg_func));
		}

		static float_type holm_sidak_func(pValue<float_type> p, int n, int i)
		{
			return 1.0 - std::pow((1.0 - p.second), n - i + 1);
		}

		/**
		 * Holm-Sidak p-value adjustment.
		 *
		 * Reference:
		 *
		 * @param pvalues P-values
		 * @return Adjusted p-values
		 */
		static std::vector<pValue<float_type>>
		holm_sidak(const std::vector<pValue<float_type>>& pvalues)
		{
			return stepDown(
			    adjustByOrder(pvalues, holm_sidak_func));
		}

		static float_type finner_func(pValue<float_type> p, int n, int i)
		{
			return 1.0 - std::pow((1.0 - p.second),
			                      ((float_type)n) / ((float_type)i));
		}

		/**
		 *  Finner p-value adjustment.
		 *
		 *  Reference:
		 *
		 *  @param pvalues P-values
		 *  @return Adjusted p-values
		 */
		static std::vector<pValue<float_type>>
		finner(const std::vector<pValue<float_type>>& pvalues)
		{
			return stepDown(adjustByOrder(pvalues, finner_func));
		}

		static float_type fdr_func(pValue<float_type> p, int n, int i)
		{
			return p.second * ((float_type)n) / ((float_type)i);
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		// Step-up adjustments
		////////////////////////////////////////////////////////////////////////////////////////////////////


		/**
		 * Benjamini-Hochberg p-value adjustment.
		 *
		 * Reference: 	Benjamini, Y., and Hochberg, Y. (1995).  Controlling the
		 *false discovery rate: a practical and powerful approach to multiple
		 *testing.
		 *
		 * @param pvalues P-values
		 * @return Adjusted p-values
		 */
		static std::vector<pValue<float_type>>
		benjamini_hochberg(const std::vector<pValue<float_type>>& pvalues)
		{
			return stepUp(adjustByOrder(pvalues, fdr_func));
		}

		/**
		 * Benjamini-Yekutieli p-value adjustment.
		 *
		 * Reference:	Benjamini, Y., and Yekutieli, D. (2001).  The control of
		 *the false discovery rate in multiple testing under dependency.
		 *
		 * @param pvalues P-values
		 * @return Adjusted p-values
		 */
		static std::vector<pValue<float_type>>
		benjamini_yekutieli(const std::vector<pValue<float_type>>& pvalues)
		{
			std::vector<pValue<float_type>> adj =
			    adjustByOrder(pvalues, fdr_func);
			float_type q = 0.0;
			for(size_t i = 0; i < adj.size(); ++i) {
				q += 1 / ((float_type)i + 1.0);
			}
			for(size_t i = 0; i < adj.size(); ++i) {
				adj[i].second *= q;
			}
			return stepUp(adj);
		}

		/**
		 * Hochberg p-value adjustment.
		 *
		 * Reference: Hochberg, Y. (1988).  A sharper Bonferroni procedure for
		 *multiple tests of significance.
		 *
		 * @param pvalues P-values
		 * @return Adjusted p-values
		 */
		static std::vector<pValue<float_type>>
		hochberg(const std::vector<pValue<float_type>>& pvalues)
		{
			return stepUp(adjustByOrder(pvalues, hochberg_func));
		}

		/**
		 * Simes p-value adjustment.
		 *
		 * Reference: Simes RJ (1986) "An improved Bonferroni procedure for
		 *multiple tests of significance"
		 *
		 * @param pvalues P-values
		 * @return Adjusted p-values
		 */
		static std::vector<pValue<float_type>> simes(const std::vector<pValue<float_type>>& pvalues)
		{
			return stepUp(adjustByOrder(pvalues, fdr_func));
		}

		static std::vector<pValue<float_type>> noop(const std::vector<pValue<float_type>>& pvalues)
		{
			return pvalues;
		}

		using CorrectionMethod = std::vector<pValue<float_type>>(*)(const std::vector<pValue<float_type>>&);

		static boost::optional<MultipleTestingCorrection> getCorrectionMethod(const std::string& method) {
			if(method == "bonferroni") {
				return boost::make_optional(MultipleTestingCorrection::Bonferroni);
			}

			if(method == "sidak") {
				return boost::make_optional(MultipleTestingCorrection::Sidak);
			}

			if(method == "holm") {
				return boost::make_optional(MultipleTestingCorrection::Holm);
			}

			if(method == "holm_sidak" || method == "holm-sidak") {
				return boost::make_optional(MultipleTestingCorrection::HolmSidak);
			}

			if(method == "finner") {
				return boost::make_optional(MultipleTestingCorrection::Finner);
			}

			if(method == "benjamini_hochberg" || method == "benjamini-hochberg") {
				return boost::make_optional(MultipleTestingCorrection::BenjaminiHochberg);
			}

			if(method == "benjamini_yekutieli" || method == "benjamini-yekutieli") {
				return boost::make_optional(MultipleTestingCorrection::BenjaminiYekutieli);
			}

			if(method == "hochberg") {
				return boost::make_optional(MultipleTestingCorrection::Hochberg);
			}

			if(method == "simes") {
				return boost::make_optional(MultipleTestingCorrection::Simes);
			}

			if(method == "gsea" || method == "GSEA") {
				return boost::make_optional(MultipleTestingCorrection::GSEA);
			}

			return boost::none;
		}

		static CorrectionMethod getCorrectionMethod(MultipleTestingCorrection method) {
			switch(method) {
				case MultipleTestingCorrection::Bonferroni: return bonferroni;
				case MultipleTestingCorrection::Sidak: return sidak;
				case MultipleTestingCorrection::Holm: return holm;
				case MultipleTestingCorrection::HolmSidak: return holm_sidak;
				case MultipleTestingCorrection::Finner: return finner;
				case MultipleTestingCorrection::BenjaminiHochberg: return benjamini_hochberg;
				case MultipleTestingCorrection::BenjaminiYekutieli: return benjamini_yekutieli;
				case MultipleTestingCorrection::Hochberg: return hochberg;
				case MultipleTestingCorrection::Simes: return simes;
				case MultipleTestingCorrection::GSEA: return noop;
				default:
					throw NotImplemented(__FILE__, __LINE__, "The requested correction method has not yet been implemented.");
			}
		}

		static std::vector<pValue<float_type>>
		adjustPValues(const std::vector<pValue<float_type>>& pvalues, MultipleTestingCorrection method)
		{
			return getCorrectionMethod(method)(pvalues);
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		// P-value aggregation
		////////////////////////////////////////////////////////////////////////////////////////////////////

		/**
		 * Fisher method to aggregate p-values.
		 *
		 * @param pvalues P-values
		 * @return Aggregated p-value
		 */
		static float_type fisher(const std::vector<pValue<float_type>>& pvalues)
		{
			std::vector<pValue<float_type>> pvals(pvalues);
			float_type x = 0.0;
			for(const auto& pval : pvals) {
				x += std::log(pval.second);
			}
			boost::math::chi_squared dist(2 * pvalues.size());
			return boost::math::cdf(boost::math::complement(dist, -2.0 * x));
		}

		/**
		 * Stouffer method to aggregate p-values.
		 *
		 * @param pvalues P-values
		 * @return Aggregated p-value
		 */
		static float_type stouffer(const std::vector<pValue<float_type>>& pvalues,
		                           const std::vector<float_type>& weights)
		{
			std::vector<pValue<float_type>> pvals(pvalues);
			boost::math::normal dist(0.0, 1.0);
			float_type zi = 0.0;
			float_type ws = 0.0;
			for(size_t i = 0; i < pvalues.size(); ++i) {
				ws += weights[i] * weights[i];
				zi +=
				    quantile(complement(dist, pvalues[i].second)) * weights[i];
			}
			float_type z = zi / std::sqrt(ws);
			return boost::math::cdf(boost::math::complement(dist, z));
		}
	};
}

#endif // GT2_CORE_PVALUE_H
