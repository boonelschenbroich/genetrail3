#ifndef GT2_CORE_PVALUE_H
#define GT2_CORE_PVALUE_H

#include <algorithm>
#include <cmath>
#include <map>
#include <functional>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>

namespace GeneTrail
{
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
			for(int i = 1; i < adjustedPValues.size(); ++i) {
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
			return stepDown(adjustByOrder(pvalues, fdr_func));
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
			for(int i = 0; i < adj.size(); ++i) {
				q += 1 / ((float_type)i + 1.0);
			}
			for(int i = 0; i < adj.size(); ++i) {
				adj[i].second *= q;
			}
			return stepDown(adj);
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		// Step-up adjustments
		////////////////////////////////////////////////////////////////////////////////////////////////////

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

		static std::vector<pValue<float_type>>
		adjustPValues(const std::vector<pValue<float_type>>& pvalues, std::string method)
		{
			std::map<std::string, std::function<std::vector<pValue<float_type>>(
			                          const std::vector<pValue<float_type>>&)>> methods;
			methods["bonferroni"] = bonferroni;
			methods["sidak"] = sidak;
			methods["holm"] = holm;
			methods["holm_sidak"] = holm_sidak;
			methods["finner"] = finner;
			methods["benjamini_hochberg"] = benjamini_hochberg;
			methods["benjamini_yekutieli"] = benjamini_yekutieli;
			methods["hochberg"] = hochberg;
			methods["simes"] = simes;
			auto it = methods.find(method);
			if(it != methods.end()) {
				return methods[method](pvalues);
			}
			return pvalues;
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
			for(int i = 0; i < pvals.size(); ++i) {
				x += std::log(pvals[i].second);
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
			for(int i = 0; i < pvalues.size(); ++i) {
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
