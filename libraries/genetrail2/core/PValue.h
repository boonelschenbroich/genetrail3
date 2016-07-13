/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2014-2016 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#include "macros.h"

#include <algorithm>
#include <cmath>
#include <utility>
#include <tuple>
#include <type_traits>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/optional/optional.hpp>

namespace GeneTrail
{

/**
 * An enum for specifying a method for multiple testing correction.
 */
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

/**
 * This namespace contains methods that implement algorithms for multiple
 * testing correction and p-value aggregation.
 *
 * All algorithms are implemented as template functions that require two
 * arguments. The first argument is a container (such as std::vector)
 * containing p-values. These p-values can be stored in arbitrarly complex
 * objects. Due to this, the second argument is a callable that extracts a
 * (mutable) reference to the p-value. The adjusted p-values are returned by
 * the methods. Due to the use of move semantics no uncessessary copies are
 * being made.
 */
namespace pvalue
{

/**
 * This helper template is needed to ensure, that we return and create
 * copies of the object containing the p-values.
 */
template <typename T>
using rcref = std::remove_const_t<std::remove_reference_t<T>>;

/**
 * A helper functor that projects the second value out of a tuple.
 * We use this throughout GeneTrail2 to extract the p-values from
 * a std::pair<string, double>.
 */
struct get_second
{
	template <typename Tuple> double operator()(const Tuple& t) const
	{
		return std::get<1>(t);
	}

	// Important: the method needs to return a mutable reference!
	template <typename Tuple> double& operator()(Tuple& t) const
	{
		return std::get<1>(t);
	}
};

/**
 * Cumulative min function. This is used as a building block for other
 * adjustment methods.
 *
 * @param pvalues A container of p-values
 * @param pvalue A callable that extracts a reference to the p-value
 *               from an object in the container.
 *
 * @return Adjusted p-values
 */
template <typename PValues, typename Access>
rcref<PValues> stepUp(PValues&& pvalues, Access&& pvalue)
{
	if(pvalues.empty()) {
		return pvalues;
	}

	rcref<PValues> adjusted(std::forward<PValues>(pvalues));

	pvalue(adjusted.back()) = std::min(pvalue(adjusted.back()), 1.0);
	for(int64_t i = adjusted.size() - 2; i >= 0; --i) {
		pvalue(adjusted[i]) = std::min(
		    std::min(pvalue(adjusted[i + 1]), pvalue(adjusted[i])), 1.0);
	}

	return adjusted;
}

/**
 * Cumulative max function. This is used as a building block for other
 * adjustment methods.
 *
 * @param pvalues A container of p-values
 * @param pvalue A callable that extracts a reference to the p-value
 *               from an object in the container.
 *
 * @return Adjusted p-values
 */
template <typename PValues, typename Access>
rcref<PValues> stepDown(PValues&& pvalues, Access&& pvalue)
{
	rcref<PValues> adjusted(std::forward<PValues>(pvalues));

	pvalue(adjusted[0]) = std::min(pvalue(adjusted[0]), 1.0);

	for(size_t i = 1; i < adjusted.size(); ++i) {
		pvalue(adjusted[i]) = std::min(
		    std::max(pvalue(adjusted[i - 1]), pvalue(adjusted[i])), 1.0);
	}

	return adjusted;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// One step adjustments
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename PValues, typename Access, typename func>
rcref<PValues> adjustBySize(PValues&& pvalues, Access&& pvalue, func f)
{
	rcref<PValues> adjustedPValues(std::forward<PValues>(pvalues));

	size_t n = pvalues.size();
	for(auto& p : adjustedPValues) {
		pvalue(p) = std::min(f(pvalue(p), n), 1.0);
	}

	return adjustedPValues;
}

inline double bonferroni_func(double p, size_t n) { return p * n; }

/**
 * Bonferroni p-value adjustment.
 *
 * Reference:
 *
 * @TODO citation
 * @param pvalues A container of p-values
 * @param pvalue A callable that extracts a reference to the p-value from an
 *               object in the container.
 *
 * @return Adjusted p-values
 */
template <typename PValues, typename Access>
rcref<PValues> bonferroni(PValues&& pvalues, Access&& access)
{
	return adjustBySize(std::forward<PValues>(pvalues),
	                    std::forward<Access>(access), bonferroni_func);
}

inline double sidak_func(double p, size_t n)
{
	return 1.0 - std::pow((1.0 - p), n);
}

/**
 * Sidak p-value adjustment.
 *
 * Reference: Zbyněk Šidák (1967). Rectangular confidence regions for
 * the means of multivariate normal distributions
 *
 * @TODO citation
 * @param pvalues A container of p-values
 * @param pvalue A callable that extracts a reference to the p-value
 *               from an object in the container.
 *
 * @return Adjusted p-values
 */
template <typename PValues, typename Access>
rcref<PValues> sidak(PValues&& pvalues, Access&& access)
{
	return adjustBySize(std::forward<PValues>(pvalues),
	                    std::forward<Access>(access), sidak_func);
}

template <typename PValues, typename Access, typename Func>
rcref<PValues> adjustByOrder(PValues&& pvalues, Access&& pvalue, Func&& f)
{
	rcref<PValues> adjusted(std::forward<PValues>(pvalues));

	auto pval_less =
	    [pvalue](auto&& a, auto&& b) { return pvalue(a) < pvalue(b); };

	std::sort(adjusted.begin(), adjusted.end(), pval_less);

	size_t n = adjusted.size();
	for(unsigned int i = 0; i < n; ++i) {
		pvalue(adjusted[i]) = f(pvalue(adjusted[i]), n, i + 1);
	}

	return adjusted;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Step-down adjustments
////////////////////////////////////////////////////////////////////////////////////////////////////

inline double hochberg_func(double p, int n, int i) { return p * (n - i + 1); }

/**
 * Holm p-value adjustment.
 *
 * Reference: Holm, S. (1979).  A simple sequentially rejective multiple
 * test procedure.
 *
 * @TODO citation
 * @param pvalues A container of p-values
 * @param pvalue A callable that extracts a reference to the p-value
 *               from an object in the container.
 *
 * @return Adjusted p-values
 */
template <typename PValues, typename Access>
rcref<PValues> holm(PValues&& pvalues, Access&& access)
{
	return stepDown(
	    adjustByOrder(std::forward<PValues>(pvalues), access, hochberg_func),
	    access);
}

inline double holm_sidak_func(double p, int n, int i)
{
	return 1.0 - std::pow(1.0 - p, n - i + 1);
}

/**
 * Holm-Sidak p-value adjustment.
 *
 * Reference:
 *
 * @TODO citation
 * @param pvalues A container of p-values
 * @param pvalue A callable that extracts a reference to the p-value
 *               from an object in the container.
 *
 * @return Adjusted p-values
 */
template <typename PValues, typename Access>
rcref<PValues> holm_sidak(PValues&& pvalues, Access&& access)
{
	return stepDown(
	    adjustByOrder(std::forward<PValues>(pvalues), access, holm_sidak_func),
	    std::forward<Access>(access));
}

inline double finner_func(double p, int n, int i)
{
	return 1.0 - std::pow(1.0 - p, static_cast<double>(n) / i);
}

/**
 *  Finner p-value adjustment.
 *
 *  Reference:
 *
 * @TODO citation
 * @param pvalues A container of p-values
 * @param pvalue A callable that extracts a reference to the p-value
 *               from an object in the container.
 *
 * @return Adjusted p-values
 */
template <typename PValues, typename Access>
rcref<PValues> finner(PValues&& pvalues, Access&& access)
{
	return stepDown(
	    adjustByOrder(std::forward<PValues>(pvalues), access, finner_func),
	    std::forward<Access>(access));
}

inline double fdr_func(double p, int n, int i)
{
	// The parenthesis are important due to double <-> int conversions
	return (n * p) / i;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Step-up adjustments
////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Benjamini-Hochberg p-value adjustment.
 *
 * Reference: Benjamini, Y., and Hochberg, Y. (1995).  Controlling the
 * false discovery rate: a practical and powerful approach to multiple
 * testing.
 *
 * Benjamini, Y., and Yekutieli, D. (1999). Resampling-based false
 * discovery rate controlling multiple test procedures for correlated
 * test statistics
 *
 * @param pvalues P-values
 * @return Adjusted p-values
 */
template <typename PValues, typename Access>
rcref<PValues> benjamini_hochberg(PValues&& pvalues, Access&& access)
{
	return stepUp(
	    adjustByOrder(std::forward<PValues>(pvalues), access, fdr_func),
	    std::forward<Access>(access));
}

/**
 * Benjamini-Yekutieli p-value adjustment.
 *
 * Reference: Benjamini, Y., and Yekutieli, D. (2001).  The control of
 * the false discovery rate in multiple testing under dependency.
 *
 * @param pvalues P-values
 * @return Adjusted p-values
 */
template <typename PValues, typename Access>
rcref<PValues> benjamini_yekutieli(PValues&& pvalues, Access&& pvalue)
{
	auto adj = adjustByOrder(std::forward<PValues>(pvalues), pvalue, fdr_func);

	double q = 0.0;
	for(size_t i = 0; i < adj.size(); ++i) {
		q += 1 / (i + 1.0);
	}

	for(auto& p : adj) {
		pvalue(p) *= q;
	}

	return stepUp(std::move(adj), std::forward<Access>(pvalue));
}

/**
 * Hochberg p-value adjustment.
 *
 * Reference: Hochberg, Y. (1988).  A sharper Bonferroni procedure for
 * multiple tests of significance.
 *
 * @param pvalues P-values
 * @return Adjusted p-values
 */
template <typename PValues, typename Access>
rcref<PValues> hochberg(PValues&& pvalues, Access&& access)
{
	return stepUp(
	    adjustByOrder(std::forward<PValues>(pvalues), access, hochberg_func),
	    std::forward<Access>(access));
}

/**
 * Simes p-value adjustment.
 *
 * Reference: Simes RJ (1986) "An improved Bonferroni procedure for
 * multiple tests of significance"
 *
 * @param pvalues P-values
 * @return Adjusted p-values
 */
template <typename PValues, typename Access>
rcref<PValues> simes(PValues&& pvalues, Access&& access)
{
	return stepUp(
	    adjustByOrder(std::forward<PValues>(pvalues), access, fdr_func),
	    std::forward<Access>(access));
}

/**
 * Takes an input string and returns a matching multiple testing correction
 * method. If the requested method is unknown an empty option is returned.
 *
 * @param method A string specifying an input method.
 *
 * @return An optional containing the enum belonging to the requested
 * correction method. boost::none if the string could not be matched.
 */
GT2_EXPORT
boost::optional<MultipleTestingCorrection>
getCorrectionMethod(const std::string& method);

/**
 * Adjust a vector of p-values.
 *
 * @param pvalues The vector of pvalues that should be adjusted.
 * @param method The method for multiple testing correction.
 *
 * @return A copy of the input vector containing the adjusted p-values.
 */
template <typename PValues, typename Access>
rcref<PValues> adjustPValues(PValues&& pvalues, Access&& access, MultipleTestingCorrection method)
{
	switch(method) {
		case MultipleTestingCorrection::Bonferroni:
			return bonferroni(std::forward<PValues>(pvalues), std::forward<Access>(access));
		case MultipleTestingCorrection::Sidak:
			return sidak(std::forward<PValues>(pvalues), std::forward<Access>(access));
		case MultipleTestingCorrection::Holm:
			return holm(std::forward<PValues>(pvalues), std::forward<Access>(access));
		case MultipleTestingCorrection::HolmSidak:
			return holm_sidak(std::forward<PValues>(pvalues), std::forward<Access>(access));
		case MultipleTestingCorrection::Finner:
			return finner(std::forward<PValues>(pvalues),std::forward<Access>(access));
		case MultipleTestingCorrection::BenjaminiHochberg:
			return benjamini_hochberg(std::forward<PValues>(pvalues), std::forward<Access>(access));
		case MultipleTestingCorrection::BenjaminiYekutieli:
			return benjamini_yekutieli(std::forward<PValues>(pvalues), std::forward<Access>(access));
		case MultipleTestingCorrection::Hochberg:
			return hochberg(std::forward<PValues>(pvalues), std::forward<Access>(access));
		case MultipleTestingCorrection::Simes:
			return simes(std::forward<PValues>(pvalues), std::forward<Access>(access));
		case MultipleTestingCorrection::GSEA:
			return pvalues;
		default:
			throw NotImplemented(__FILE__, __LINE__, "The requested correction "
			                                         "method has not yet been "
			                                         "implemented.");
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// P-value aggregation
////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Fisher method to aggregate p-values.
 * @TODO citation
 *
 * @param pvalues P-values
 * @return Aggregated p-value
 */
template <typename PValues, typename Access>
double fisher(const PValues& pvalues, Access&& pvalue)
{
	double x = 0.0;
	for(const auto& p : pvalues) {
		x += std::log(pvalue(p));
	}

	boost::math::chi_squared dist(2 * pvalues.size());
	return boost::math::cdf(boost::math::complement(dist, -2.0 * x));
}

/**
 * Stouffer method to aggregate p-values.
 *
 * @todo citation
 *
 * @param pvalues The vector of pvalues that should be aggregated.
 * @param weights A container of weights that determine the influence of each
 *                p-value. If unsure use uniform weights.
 * @param pvalue A callable that extracts the pvalue from each object
 *               stored in pvalues. The returned pvalue does not need to be
 *               a mutable reference.
 *
 * @return Aggregated p-value
 */
template <typename PValues, typename Weights, typename Access>
double stouffer(const PValues& pvalues, const Weights& weights, Access&& pvalue)
{
	boost::math::normal dist(0.0, 1.0);
	double zi = 0.0;
	double ws = 0.0;

	for(size_t i = 0; i < pvalues.size(); ++i) {
		ws += weights[i] * weights[i];
		zi += quantile(complement(dist, pvalue(pvalues[i]))) * weights[i];
	}

	double z = zi / std::sqrt(ws);
	return boost::math::cdf(boost::math::complement(dist, z));
}

} // namespace pvalue
} // namespace GeneTrail

#endif // GT2_CORE_PVALUE_H
