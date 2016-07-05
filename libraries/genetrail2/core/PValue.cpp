/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2016 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#include "PValue.h"

namespace GeneTrail
{
namespace pvalue
{

pValueVec bonferroni_v(const pValueVec& pvalues)
{
	return bonferroni(pvalues, get_second());
}

pValueVec sidak_v(const pValueVec& pvalues)
{
	return sidak(pvalues, get_second());
}

pValueVec holm_v(const pValueVec& pvalues)
{
	return holm(pvalues, get_second());
}

pValueVec holm_sidak_v(const pValueVec& pvalues)
{
	return holm_sidak(pvalues, get_second());
}

pValueVec finner_v(const pValueVec& pvalues)
{
	return finner(pvalues, get_second());
}

pValueVec benjamini_hochberg_v(const pValueVec& pvalues)
{
	return benjamini_hochberg(pvalues, get_second());
}

pValueVec benjamini_yekutieli_v(const pValueVec& pvalues)
{
	return benjamini_yekutieli(pvalues, get_second());
}

pValueVec hochberg_v(const pValueVec& pvalues)
{
	return hochberg(pvalues, get_second());
}

pValueVec simes_v(const pValueVec& pvalues)
{
	return simes(pvalues, get_second());
}

pValueVec noop(const pValueVec& pvalues) { return pvalues; }

boost::optional<MultipleTestingCorrection>
getCorrectionMethod(const std::string& method)
{
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
		return boost::make_optional(
		    MultipleTestingCorrection::BenjaminiHochberg);
	}

	if(method == "benjamini_yekutieli" || method == "benjamini-yekutieli") {
		return boost::make_optional(
		    MultipleTestingCorrection::BenjaminiYekutieli);
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

using CorrectionMethod = pValueVec (*)(const pValueVec&);

CorrectionMethod getCorrectionMethod(MultipleTestingCorrection method)
{
	switch(method) {
		case MultipleTestingCorrection::Bonferroni:
			return bonferroni_v;
		case MultipleTestingCorrection::Sidak:
			return sidak_v;
		case MultipleTestingCorrection::Holm:
			return holm_v;
		case MultipleTestingCorrection::HolmSidak:
			return holm_sidak_v;
		case MultipleTestingCorrection::Finner:
			return finner_v;
		case MultipleTestingCorrection::BenjaminiHochberg:
			return benjamini_hochberg_v;
		case MultipleTestingCorrection::BenjaminiYekutieli:
			return benjamini_yekutieli_v;
		case MultipleTestingCorrection::Hochberg:
			return hochberg_v;
		case MultipleTestingCorrection::Simes:
			return simes_v;
		case MultipleTestingCorrection::GSEA:
			return noop;
		default:
			throw NotImplemented(__FILE__, __LINE__, "The requested correction "
			                                         "method has not yet been "
			                                         "implemented.");
	}
}

pValueVec adjustPValues(const pValueVec& pvalues,
                        MultipleTestingCorrection method)
{
	return getCorrectionMethod(method)(pvalues);
}

} // namespace pvalue
} // namespace GeneTrail2
