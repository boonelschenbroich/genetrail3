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

} // namespace pvalue
} // namespace GeneTrail2
