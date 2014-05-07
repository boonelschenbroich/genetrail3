/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_ORA_RESULT_H
#define GT2_CORE_ORA_RESULT_H

#include "EnrichmentResult.h"

#include "macros.h"

namespace GeneTrail {

    /**
     * EnrichmentResult
     */
    struct GT2_EXPORT ORAResult : public EnrichmentResult {

		ORAResult() : expected_hits(0.0)
		{}

		double expected_hits;

		std::string serialize() const{
			std::string result = "";
			result += boost::lexical_cast<std::string>(name) + "\t";
			result += boost::lexical_cast<std::string>(reference) + "\t";
			result += boost::lexical_cast<std::string>(expected_hits) + "\t";
			result += boost::lexical_cast<std::string>(hits) + "\t";
			result += boost::lexical_cast<std::string>(pvalue) + "\t";
			result += boost::lexical_cast<std::string>(info) + "\t";
			result += boost::lexical_cast<std::string>(enriched) + "\t";
			return result;
		}
	};
}

#endif //GT2_CORE_ORA_RESULT_H

