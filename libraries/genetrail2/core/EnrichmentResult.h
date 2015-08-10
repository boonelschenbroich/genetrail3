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
#ifndef GT2_CORE_ENRICHMENT_RESULT_H
#define GT2_CORE_ENRICHMENT_RESULT_H

#include "macros.h"
#include "multiprecision.h"

#include <string>
#include "Category.h"

#include <boost/lexical_cast.hpp>

namespace GeneTrail {

    /**
     * GeneralEnrichmentResult
     */
    struct GT2_EXPORT EnrichmentResult {
		EnrichmentResult(const std::shared_ptr<Category>& c) : category(c), hits(0), pvalue(-1.0), enriched(false), score(0.0)
		{}

		std::shared_ptr<Category> category;
		unsigned int hits;
		big_float pvalue;
		std::string info;
		bool enriched;
		double score;

		virtual std::string header() const{
			std::string header = "#";
			header += "Name\t";
			header += "Reference\t";
			header += "Hits\t";
			header += "Score\t";
			header += "P-value\t";
			header += "Info\t";
			header += "Regulation_direction";
			return header;	
		}

		virtual std::string serialize() const{
			std::string result = "";
			result += category->name() + "\t";
			result += category->reference() + "\t";
			result += boost::lexical_cast<std::string>(hits) + "\t";
			result += boost::lexical_cast<std::string>(score) + "\t";
			result += boost::lexical_cast<std::string>(pvalue) + "\t";
			result += info + "\t";
			result += boost::lexical_cast<std::string>(enriched);
			return result;
		}
	};
}

#endif //GT2_CORE_ENRICHMENT_RESULT_H

