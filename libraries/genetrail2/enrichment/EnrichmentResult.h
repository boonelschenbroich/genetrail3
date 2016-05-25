/*
 * GeneTrail2 - An efficient library for interpreting genetic data
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

#include <genetrail2/core/Category.h>
#include <genetrail2/core/macros.h>
#include <genetrail2/core/multiprecision.h>

#include <ostream>
#include <string>

namespace GeneTrail
{
	/**
	 * GeneralEnrichmentResult
	 */
	struct GT2_EXPORT EnrichmentResult
	{
		EnrichmentResult(const std::shared_ptr<Category>& c)
		    : category(c),
		      hits(0),
		      pvalue(1.0),
		      enriched(false),
		      score(0.0),
		      expected_score(0.0)
		{
		}

		std::shared_ptr<Category> category;
		unsigned int hits;
		big_float pvalue;
		std::string info;
		bool enriched;
		double score;
		double expected_score;

		virtual std::string header() const
		{
			std::string header = "#";
			header += "Name\t";
			header += "Reference\t";
			header += "Hits\t";
			header += "Score\t";
			header += "Expected Score\t";
			header += "P-value\t";
			header += "Info\t";
			header += "Regulation_direction";
			return header;
		}

		virtual void serialize(std::ostream& strm) const
		{
			strm << category->name() << '\t'
			     << category->reference() << '\t'
			     << hits << '\t'
			     << score << '\t'
			     << expected_score << '\t'
			     << pvalue << '\t'
			     << info << '\t'
			     << enriched;
		}
	};

	inline std::ostream& operator<<(std::ostream& strm, const EnrichmentResult& result)
	{
		result.serialize(strm);
		return strm;
	}

	using EnrichmentResultPtr = std::shared_ptr<EnrichmentResult>;
	using EnrichmentResults = std::vector<EnrichmentResultPtr>;
}

#endif // GT2_CORE_ENRICHMENT_RESULT_H
