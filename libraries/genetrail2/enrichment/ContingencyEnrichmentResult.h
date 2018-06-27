/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2018 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_CONTINGENCE_ORA_RESULT_H
#define GT2_CONTINGENCE_ORA_RESULT_H

#include <string>
#include <vector>
#include <tuple>

#include <genetrail2/core/ConfidenceInterval.h>
#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/macros.h>

namespace GeneTrail
{
struct GT2_EXPORT ContingencyEnrichmentResult
{
	std::string name = "";
	size_t hits = 0;
	double expected_hits = 0;
	double p_value = 1.0;
	double corrected_p_value = 1.0;

	template <typename Writer> void serializeJSON(Writer& writer)
	{
		writer.StartObject();

		writer.String("name");
		writer.String(name.c_str());

		writer.String("hits");
		writer.Int(hits);

		writer.String("expectedHits");
		writer.Double(expected_hits);

		writer.String("pValue");
		writer.Double(p_value);

		writer.String("correctedPValue");
		writer.Double(corrected_p_value);
		
		writer.EndObject();
	}
};
}

#endif // GT2_CONTINGENCE_ORA_RESULT_H
