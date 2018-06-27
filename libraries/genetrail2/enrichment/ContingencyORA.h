/**
* GeneTrail2 - An efficent library for interpreting genetic data
* Copyright (C) 2018 Tim Kehl tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_CONTINGENCY_ORA_H
#define GT2_CONTINGENCY_ORA_H

#include <genetrail2/core/macros.h>
#include <genetrail2/core/Exception.h>
#include <genetrail2/core/OverRepresentationAnalysis.h>
#include <genetrail2/core/FishersExactTest.h>
#include <genetrail2/core/multiprecision.h>

#include "ContingencyEnrichmentResult.h"

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <vector>
#include <functional>
#include <tuple>

namespace GeneTrail
{
template <typename IntType, typename FloatType> class GT2_EXPORT ContingencyORA
{
  public:
	using int_type = IntType;
	using float_type = FloatType;
	using contingency_table = std::tuple<std::string, int_type, int_type, int_type, int_type>;
	using contingency_tables = std::vector<contingency_table>;

	ContingencyORA(const contingency_tables& ctables)
	:ctables_(ctables)
	{}

	std::vector<ContingencyEnrichmentResult> run(bool upperTailed = true){
		std::vector<ContingencyEnrichmentResult> results;
		results.reserve(ctables_.size());
		for(auto& ctable : ctables_){
			ContingencyEnrichmentResult result;
			std::string name;
			int_type a, b, c, d;
			std::tie(name, a, c, b ,d) = ctable;
			
			// Check if upper-tailed or lower-tailed p-value should be calculated
			auto k_ = (((double)(a+c)*b)/((double)(b+d)));
			float_type p;
			if(upperTailed || (k_ < a)) {
				p = test_.upperTailedPValue(b+d, b, a+c, a);
			} else {
				p = test_.lowerTailedPValue(b+d, b, a+c, a);
			}

			result.name = name;
			result.hits = a;
			result.expected_hits = k_;
			result.p_value = p;
			results.emplace_back(result);
		}
		return std::move(results);
	}

  private:
	const contingency_tables& ctables_;
	FishersExactTest<int_type, float_type> test_;
};
}

#endif // GT2_CONTINGENCY_ORA_H