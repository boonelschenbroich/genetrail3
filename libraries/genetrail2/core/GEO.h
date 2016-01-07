/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_GEO_H
#define GT2_GEO_H

#include <vector>
#include <set>
#include <string>
#include <map>
#include <list>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <functional>
#include <initializer_list>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>

#include "Statistic.h"
#include "macros.h"

namespace GeneTrail
{

	struct GT2_EXPORT GEOMap
	{
		std::map<std::string, std::vector<double>> gene2exprs;
		std::vector<std::string> sampleNames;
		std::string dataset = "";
		std::string platform = "";
	};

	class GT2_EXPORT GEO
	{
		public:

		GEO();

		~GEO();

		std::vector<double> mergeDuplicatedVectors(std::vector<std::vector<double>>, std::string);

		std::map<std::string, std::vector<double>> mapAndRemoveDuplicates(std::map<std::string, std::vector<double>>, std::map<std::string, std::string>, std::string);

		double apply(std::string method, std::vector<double> tmp);

		void writeGEOMap(const std::string& filename, const GEOMap& map);

		private:
		typedef std::vector<double>::iterator _viter;
	};
}

#endif //GT2_GEO_H
