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

#ifndef GT2_Path_H
#define GT2_Path_H

#include "config.h"

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

namespace GeneTrail
{
	class GT2_EXPORT Path
	{
		public:
			Path(){};
			
			void writeToSIFFile(const std::string& file);
			
			void addVertex(const std::string& vertex);
			
			std::string getVertex(const int& position);
			
			void addRegulation(const std::string& v1, const std::string& v2, const std::string& regulation);
			
			std::string getRegulation(const std::string& v1, const std::string& v2);
			
			int length();
			
			int runningSum();
			
			double pValue();

		protected:
                    std::vector<std::string> identifier_;
                    std::map<std::string, std::string> regulations_;
                    int runningSum_;
                    double pValue_;
	};
}

#endif //GT2_Path_H

