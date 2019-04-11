 
/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2014 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#ifndef GT2_REGULATION_REGULATION_FILE_WRITER_H
#define GT2_REGULATION_REGULATION_FILE_WRITER_H

#include <fstream>
#include <iostream>
#include <string>


namespace GeneTrail
{

	class GT2_EXPORT RegulationFileWriter2
	{
		public:

		RegulationFileWriter2(){};

		void write(const std::vector<std::tuple<std::string,std::string,double>> regulations, const std::string& path,
		           const std::string& delimiter) {
		
		std::ofstream output;
		output.open(path);

		if(!output) {
			throw IOError("File (" + path + ") is not open for reading");
		}

		for(const auto& p : regulations) {
			output << std::get<0>(p) << delimiter << std::get<1>(p) << delimiter << std::get<2>(p)<< std::endl;
		}

		output.close();  
			  }

	};
}

#endif //GT2_CORE_GENE_SET_WRITER_H

