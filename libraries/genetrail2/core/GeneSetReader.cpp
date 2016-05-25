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
#include "GeneSetReader.h"

#include "Exception.h"
#include "GeneSet.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>

namespace GeneTrail
{
	template <typename Processor>
	GeneSet GeneSetReader::read(const std::string& path, Processor p,
	                            size_t numberOfElementPerLine,
	                            const std::string& delimiter) const
	{
		GeneSet gene_set;
		std::ifstream input(path);
		if(!input) {
			throw IOError("File (" + path + ") is not open for reading");
		}

		std::vector<std::string> sline;
		int l = 1;
		for(std::string line; getline(input, line); ++l) {
			boost::trim(line);

			if(line == "") {
				continue;
			}

			if(line.find("(class=") != std::string::npos) {
				continue;
			}

			boost::split(sline, line, boost::is_any_of(delimiter),
			             boost::token_compress_on);
			if(sline.size() == numberOfElementPerLine) {
				gene_set.insert(p(sline));
			} else {
				std::string err =
				    (sline.size() > numberOfElementPerLine) ? "many" : "few";
				throw IOError("Wrong file format: Line " +
				              boost::lexical_cast<std::string>(l) +
				              " contains too " + err + " elements");
			}
		}
		return gene_set;
	}

	GeneSet::Element scoringFileProcessor(std::vector<std::string>& sline)
	{
		boost::trim(sline[0]);
		boost::trim(sline[1]);
		return std::make_pair(sline[0],
		                      boost::lexical_cast<double>(sline[1]));
	}

	GeneSet::Element geneListProcessor(std::vector<std::string>& sline)
	{
		return std::make_pair(sline[0], 0.0);
	}

	GeneSet GeneSetReader::readScoringFile(const std::string& path) const
	{
		return read(path, scoringFileProcessor, 2, " \t");
	}

	GeneSet GeneSetReader::readGeneList(const std::string& path) const
	{
		return read(path, geneListProcessor, 1, " \t");
	}

	GeneSet GeneSetReader::readNAFile(const std::string& path) const
	{
		return read(path, scoringFileProcessor, 2, "=");
	}
}

