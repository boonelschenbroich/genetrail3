/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#include "GeneSetWriter.h"

#include "Exception.h"
#include "GeneSet.h"
#include "Scores.h"

namespace GeneTrail
{
	void GeneSetWriter::write(const Scores& scores, const std::string& path,
	                          const std::string& delimiter,
	                          const std::string& header) const
	{
		std::ofstream output;
		output.open(path);

		if(!output) {
			throw IOError("File (" + path + ") is not open for reading");
		}

		if(header != "") {
			output << header << std::endl;
		}

		for(const auto& p : scores) {
			output << p.name(*scores.db()) << delimiter << p.score() << std::endl;
		}

		output.close();
	}

	void GeneSetWriter::writeNAFile(const Scores& gene_set,
	                                const std::string& path) const
	{
		write(gene_set, path, " = ", "Scores (class=Double)");
	}

	void GeneSetWriter::writeScoringFile(const Scores& gene_set,
	                                     const std::string& path) const
	{
		write(gene_set, path, "\t");
	}
}

