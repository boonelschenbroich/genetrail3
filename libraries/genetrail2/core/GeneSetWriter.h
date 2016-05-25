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
#ifndef GT2_CORE_GENE_SET_WRITER_H
#define GT2_CORE_GENE_SET_WRITER_H

#include <fstream>
#include <iostream>
#include <string>

#include "macros.h"

namespace GeneTrail
{
	class GeneSet;
	class Scores;

	class GT2_EXPORT GeneSetWriter
	{
		public:

		GeneSetWriter(){};

		void write(const Scores& gene_set, const std::string& path,
		           const std::string& delimiter,
		           const std::string& header = "") const;

		void writeNAFile(const Scores& gene_set,
		                 const std::string& path) const;

		void writeScoringFile(const Scores& gene_set,
		                      const std::string& path) const;
	};
}

#endif //GT2_CORE_GENE_SET_WRITER_H

