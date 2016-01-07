/*
 * GeneTrail2 - An efficent library for interpreting genetic data
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
#ifndef GT2_CORE_GENE_SET_READER_H
#define GT2_CORE_GENE_SET_READER_H

#include <string>

#include "macros.h"

namespace GeneTrail
{
	class GeneSet;

	class GT2_EXPORT GeneSetReader
	{
		public:

		GeneSetReader(){};

		/**
		 * Generic function to parse all kinds of files.
		 *
		 * @param path Path to file
		 * @param p Function to process each line of the file
		 * @param numberOfElementPerLine Allowed number of elements per line
		 * @param delimiter The delimiter to split each line (default is " \t")
		 */
		template <typename Processor>
		GeneSet read(const std::string& path, Processor p,
		             size_t numberOfElementPerLine,
		             const std::string& delimiter) const;

		/**
		 * Specialization of the generic read function for scoring files
		 *
		 * @param path Path to the file
		 */
		GeneSet readScoringFile(const std::string& path) const;

		/**
		 * Specialization of the generic read function for gene lists
		 *
		 * @param path Path to the file
		 */
		GeneSet readGeneList(const std::string& path) const;

		/**
		 * Specialization of the generic read function for na files
		 *
		 * @param path Path to the file
		 */
		GeneSet readNAFile(const std::string& path) const;
	};
}

#endif // GT2_CORE_GENE_SET_READER_H

