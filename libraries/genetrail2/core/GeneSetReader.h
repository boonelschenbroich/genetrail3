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

