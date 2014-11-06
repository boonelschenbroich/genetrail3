#ifndef GT2_CORE_GENE_SET_WRITER_H
#define GT2_CORE_GENE_SET_WRITER_H

#include <fstream>
#include <iostream>
#include <string>

#include "macros.h"

namespace GeneTrail
{
	class GeneSet;

	class GT2_EXPORT GeneSetWriter
	{
		public:

		GeneSetWriter(){};

		void write(const GeneSet& gene_set, const std::string& path,
		           const std::string& delimiter,
		           const std::string& header = "") const;

		void writeNAFile(const GeneSet& gene_set,
		                 const std::string& path) const;

		void writeScoringFile(const GeneSet& gene_set,
		                      const std::string& path) const;
	};
}

#endif //GT2_CORE_GENE_SET_WRITER_H

