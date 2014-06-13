#ifndef GT2_CORE_GENE_SET_WRITER_H
#define GT2_CORE_GENE_SET_WRITER_H

#include <fstream>
#include <iostream>
#include <string>

#include "Exception.h"
#include "GeneSet.h"

#include "macros.h"

namespace GeneTrail
{
	template<typename value_type>
	class GT2_EXPORT GeneSetWriter
	{
		public:

		GeneSetWriter(){};

		void write(GeneSet<value_type> gene_set, const std::string& path, const std::string& delimiter, const std::string& header = "")
		{
			std::ofstream output;
			output.open(path);

			if(!output)
			{
				throw IOError("File (" + path + ") is not open for reading");
			}

			if(header != "")
			{
				output << header << std::endl;
			}

			for(const auto& p : gene_set.getScores())
			{
				output << p.first << delimiter << p.second << std::endl;
			}

			output.close();
		}

		void writeNAFile(GeneSet<value_type> gene_set, const std::string& path)
		{
			write(gene_set, path, " = ", "Scores (class=Double)");
		}

		void writeScoringFile(GeneSet<value_type> gene_set, const std::string& path)
		{
			write(gene_set, path, "\t");
		}
	};
}

#endif //GT2_CORE_GENE_SET_WRITER_H

