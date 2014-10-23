#ifndef GT2_CORE_GENE_SET_READER_H
#define GT2_CORE_GENE_SET_READER_H

#include <fstream>
#include <iostream>
#include <utility>
#include <string>
#include <functional>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

#include "Exception.h"
#include "Category.h"
#include "GeneSet.h"

#include "macros.h"

namespace GeneTrail
{
	template<typename value_type>
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
		GeneSet<value_type> read(const std::string& path, Processor p,
		                         size_t numberOfElementPerLine,
		                         const std::string& delimiter)
		{
			GeneSet<value_type> gene_set;
			std::ifstream input(path);
			if(!input) {
				throw IOError("File (" + path + ") is not open for reading");
			}
			int l = 1;
			for(std::string line; getline(input, line);) {
				if(line.find("(class=") != std::string::npos)
				{
					continue;
				}
				std::vector<std::string> sline;
				boost::split(sline, line, boost::is_any_of(delimiter), boost::token_compress_on);
				if(sline.size() == numberOfElementPerLine) {
					gene_set.insert(p(sline));
				} else {
					std::string err = (sline.size() > numberOfElementPerLine)
					                      ? "many"
					                      : "less";
					throw IOError("Wrong file format: Line " + boost::lexical_cast<std::string>(l) +
					              " contains too " + err + " elements");
				}
				++l;
			}
			return gene_set;
		}

		/**
		 * Processor function for a line in a scoring file.
		 *
		 * @param sline Splitted line to parse
		 */
		static std::pair<std::string, value_type>
		scoringFileProcessor(std::vector<std::string>& sline)
		{
			boost::trim(sline[0]);
			boost::trim(sline[1]);
			return std::make_pair(sline[0],
			                      boost::lexical_cast<value_type>(sline[1]));
		}

		/**
		 * Specialization of the generic read function for scoring files
		 *
		 * @param path Path to the file
		 */
		GeneSet<value_type> readScoringFile(const std::string& path)
		{
			return read(path, scoringFileProcessor, 2, " \t");
		}

		/**
		 * Processor function for a line in a gene list.
		 *
		 * @param sline Splitted line to parse
		 */
		static std::pair<std::string, value_type>
		geneListProcessor(std::vector<std::string>& sline)
		{
			boost::trim(sline[0]);
			return std::make_pair(sline[0], boost::lexical_cast<value_type>(0.0));
		}

		/**
		 * Specialization of the generic read function for gene lists
		 *
		 * @param path Path to the file
		 */
		GeneSet<value_type> readGeneList(const std::string& path)
		{
			return read(path, geneListProcessor, 1, " \t");
		}

		/**
		 * Specialization of the generic read function for na files
		 *
		 * @param path Path to the file
		 */
		GeneSet<value_type> readNAFile(const std::string& path)
		{
			return read(path, scoringFileProcessor, 2, "=");
		}
	};
}

#endif //GT2_CORE_GENE_SET_READER_H

