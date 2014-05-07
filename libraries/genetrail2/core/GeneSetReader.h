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

		template<typename Processor>
		GeneSet<value_type> read(const std::string& path, Processor p, int numberOfElementPerLine, std::string delimiter = " \t")
		{
			GeneSet<value_type> gene_set;
			std::ifstream input(path);
			if(!input) {
				throw IOError("File (" + path + ") is not open for reading");
			}
			for(std::string line; getline(input, line);) {
				std::vector<std::string> sline;
				boost::split(sline, line, boost::is_any_of(delimiter));
				if(sline.size() == numberOfElementPerLine) {
				gene_set.insert(p(sline));
				}else{
					throw IOError("Wrong file format.");
				}
			}
			return gene_set;
		}

		static std::pair<std::string,value_type> scoringFileProcessor(std::vector<std::string> sline)
		{
			boost::trim(sline[0]);
			boost::trim(sline[1]);
			return std::make_pair(sline[0], boost::lexical_cast<value_type>(sline[1]));
		}

		GeneSet<value_type> readScoringFile(const std::string& path)
		{
			return read(path, scoringFileProcessor, 2);
		}

		static std::pair<std::string,value_type> geneListProcessor(std::vector<std::string> sline)
		{
			boost::trim(sline[0]);
			return std::make_pair(sline[0], (value_type) 0.0);
		}

		GeneSet<value_type> readGeneList(const std::string& path)
		{
			return read(path, geneListProcessor, 1);
		}

		GeneSet<value_type> readNAFile(const std::string& path)
		{
			return read(path, scoringFileProcessor, 2, "=");
		}
	};
}

#endif //GT2_CORE_GENE_SET_READER_H

