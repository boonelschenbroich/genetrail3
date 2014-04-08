#ifndef GT2_CORE_GENE_SET_READER_H
#define GT2_CORE_GENE_SET_READER_H

#include <fstream>
#include <iostream>
#include <utility>
#include <string>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

#include "Exception.h"
#include "Category.h"
#include "ScoringFile.h"

#include "config.h"

namespace GeneTrail
{
	class GT2_EXPORT GeneSetReader
	{
		public:
		GeneSetReader(){};

		Category readGeneSet(const std::string& path, const std::string name);

		std::vector<std::string> readGeneList(const std::string& path);

		template<typename value_type>
		ScoringFile<value_type> readScoringFile(const std::string& path){
			ScoringFile<value_type> scores;
			std::ifstream input(path);
			if(!input) {
				std::cerr << "Could not open " << path << " for reading." << std::endl;
			}
			for(std::string line; getline(input, line);) {
				std::vector<std::string> sline;
				boost::split(sline, line, boost::is_any_of(" \t"));
				if(sline.size() == 2) {
					boost::trim(sline[0]);
					boost::trim(sline[1]);
					scores.add(std::make_pair(sline[0], boost::lexical_cast<value_type>(sline[1])));
				}
			}
			return scores;
		}

		template <typename value_type>
		ScoringFile<value_type> readNAFile(const std::string& path)
		{
			ScoringFile<value_type> scores;
			std::ifstream input(path);
			if(!input) {
				std::cerr << "Could not open " << path << " for reading."
				          << std::endl;
			}
			for(std::string line; getline(input, line);) {
				if(line.find("class=") != std::string::npos)
					continue;
				std::vector<std::string> sline;
				boost::split(sline, line, boost::is_any_of("="));
				if(sline.size() == 2) {
					boost::trim(sline[0]);
					boost::trim(sline[1]);
					scores.add(std::make_pair(
					    sline[0], boost::lexical_cast<value_type>(sline[1])));
				}
			}
			return scores;
		}
	};
}

#endif //GT2_CORE_GENE_SET_READER_H

