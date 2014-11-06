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

