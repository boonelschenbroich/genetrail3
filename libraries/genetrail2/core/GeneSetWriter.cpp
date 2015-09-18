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

