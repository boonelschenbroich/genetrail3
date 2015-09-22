#ifndef GT2_ENRICHMENT_PARAMETERS_H
#define GT2_ENRICHMENT_PARAMETERS_H

#include <genetrail2/core/MatrixHTest.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/PermutationTest.h>

#include <boost/optional/optional.hpp>

#include <string>

namespace GeneTrail {
	struct GT2_EXPORT FilePath
	{
		FilePath() {}
		explicit FilePath(const std::string& path) : filePath(path) {
		}

		std::string filePath;
	};

	struct GT2_EXPORT DirectoryPath
	{
		DirectoryPath() {}
		explicit DirectoryPath(const std::string& path) : directoryPath(path) {
		}

		std::string directoryPath;
	};

	struct GT2_EXPORT Params
	{
		Params();

		const std::string& categories() const {
			return categories_.filePath;
		}

		const std::string& scores() const {
			return scores_.filePath;
		}

		const std::string& identifier() const {
			return identifier_.filePath;
		}

		const std::string& dataMatrixPath() const {
			return dataMatrixPath_.filePath;
		}

		const std::string& groups() const {
			return groups_.filePath;
		}

		const std::string& out() const {
			return out_.directoryPath;
		}

		double significance;

		FilePath categories_;
		FilePath scores_;
		FilePath identifier_;
		FilePath dataMatrixPath_;
		FilePath groups_;

		DirectoryPath out_;

		size_t minimum;
		size_t maximum;
		size_t numPermutations;
		size_t randomSeed;

		bool adjustSeparately;

		boost::optional<MultipleTestingCorrection> adjustment;
		PValueMode pValueMode;
		boost::optional<MatrixHTests> scoringMethod;
	};
}

#endif // GT2_ENRICHMENT_PARAMETERS_H
