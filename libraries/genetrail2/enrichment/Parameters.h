/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *				 2018 Tim Kehl <tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_ENRICHMENT_PARAMETERS_H
#define GT2_ENRICHMENT_PARAMETERS_H

#include "PermutationTest.h"

#include <genetrail2/core/MatrixHTest.h>
#include <genetrail2/core/PValue.h>

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

		const std::string& category() const {
			return category_.filePath;
		}

		const std::string& category_name() const {
			return category_name_;
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
		FilePath category_;
		FilePath scores_;
		FilePath identifier_;
		FilePath dataMatrixPath_;
		FilePath groups_;

		std::string category_name_;

		DirectoryPath out_;

		size_t minimum;
		size_t maximum;
		size_t numPermutations;
		size_t randomSeed;

		bool adjustSeparately;
		bool reducedOutput;
		bool verbose;

		boost::optional<MultipleTestingCorrection> adjustment;
		PValueMode pValueMode;
		boost::optional<MatrixHTests> scoringMethod;
	};
}

#endif // GT2_ENRICHMENT_PARAMETERS_H
