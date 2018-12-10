/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *               2014-2018 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_ENRICHMENT_COMMANDLINEINTERFACE_H
#define GT2_ENRICHMENT_COMMANDLINEINTERFACE_H

#include "PermutationTest.h"

#include <genetrail2/core/macros.h>
#include <genetrail2/core/PValue.h>

#include <boost/optional/optional.hpp>

#include <string>
#include <vector>

namespace boost {
	class any;
	namespace program_options {
		class options_description;
	}
}

namespace GeneTrail {
	struct DirectoryPath;
	struct FilePath;
	struct Params;

	/**
	 * This function adds common arguments to the BOOST commandline parser.
	 *
	 * @param desc BOOST program-options description
	 * @param p Parameter object
	 */
	GT2_EXPORT void addCommonCLIArgs(boost::program_options::options_description& desc, Params& p);

	GT2_EXPORT bool checkCLIArgs(const Params& p);

	GT2_EXPORT void validate(boost::any&, const std::vector<std::string>&, PValueMode*,
	              int);

	GT2_EXPORT void validate(boost::any&, const std::vector<std::string>&, boost::optional<MultipleTestingCorrection>*,
	              int);

	GT2_EXPORT void validate(boost::any&, const std::vector<std::string>&, boost::optional<MatrixHTests>*,
	              int);

	GT2_EXPORT void validate(boost::any&, const std::vector<std::string>&, FilePath*,
	              int);

	GT2_EXPORT void validate(boost::any&, const std::vector<std::string>&, DirectoryPath*,
	              int);
}

#endif // GT2_ENRICHMENT_COMMANDLINEINTERFACE_H
