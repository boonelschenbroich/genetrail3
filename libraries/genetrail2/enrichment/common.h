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
#ifndef GT2_ENRICHMENT_COMMON_H
#define GT2_ENRICHMENT_COMMON_H

#include "CommandLineInterface.h"
#include "EnrichmentAlgorithm.h"

#include <genetrail2/core/MatrixHTest.h>
#include <genetrail2/core/macros.h>
#include <genetrail2/core/PValue.h>

#include <boost/program_options.hpp>

#include <list>
#include <utility>
#include <string>
#include <memory>
#include <random>

namespace bpo = boost::program_options;

using namespace GeneTrail;

namespace GeneTrail
{
	class GeneSet;

	struct DirectoryPath;
	struct EnrichmentResult;
	struct FilePath;
	struct Params;

	using EnrichmentAlgorithmPtr = std::unique_ptr<EnrichmentAlgorithm>;
}

typedef std::map<std::string, std::shared_ptr<EnrichmentResult>> Results;
typedef std::map<std::string, Results> AllResults;
typedef std::vector<std::pair<std::string, double>> PValueList;

using CategoryList = std::list<std::pair<std::string, std::string>>;


/**
 * This function initializes the needed attributes.
 *
 * @param test_set Test set to be filled
 * @param cat_list CategoryList to be filled
 * @param p Parameter object
 * @return -1 if an error occurred and 0 if not
 */
GT2_EXPORT int init(GeneSet& test_set, CategoryList& cat_list, const Params& p);

/**
 * This function runs the entire pipeline.
 *
 * @param test_set Test set for which the computation should be started
 * @param cat_list List of categories for the computation
 * @param p
 */
GT2_EXPORT void run(Scores& test_set, CategoryList& cat_list, EnrichmentAlgorithmPtr& algorithm, const Params& p, bool computePValues);

#endif // GT2_APPLICATIONS_ENRICHMENT_COMMON_H
