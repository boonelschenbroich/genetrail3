/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef FIDEPA_H
#define FIDEPA_H

#include "config.h"

#include "ScoringFile.h"
#include "GeneSetReader.h"
#include "BoostGraphParser.h"
#include "BoostGraphProcessor.h"
#include "Path.h"
#include "Pathfinder.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <string>
#include <cstring>
#include <tuple>
#include <fstream>
#include <set>
#include <vector>
#include <iostream>
#include <ostream>
#include <sstream>
#include <algorithm>


namespace GeneTrail {

    class GT2_EXPORT FiDePaRunner {
    public:

        FiDePaRunner() {
        };

        ~FiDePaRunner() {
        };

        std::string convertInt(int number) const;

        void writeSifFiles(std::vector<std::vector<std::string> > best_paths, std::map<std::string, std::string> regulations, std::string scores) const;

        void computeDeregulatedPaths(std::string kegg, std::string scores, int pathlength, bool descending, bool absolute) const;
    };
}

#endif

