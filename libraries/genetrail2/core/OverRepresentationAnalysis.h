/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_OVER_REPRESENTATION_ANALYSIS_H
#define GT2_OVER_REPRESENTATION_ANALYSIS_H

#include "macros.h"

#include "Category.h"
#include "FishersExactTest.h"
#include "HypergeometricTest.h"

#include <utility>
#include <tuple>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/binomial.hpp>

namespace GeneTrail {

    /**
     * Over-representation analysis (ORA)
     */
    class GT2_EXPORT OverRepresentationAnalysis {

		public:

			/**
			 * This method checks if all genes in set are contained in category.
			 *
			 * @param category
			 * @param set
			 * @return 
			 */
			bool categoryContainsAllGenes(const Category& category, const Category& set);

			/**
			 * This method computes a one-sided p-value.
			 *
			 * @param category The category for which a p-value should be computed.
			 * @param reference_set Reference set
			 * @param test_set Test set
			 * @return P-value
			 */
			std::tuple<double, double, std::string> computePValue(const Category& category, const Category& reference_set, const Category& test_set);

	};
}

#endif // GT2_OVER_REPRESENTATION_ANALYSIS_H

