/*
 * GeneTrail2 - An efficient library for interpreting genetic data
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
#ifndef GT2_CORE_OVER_REPRESENTATION_ANALYSIS_H
#define GT2_CORE_OVER_REPRESENTATION_ANALYSIS_H

#include "macros.h"

#include "Category.h"
#include "FishersExactTest.h"
#include "HypergeometricTest.h"
#include "multiprecision.h"

#include <utility>
#include <tuple>

#include <boost/math/special_functions/binomial.hpp>

namespace GeneTrail {

    /**
     * Over-representation analysis (ORA)
     */
    class GT2_EXPORT OverRepresentationAnalysis {

		private:

			/**
			 * This method checks if all genes in set are contained in category.
			 *
			 * @param category
			 * @param set
			 * @return 
			 */
			bool categoryContainsAllGenes(const Category& category, const Category& set);

		public:

			OverRepresentationAnalysis() = default;

			OverRepresentationAnalysis(const Category& reference_set,const Category& test_set);

			/**
			 * This method computes a one-sided p-value.
			 *
			 * @param category The category for which a p-value should be computed.
			 *
			 * @return P-value
			 */
			double computePValue(const Category& category) const;

			double computeScore(const Category& category) const;

			double numberOfHits(const Category& category) const;

			double expectedNumberOfHits(const Category& category) const;
		private:

			Category reference_set_;
			//Size of reference set
			size_t m_;

			Category test_set_;
			//Size of test set
			size_t n_;

			bool useHypergeometricTest_;

			FishersExactTest<uint64_t, big_float> fisherTest_;
			HypergeometricTest<uint64_t, big_float> hyperTest_;

	};
}

#endif // GT2_CORE_OVER_REPRESENTATION_ANALYSIS_H

