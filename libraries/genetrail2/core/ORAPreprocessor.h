/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2014-2019 Tim Kehl <tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_ORA_PREPROCESSOR_H
#define GT2_CORE_ORA_PREPROCESSOR_H

#include "DenseMatrix.h"
#include "macros.h"
#include "multiprecision.h"

#include <boost/math/special_functions/binomial.hpp>

namespace GeneTrail {

    /**
     * ORAPreprocessor
     */
    class GT2_EXPORT ORAPreprocessor {
		public:

			ORAPreprocessor() = default;

			ORAPreprocessor(uint64_t m, uint64_t n, uint64_t l_max);

            const DenseMatrix& getPValues();

        private:

            void fill_matrix();

            // Reference size
            uint64_t m_;
            // Test set size
            uint64_t n_;
            uint64_t l_max_;

            DenseMatrix p_values_;
    };
}

#endif // GT2_CORE_ORA_PREPROCESSOR_H