/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2019 Tim Kehl <tkehl@bioinf.uni-sb.de>
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
#include "ORAPreprocessor.h"
#include <iostream>

using namespace GeneTrail;

ORAPreprocessor::ORAPreprocessor(uint64_t m, uint64_t n, uint64_t l_max)
: m_(m), 
  n_(n),
  l_max_(l_max),
  p_values_(l_max , l_max)
{
    for(uint64_t i=0; i<l_max; ++i) {
        for(uint64_t j=0; j<l_max; ++j) {
            p_values_(i,j) = 1.0;
        }
    }
    fill_matrix();

}

const DenseMatrix& ORAPreprocessor::getPValues() {
    return p_values_;
}

void ORAPreprocessor::fill_matrix() {
    big_float denom = boost::math::binomial_coefficient<big_float>(m_, n_);
    for(uint64_t l=1; l<l_max_; ++l) {
        std::cout << "INFO: Calculating " << l << "/" << l_max_ << std::endl;
        big_float p_val = 0.0;
        uint64_t d = std::min(n_, l);
        for(uint64_t k=d; k>0; --k) {
            p_val += boost::math::binomial_coefficient<big_float>(l, k) *
				    boost::math::binomial_coefficient<big_float>(m_ - l, n_ - k);
            big_float p = (p_val / denom);
            p_values_(l,k) = p.convert_to<double>();;
        }
    }
}