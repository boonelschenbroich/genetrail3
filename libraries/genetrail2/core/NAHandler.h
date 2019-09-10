/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_NA_HANDLER_H
#define GT2_CORE_NA_HANDLER_H

#include "DenseMatrix.h"

#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <iostream>

namespace GeneTrail{

struct GT2_EXPORT NAStrategyZero{
	template <typename Matrix> bool handle(Matrix& m, unsigned int r) const{
		using T = typename std::remove_reference<decltype (m(0,0))>::type;
		for(unsigned int i=0; i < m.cols(); i++){
			if(std::isnan(m(r,i))) m(r,i) = static_cast<T>(0.0);
		}
		return true;
	}
};

struct GT2_EXPORT NAStrategyMean{
	template <typename Matrix> bool handle(Matrix& m, unsigned int r) const{
		using T = typename std::remove_reference<decltype (m(0,0))>::type;
		T mean = static_cast<T>(0.0);
		unsigned int c = m.cols();
		unsigned int c_na = 0;
		for(unsigned int i=0; i < c; i++){
			if(!std::isnan(m(r,i))){
				mean += m(r,i);
				c_na++;
			}
		}
		mean /= c_na;
		
		for(unsigned int i=0; i < c; i++){
			if(std::isnan(m(r,i))) m(r,i) = mean;
		}
		return true;
	}
};

struct GT2_EXPORT NAStrategyMedian{
	template <typename Matrix> bool handle(Matrix& m, unsigned int r) const{
		using T = typename std::remove_reference<decltype (m(0,0))>::type;
		std::vector<T> tmp;
		for(unsigned int i=0; i < m.cols(); i++){
			if(!std::isnan(m(r,i))) tmp.push_back(m(r,i));
		}
		std::sort(tmp.begin(), tmp.end());
		size_t size = tmp.size();
		T median;
		if(size % 2 == 0){
		median = (tmp[size/2-1] + tmp[size/2]) / static_cast<T>(2.0);
		} else{ 
		median = tmp[size/2];
		}
		for(unsigned int i=0; i < m.cols(); i++){
			if(std::isnan(m(r,i))) m(r,i) = median;
		}
		return true;
	}
};

struct GT2_EXPORT NAStrategyRemove{
	template <typename Matrix> bool handle(Matrix& m, unsigned int r) const{
		for(unsigned int i=0; i < m.cols(); i++){
			if(std::isnan(m(r,i))) return false;
		}
		return true;
	}
};

struct GT2_EXPORT NAHandler{
	/**
	 * Removes NAs in a file by applying the given removal
	 * strategy.
	 */
	template <typename Matrix, typename Strategy>
	void handle(Matrix& input, const Strategy& strategy){
		std::vector<unsigned int> toBeRemoved;
		for(unsigned int r=0; r < input.rows(); r++){
			if(!strategy.handle(input, r)){
				toBeRemoved.push_back(r);
			}
		}
		if(!toBeRemoved.empty()) input.removeRows(toBeRemoved);
	}

};
} // namespace GeneTrail

#endif // GT2_CORE_NA_HANDLER_H
