 
/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2017 Lea Eckhart <leckhart@bioinf.uni-sb.de>
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

#ifndef GT2_MATRIX_TRANSFORMATION_H
#define GT2_MATRIX_TRANSFORMATION_H

#include <numeric> 
#include <cmath>

#include "macros.h"
#include "Matrix.h"
#include "Statistic.h"
#include "MatrixIterator.h"

namespace GeneTrail
{	
	  
	/**
	 * This function transforms values in a matrix to ranks for each column.
	 * Ranks start with 1 for the largest value in the column.
	 *
	 * @param matrix The Matrix whose values should be transformed into ranks.
	 * 
	 */
	  template <typename Matrix>
	  void valuesToRanks(Matrix& in_matrix) {

	     for(size_t c=0; c<in_matrix.cols(); ++c) {
	      std::vector<size_t> rankVec(in_matrix.rows(),0);
 	      std::iota(rankVec.begin(), rankVec.end(), 0);
	      std::sort(rankVec.begin(), rankVec.end(),
  	     [&](const size_t a,const size_t b)-> bool {return in_matrix(a,c) > in_matrix(b,c);});
	      for(size_t v=0; v<rankVec.size(); ++v) {
		in_matrix.set(rankVec[v],c,v+1);
	      }
	     }
	     return;
	  }
	  
	/**
	 * This function applies a function to each value of the Matrix
	 *
	 * @param matrix The Matrix to which the function should be applied.
	 * @param f The function that should be applied.
	 * 
	 */
	  template <typename Matrix, typename Func>
	  void transformMatrix(Matrix& matrix, Func f){
	    
	      for(size_t i=0; i<matrix.rows();++i) {
		for(size_t j=0; j<matrix.cols();++j) {
		  matrix.set(i,j,f(matrix(i,j)));
		}
	      }
	  }
	  
 
	/**
	 * This function applies |r/2 + 0.5 - value| to each value of the Matrix,
	 * where r denotes the number of rows.
	 * 
	 * The approach is based on the paper
	 * GSVA: gene set variation analysis for microarray and RNA-Seq data
	 * by Hänzelmann et al. (DOI: 10.1186/1471-2105-14-7)
	 *
	 * @param matrix The Matrix to which the function should be applied.
	 * 
	 */
	template <typename Matrix> 
	 void upweightEnds(Matrix& matrix) {
	   transformMatrix(matrix, [&](double d) {return std::abs((matrix.rows()/2) + 0.5 -d); });
	 }
	 
	/**
	 * This function applies r/2 + 0.5 - value to each value of the Matrix,
	 * where r denotes the number of rows.
	 * 
	 * The approach is based on the paper
	 * GSVA: gene set variation analysis for microarray and RNA-Seq data
	 * by Hänzelmann et al. (DOI: 10.1186/1471-2105-14-7)
	 *
	 * @param matrix The Matrix to which the function should be applied.
	 * 
	 */
	template <typename Matrix> 
	 void upweightTail(Matrix& matrix) {
	   transformMatrix(matrix, [&](double d) {return (matrix.rows()/2) + 0.5 -d; });
	 }

  
	/**
	 * This function applies std::abs to all values of the Matrix.
	 *
	 * @param matrix The Matrix to which the function should be applied.
	 * 
	 */
	template <typename Matrix>
	void abs(Matrix& matrix) {
	  transformMatrix(matrix, [](double d) { return std::abs(d); });
	}


	/**
	 * This function applies std::sqrt to all values of the Matrix.
	 *
	 * @param matrix The Matrix to which the function should be applied.
	 */
	template <typename Matrix>
	void sqrt(Matrix& matrix) {
	  transformMatrix(matrix, [](double d) { return std::sqrt(d); });
	}

	
	/**
	 * This function applies std::log to all values of the Matrix.
	 *
	 * @param matrix The Matrix to which the function should be applied.
	 */
	template <typename Matrix>
	void log(Matrix& matrix) {
	  transformMatrix(matrix, [](double d) { return std::log(d); });
	}

	/**
	 * This function applies std::log2 to all values of the Matrix.
	 *
	 * @param matrix The Matrix to which the function should be applied.
	 */
	template <typename Matrix>
	void log2(Matrix& matrix) {
	  transformMatrix(matrix, [](double d) { return std::log2(d); });
	}

	/**
	 * This function applies std::log10 to all values of the Matrix.
	 *
	 * @param matrix The Matrix to which the function should be applied.
	 */
	template <typename Matrix>
	void log10(Matrix& matrix) {
	  transformMatrix(matrix, [](double d) { return std::log10(d); });
	}
	
	/**
	 * This function applies std::pow to all values of the Matrix.
	 *
	 * @param matrix The Matrix to which the function should be applied.
	 * @param n The exponent that should be used in the power function.
	 */
	template <typename Matrix>
	void pow(Matrix& matrix, int n) {
	 transformMatrix(matrix, [&](double d) { return std::pow(d,n); });
	}

	/**
	 * This function applies std::pow (power = 2) to all values of the
	 * Matrix.
	 *
	 * @param matrix The Matrix to which the function should be applied.
	 */
	template <typename Matrix>
	void pow2(Matrix& matrix) {
	  pow(matrix,2);
	}
  
  
  
}

#endif //GT2_MATRIX_TRANSFORMATION_H