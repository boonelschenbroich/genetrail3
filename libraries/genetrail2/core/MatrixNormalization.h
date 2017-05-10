 
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
#ifndef GT2_CORE_MATRIX_NORMALIZATION_H
#define GT2_CORE_MATRIX_NORMALIZATION_H

#include "macros.h"
#include "Matrix.h"
#include "Statistic.h"
#include "MatrixIterator.h"
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>

namespace GeneTrail
{
  	/**
	 * A class that performs gaussian kernel estimation 
	 * of the cumulative density function
	 */
	class GT2_EXPORT GaussEstimator 
	{
		public:
		
		 /**
		 * This method normalizes a value by using a gaussian kernel estimation 
		 * of the cumulative density function
		 * 
		 * The approach is based on the paper
		 * GSVA: gene set variation analysis for microarray and RNA-Seq data
		 * by Hänzelmann et al. (DOI: 10.1186/1471-2105-14-7)
		 *
		 * @param value InputIterator corresponding to the value to be normalized.
		 * @param begin	InputIterator corresponding to the begin of expression data also containing 'value'.
		 * @param end	InputIterator corresponding to the end of expression data.
		 *
		 * @return normalized value
		 */
		  template<typename InputIterator>
		  double normalizeValue (InputIterator value, InputIterator begin, InputIterator end) const {
		    
		    boost::math::normal normal;
		    double normalizedValue = 0.0;
		    double stddeviation = statistic::sd<double>(begin,end);
		    
		    if(stddeviation == 0) {
		      normalizedValue = std::distance(begin,end)*cdf(normal,0);
		    }
		    
		    else {
		      for(InputIterator iter = begin; iter != end; ++iter) {
		       normalizedValue += cdf(normal, ((*value)-(*iter))/(stddeviation/4));
		      }
		    }
		    
		    return normalizedValue/std::distance(begin,end);
		  }
	};
	
	 /**
	 * A class that performs poisson kernel estimation 
	 * of the cumulative density function
	 */
	class GT2_EXPORT PoissonEstimator 
	{
		public:
		
		 /**
		 * This method normalizes a value by using a poisson kernel estimation 
		 * of the cumulative density function
		 * 
		 * The approach is based on the paper
		 * GSVA: gene set variation analysis for microarray and RNA-Seq data
		 * by Hänzelmann et al. (DOI: 10.1186/1471-2105-14-7)
		 *
		 * @param value InputIterator corresponding to the value to be normalized.
		 * @param begin	InputIterator corresponding to the begin of discrete expression data also containing 'value'.
		 * @param end	InputIterator corresponding to the end of discrete expression data.
		 *
		 * @return normalized value
		 */
		template<typename InputIterator>
		double normalizeValue (InputIterator value, InputIterator begin, InputIterator end) const {
		    
		    double normalizedValue = 0.0;
		    
		    for(InputIterator iter = begin; iter != end; ++iter) {
		      boost::math::poisson poisson((*iter)+0.5);
		      normalizedValue += cdf(poisson, *value);
		    }
	
		    return normalizedValue/std::distance(begin,end);
		 }
	};
	
	/**
	 * A class that normalizes by applying an exclusive Z-Score.
	 */
	class GT2_EXPORT ExclusiveZScore
	{
		public:

		/**
		 * This method normalizes a value by using an exclusive Z-Score
		 * Score for x:
		 * The mean of all samples, excluding x, gets substracted from x.
		 * The result is then divided by the standard deviation of all samples, excluding x.
		 *
		 * @param value InputIterator corresponding to the value to be normalized.
		 * @param begin	InputIterator corresponding to the begin of expression data also containing 'value'.
		 * @param end	InputIterator corresponding to the end of expression data.
		 *
		 * @return exclusive Z-Score for 'value'
		 */
		template <typename InputIterator>
		double normalizeValue (InputIterator value, InputIterator begin, InputIterator end) const {
		  
		  int distance = std::distance(begin,end);
		  
		  double normalizedValue = *value;
		  
		  double mean = statistic::mean<double>(begin, end);
		  mean *= distance;
		  mean -= *value;
		  //mean of all samples but the considered one
		  mean /= (distance -1);
		  
		  double var = 0.0;
		  
		  for(InputIterator iter = begin; iter!= end; ++iter) {
		    var += (*iter - mean) * (*iter - mean);
		  }
		  
		  //sum of differences should not include considered sample
		  var -= (*value - mean) * (*value - mean);
		  var /= (distance -2);
		  double sd = std::sqrt(var);

		  return (sd == 0) ? 0 : (normalizedValue - mean)/sd;
		}

	};
	
	/**
	 * A class that normalizes rows of a matrix
	 * - by performing kernel estimation of the cumulative density function
	 *   for different (non-)discrete distributions
	 * - or by computing exclusive Z-Scores
	 */
	class GT2_EXPORT MatrixNormalization
	{
		public:
		
		 /**
		 * This method normalizes data by using
		 * - a kernel estimation of the cumulative density function
		 * - or an exclusive Z-Score
		 * and writes the result in a given row of a given Matrix
		 * 
		 * The approach (excluding the Z-Score) is based on the paper
		 * GSVA: gene set variation analysis for microarray and RNA-Seq data
		 * by Hänzelmann et al. (DOI: 10.1186/1471-2105-14-7)
		 *
		 * @param in_begin	InputIterator corresponding to the begin of the data to be normalized
		 * @param in_end	InputIterator corresponding to the end of the data to be normalized
		 * @param matrix	Matrix in which normalized values will be written
		 * @param row_index	Row of 'matrix' in which normalized values will be written
		 * @param norm		The distribution (gaussian, poisson) used for kernel density estimation,
		 * 			or the exclusive Z-Score
		 */
		 template<typename InputIterator, typename Matrix, typename Normalization> 
		 void normalizeRange(InputIterator in_begin, InputIterator in_end, Matrix& matrix, int row_index,  Normalization& norm) {
		   
		   int col_index = 0;
		   
		   // Iterate over one row (samples of one gene)
		   for(InputIterator iter = in_begin ;iter != in_end; ++iter) {
		     
		     matrix.set(row_index, col_index,norm.normalizeValue(iter, in_begin, in_end));
		     ++col_index;
		   }
		 }
		 
 		 /**
		 * This method normalizes a matrix by using
		 * - a kernel estimation of the cumulative density function for each row
		 * - or an exclusive Z-Score
		 * 
		 * The approach (excluding the Z-Score) is based on the paper
		 * GSVA: gene set variation analysis for microarray and RNA-Seq data
		 * by Hänzelmann et al. (DOI: 10.1186/1471-2105-14-7)
		 *
		 * @param matrix	The matrix to be normalized.
		 * @param norm		The distribution (gaussian, poisson) used for kernel density estimation,
		 * 			or the exclusive Z-Score
		 *
		 * @return normalized matrix
		 */
		 template<typename Normalization, typename Matrix> 
		 Matrix normalizeMatrix( const Matrix& matrix,  Normalization& norm) {
		   
		  Matrix normalizedMatrix(matrix.rowNames(), matrix.colNames());
		   
		  RowMajorMatrixIterator<Matrix> row_it(&matrix, 0);
		  RowMajorMatrixIterator<Matrix> row_end(&matrix, matrix.rows());
		   
		   int row_index = 0;
		   
		  // Iterate over each row (each gene)
		  for(; row_it != row_end; ++row_it) {
		     // compute normalization for each sample of the gene
		     normalizeRange(row_it->begin(), row_it->end(), normalizedMatrix, row_index, norm);
		     ++row_index;
		   }
		   
		   return normalizedMatrix;
		 }
	};
  
}

#endif //GT2_CORE_MATRIX_NORMALIZATION_H
