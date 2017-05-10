#include <gtest/gtest.h>
#include <config.h>
#include <iterator>
#include <vector>
#include <iostream>
#include <genetrail2/core/MatrixNormalization.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/MatrixIterator.h>
#include <genetrail2/core/Statistic.h>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>

using namespace GeneTrail;

double TOLERANCE = 0.00001;

template <typename InputIterator>
double exclusiveZ(double value, InputIterator begin, InputIterator end) {
  
  std::vector<double> newVec;
  bool notErased = true;
  
  for(; begin != end; ++ begin) {
    if ((*begin) == value && notErased) {
      notErased = false;
      continue;
    }
    else
      newVec.push_back(*begin);
  }

  double mean = statistic::mean<double>(newVec.begin(), newVec.end());
  double sd = statistic::sd<double>(newVec.begin(), newVec.end());
  
  return (sd==0) ? 0 : (value - mean)/sd;
}

/**
 * The results used in the test can be recreated by using the following R-code:
 * 
 * vec <- c(1,2,3,4)
 * 
 * for (i in vec) {
 *    result <- 0
 *    for (j in vec) {
 *       result <- result + ppois(i,j+0.5)
 *    }
 *    print(result/length(vec))
 * }
*/
TEST(MatrixNormalization, PoissonKDE_one_row)
{
  
   std::vector<int> values;
   std::vector<double> r_results;
   
   for(int i=1; i<5; ++i)
     values.push_back(i);
   
   r_results.push_back(0.2605277);
   r_results.push_back(0.4617713);  
   r_results.push_back(0.6427156);
   r_results.push_back(0.7825377);
   
   
   PoissonEstimator poisson;
   
   for(int j=0; j<4; ++j) {
      double norm = poisson.normalizeValue(values.begin()+j,values.begin(), values.end());
      EXPECT_NEAR(norm, r_results[j], TOLERANCE);
   } 
}



/**
 * The results used in the test can be recreated by using the following R-code:
 * 
 * vec <- c(1,1,1,1,1)
 * 
 * for (i in vec) {
 *    result <- 0
 *    for (j in vec) {
 *       result <- result + ppois(i,j+0.5)
 *    }
 *    print(result/length(vec))
 * }
*/
TEST(MatrixNormalization, PoissonKDE_all_one)
{
  
   std::vector<int> values;
   
   for(int i=0; i<5; ++i)
     values.push_back(1);  
   
   PoissonEstimator poisson;
   
   for(int j=0; j<5; ++j) {
      double norm = poisson.normalizeValue(values.begin()+j,values.begin(), values.end());
      EXPECT_NEAR(norm, 0.5578254, TOLERANCE);
   } 
}

/**
 * The results used in the test can be recreated by using the following R-code:
 * 
 * vec <- c(1,2,3,4)
 * 
 * for (i in vec) {
 *   result <- 0
 *   for (j in vec) {
 *     result <- result + pnorm(((i-j)*4)/sd(vec))
 *   }
 *   print(result/length(vec))
 * }
*/
TEST(MatrixNormalization, GaussKDE_ints)
{
  
   std::vector<int> values;
   std::vector<double> r_results;
   
   for(int i=1; i<5; ++i)
     values.push_back(i);
   
   r_results.push_back(0.1252432);
   r_results.push_back(0.375);
   r_results.push_back(0.625);
   r_results.push_back(0.8747568); 
   
   
   GaussEstimator gauss;
   
   for(int j=0; j<4; ++j) {
      double norm = gauss.normalizeValue(values.begin()+j,values.begin(), values.end());
      EXPECT_NEAR(norm, r_results[j], TOLERANCE);
   } 
}

/**
 * The results used in the test can be recreated by using the following R-code:
 * 
 * vec <- c(1.453,0.4243,3.42111,0.5332,1.111)
 * 
 * for (i in vec) {
 *   result <- 0
 *   for (j in vec) {
 *     result <- result + pnorm(((i-j)*4)/sd(vec))
 *   }
 *   print(result/length(vec))
 * }
*/
TEST(MatrixNormalization, GaussKDE_doubles)
{
  
   std::vector<double> values;
   std::vector<double> r_results;
   
   
   values.push_back(1.453);
   values.push_back(0.4243);
   values.push_back(3.42111);
   values.push_back(0.5332);
   values.push_back(1.111);
   
   r_results.push_back(0.6737976);
   r_results.push_back(0.1743342);
   r_results.push_back(0.9);
   r_results.push_back(0.2339638);
   r_results.push_back(0.5179044); 
   
   
   GaussEstimator gauss;
   
   for(int j=0; j<5; ++j) {
      double norm = gauss.normalizeValue(values.begin()+j,values.begin(), values.end());
      EXPECT_NEAR(norm, r_results[j], TOLERANCE);
   } 
}

/**
 * The results used in the test can be recreated by using the following R-code:
 * 
 * If standard-deviation is 0, we use pnorm(0) instead of pnorm(((i-j)*4)/sd(vec)), as we cannot divide by 0
 * 
 * vec <- c(1,1,1,1,1)
 * 
 * for (i in vec) {
 *   result <- 0
 *   for (j in vec) {
 *     result <- result + pnorm(0)  
 *   }
 *   print(result/length(vec))
 * }
*/
TEST(MatrixNormalization, GaussKDE_all_one)
{
   std::vector<int> values;
   
   for(int i=0; i<5; ++i)
     values.push_back(1);
     
   GaussEstimator gauss;
   
   for(int j=0; j<5; ++j) {
      double norm = gauss.normalizeValue(values.begin()+j,values.begin(), values.end());
      EXPECT_NEAR(norm, 0.5, TOLERANCE);
   } 
}

/**
 * The results used in the test can be recreated by using the following R-code:
 * 
 * vec <- c(-15.0,-0.0564,3.42111,8.0,-1.111,0.0)
 * 
 * for (i in vec) {
 *   result <- 0
 *   for (j in vec) {
 *     result <- result + pnorm(((i-j)*4)/sd(vec))
 *   }
 *   print(result/length(vec))
 * }
*/
TEST(MatrixNormalization, GaussKDE_negative_values)
{
  
   std::vector<double> values;
   std::vector<double> r_results;
   
   
   values.push_back(-15.0);
   values.push_back(-0.0564);
   values.push_back(3.42111);
   values.push_back(8.0);
   values.push_back(-1.111);
   values.push_back(0.0);
   
   r_results.push_back(0.08333333);
   r_results.push_back(0.4552912);
   r_results.push_back(0.7375622);
   r_results.push_back(0.9151856);
   r_results.push_back(0.3473984);
   r_results.push_back(0.4612292); 
   
   
   GaussEstimator gauss;
   
   for(int j=0; j<5; ++j) {
      double norm = gauss.normalizeValue(values.begin()+j,values.begin(), values.end());
      EXPECT_NEAR(norm, r_results[j], TOLERANCE);
   } 
}

TEST(MatrixNormalization, ExclusiveZ_doubles)
{
  
   std::vector<double> values;
   std::vector<double> r_results;
   
   
   values.push_back(1.453);
   values.push_back(0.4243);
   values.push_back(3.42111);
   values.push_back(0.5332);
   values.push_back(1.111);
   
   r_results.push_back(exclusiveZ(1.453, values.begin(), values.end()));
   r_results.push_back(exclusiveZ(0.4243, values.begin(), values.end()));
   r_results.push_back(exclusiveZ(3.42111, values.begin(), values.end()));
   r_results.push_back(exclusiveZ(0.5332, values.begin(), values.end()));
   r_results.push_back(exclusiveZ(1.111, values.begin(), values.end()));

   ExclusiveZScore zscore;
   
   for(int j=0; j<5; ++j) {
      double norm = zscore.normalizeValue(values.begin()+j,values.begin(), values.end());
      EXPECT_NEAR(norm, r_results[j], TOLERANCE);
   } 
}


TEST(MatrixNormalization, ExclusiveZ_negative_values)
{
  
   std::vector<double> values;
   std::vector<double> r_results;
   
   
   values.push_back(-15.0);
   values.push_back(-0.0564);
   values.push_back(3.42111);
   values.push_back(8.0);
   values.push_back(-1.111);
   values.push_back(0.0);
   
   r_results.push_back(exclusiveZ(-15.0, values.begin(), values.end()));
   r_results.push_back(exclusiveZ(-0.0564, values.begin(), values.end()));
   r_results.push_back(exclusiveZ(3.42111, values.begin(), values.end()));
   r_results.push_back(exclusiveZ(8.0, values.begin(), values.end()));
   r_results.push_back(exclusiveZ(-1.111, values.begin(), values.end()));
   r_results.push_back(exclusiveZ(0.0, values.begin(), values.end()));

   ExclusiveZScore zscore;
   
   for(int j=0; j<5; ++j) {
      double norm = zscore.normalizeValue(values.begin()+j,values.begin(), values.end());
      EXPECT_NEAR(norm, r_results[j], TOLERANCE);
   } 
}

TEST(MatrixNormalization, ExclusiveZ_zero_sd)
{
  
   std::vector<double> values;
   std::vector<double> r_results;
   
   
   values.push_back(1.5);
   values.push_back(1.5);
   values.push_back(1.5);
   values.push_back(1.5);
   values.push_back(2.5);
   
   r_results.push_back(exclusiveZ(1.5, values.begin(), values.end()));
   r_results.push_back(exclusiveZ(1.5, values.begin(), values.end()));
   r_results.push_back(exclusiveZ(1.5, values.begin(), values.end()));
   r_results.push_back(exclusiveZ(1.5, values.begin(), values.end()));
   r_results.push_back(exclusiveZ(2.5, values.begin(), values.end()));

   ExclusiveZScore zscore;
   
   for(int j=0; j<5; ++j) {
      double norm = zscore.normalizeValue(values.begin()+j,values.begin(), values.end());
      EXPECT_NEAR(norm, r_results[j], TOLERANCE);
   } 
}


/**
 * The results used in the test can be recreated by using the following R-code:
 * 
 * Gene1 <- c(1,2,3)
 * Gene2 <- c(4,5,6)
 * Gene3 <- c(7,8,9)
 * Gene4 <- c(10,11,12)
 * Gene5 <- c(13,14,15)
 * 
 * valueMatrix <- rbind(Gene1,Gene2,Gene3,Gene4,Gene5)
 * resultMatrix <- matrix(nrow=nrow(valueMatrix),ncol=ncol(valueMatrix))
 * 
 * for (i in 1:nrow(valueMatrix)) {
 *   for (j in 1:ncol(valueMatrix)) {
 *     result <- 0
 *     for (k in valueMatrix[i, ]) {
 *       result <- result + ppois(valueMatrix[i,j],k+0.5)
 *    }
 *    resultMatrix[i,j] <- result/ncol(valueMatrix)
 *  }
 * }
 * print(resultMatrix)
 * 
*/
TEST(MatrixNormalization, Poisson_simple_matrix)
{
	MatrixNormalization norm;
	PoissonEstimator poisson;
  
	const unsigned int num_rows = 5;
	const unsigned int num_cols = 3;

	DenseMatrix mat(num_rows, num_cols);

	// Fill the matrix
	mat(0, 0) = 1;
	mat(0, 1) = 2;
	mat(0, 2) = 3;
	mat(1, 0) = 4;
	mat(1, 1) = 5;
	mat(1, 2) = 6;
	mat(2, 0) = 7;
	mat(2, 1) = 8;
	mat(2, 2) = 9;
	mat(3, 0) = 10;
	mat(3, 1) = 11;
	mat(3, 2) = 12;
	mat(4, 0) = 13;
	mat(4, 1) = 14;
	mat(4, 2) = 15;
	
	DenseMatrix results(num_rows, num_cols);
	
	results(0, 0) = 0.3270037;
	results(0, 1) = 0.5578357;
	results(0, 2) = 0.7428554;
	results(1, 0) = 0.3710978;
	results(1, 1) = 0.5336299;
	results(1, 2) = 0.6812034;
	results(2, 0) = 0.3929663;
	results(2, 1) = 0.5256319;
	results(2, 2) = 0.6504024;
	results(3, 0) = 0.4065142;
	results(3, 1) = 0.5214280;
	results(3, 2) = 0.6312954;
	results(4, 0) = 0.4159518;
	results(4, 1) = 0.5187583;
	results(4, 2) = 0.6179891;
		
	DenseMatrix normalized = norm.normalizeMatrix(mat,poisson); //Matrix for results using "normalizeMatrix"
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized(i,j), results(i,j), TOLERANCE);
	  }
	}
	
	//Testing "normalizeRange" on some rows of mat
	DenseMatrix normalized2(num_rows, num_cols); // Matrix for results using "normalizeRange" on each row
	
	RowMajorMatrixIterator<Matrix> row_it(&mat, 0);
	norm.normalizeRange(row_it->begin(), row_it->end(), normalized2, 0, poisson);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(0,j), results(0,j), TOLERANCE);
	  }
	
	RowMajorMatrixIterator<Matrix> row_it2(&mat, 2);
	norm.normalizeRange(row_it2->begin(), row_it2->end(), normalized2, 2, poisson);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(2,j), results(2,j), TOLERANCE);
	  }
  

  
}

/**
 * The results used in the test can be recreated by using the following R-code:
 * 
 * Gene1 <- c(0,0,70,1)
 * Gene2 <- c(2,3,4,11)
 * Gene3 <- c(0,0,0,0)
 * Gene4 <- c(12,1,87,91)
 * 
 * valueMatrix <- rbind(Gene1,Gene2,Gene3,Gene4)
 * resultMatrix <- matrix(nrow=nrow(valueMatrix),ncol=ncol(valueMatrix))
 * 
 * for (i in 1:nrow(valueMatrix)) {
 *   for (j in 1:ncol(valueMatrix)) {
 *     result <- 0
 *     for (k in valueMatrix[i, ]) {
 *       result <- result + ppois(valueMatrix[i,j],k+0.5)
 *    }
 *    resultMatrix[i,j] <- result/ncol(valueMatrix)
 *  }
 * }
 * print(resultMatrix)
 * 
*/
TEST(MatrixNormalization, Poisson_matrix_with_zero)
{
	MatrixNormalization norm;
	PoissonEstimator poisson;
  
	const unsigned int num_rows = 4;
	const unsigned int num_cols = 4;

	DenseMatrix mat(num_rows, num_cols);

	// Fill the matrix
	mat(0, 0) = 0;
	mat(0, 1) = 0;
	mat(0, 2) = 70;
	mat(0, 3) = 1;
	mat(1, 0) = 2;
	mat(1, 1) = 3;
	mat(1, 2) = 4;
	mat(1, 3) = 111;
	mat(2, 0) = 0;
	mat(2, 1) = 0;
	mat(2, 2) = 0;
	mat(2, 3) = 0;
	mat(3, 0) = 12;
	mat(3, 1) = 1;
	mat(3, 2) = 87;
	mat(3, 3) = 91;

	
	DenseMatrix results(num_rows, num_cols);
	
	results(0, 0) = 0.3590479;
	results(0, 1) = 0.3590479 ;
	results(0, 2) = 0.8769829;
	results(0, 3) = 0.5943543;
	results(1, 0) = 0.2595596;
	results(1, 1) = 0.4091262;
	results(1, 2) = 0.5371816;
	results(1, 3) = 0.8765758;
	results(2, 0) = 0.6065307;
	results(2, 1) = 0.6065307;
	results(2, 2) = 0.6065307;
	results(2, 3) = 0.6065307;
	results(3, 0) = 0.3797438;
	results(3, 1) = 0.1394689;
	results(3, 2) = 0.7125841;
	results(3, 3) = 0.7944403;
	
	DenseMatrix normalized = norm.normalizeMatrix(mat,poisson);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized(i,j), results(i,j), TOLERANCE);
	  }
	}

	
	//Testing "normalizeRange" on some rows of mat
	DenseMatrix normalized2(num_rows, num_cols); // Matrix for results using "normalizeRange" on each num_rows
	
	RowMajorMatrixIterator<Matrix> row_it(&mat, 0);
	norm.normalizeRange(row_it->begin(), row_it->end(), normalized2, 0, poisson);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(0,j), results(0,j), TOLERANCE);
	}
	
	RowMajorMatrixIterator<Matrix> row_it2(&mat, 2);
	norm.normalizeRange(row_it2->begin(), row_it2->end(), normalized2, 2, poisson);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(2,j), results(2,j), TOLERANCE);
	}
  
}

/**
 * The results used in the test can be recreated by using the following R-code:
 * 
 * Gene1 <- c(7,7)
 * Gene2 <- c(7,7)
 * Gene3 <- c(7,7)
 * 
 * valueMatrix <- rbind(Gene1,Gene2,Gene3)
 * resultMatrix <- matrix(nrow=nrow(valueMatrix),ncol=ncol(valueMatrix))
 * 
 * for (i in 1:nrow(valueMatrix)) {
 *   for (j in 1:ncol(valueMatrix)) {
 *     result <- 0
 *     for (k in valueMatrix[i, ]) {
 *       result <- result + ppois(valueMatrix[i,j],k+0.5)
 *    }
 *    resultMatrix[i,j] <- result/ncol(valueMatrix)
 *  }
 * }
 * print(resultMatrix)
 * 
*/
TEST(MatrixNormalization, Poisson_matrix_all_seven)
{
	MatrixNormalization norm;
	PoissonEstimator poisson;
  
	const unsigned int num_rows = 3;
	const unsigned int num_cols = 2;

	DenseMatrix mat(num_rows, num_cols);

	// Fill the matrix
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    mat(i,j) = 7;
	  }
	}
	
	DenseMatrix normalized = norm.normalizeMatrix(mat,poisson);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized(i,j), 0.5246385, TOLERANCE);
	  }
	}
	
	//Testing "normalizeRange" on some rows of mat (all rows should be the same here)
	DenseMatrix normalized2(num_rows, num_cols); // Matrix for results using "normalizeRange" on each num_rows
	
	RowMajorMatrixIterator<Matrix> row_it(&mat, 0);
	norm.normalizeRange(row_it->begin(), row_it->end(), normalized2, 0, poisson);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(0,j), 0.5246385, TOLERANCE);
	}
	
	RowMajorMatrixIterator<Matrix> row_it2(&mat, 2);
	norm.normalizeRange(row_it2->begin(), row_it2->end(), normalized2, 2, poisson);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(2,j), 0.5246385, TOLERANCE);
	}

  
}

/**
 * The results used in the test can be recreated by using the following R-code:
 * 
 * Gene1 <- c(1)
 * Gene2 <- c(0)
 * Gene3 <- c(70)
 * Gene4 <- c(22)
 * Gene5 <- c(12)
 * 
 * valueMatrix <- rbind(Gene1,Gene2,Gene3,Gene4, Gene5)
 * resultMatrix <- matrix(nrow=nrow(valueMatrix),ncol=ncol(valueMatrix))
 * 
 * for (i in 1:nrow(valueMatrix)) {
 *   for (j in 1:ncol(valueMatrix)) {
 *     result <- 0
 *     for (k in valueMatrix[i, ]) {
 *       result <- result + ppois(valueMatrix[i,j],k+0.5)
 *    }
 *    resultMatrix[i,j] <- result/ncol(valueMatrix)
 *  }
 * }
 * print(resultMatrix)
 * 
*/
TEST(MatrixNormalization, Poisson_matrix_single_sample)
{
	MatrixNormalization norm;
	PoissonEstimator poisson;
  
	const unsigned int num_rows = 5;
	const unsigned int num_cols = 1;

	DenseMatrix mat(num_rows, num_cols);

	// Fill the matrix
	mat(0, 0) = 1;
	mat(1, 0) = 0;
	mat(2, 0) = 70;
	mat(3, 0) = 22;
	mat(4, 0) = 12;
	
	DenseMatrix results(num_rows, num_cols);
	
	results(0, 0) = 0.5578254;
	results(1, 0) = 0.6065307;
	results(2, 0) = 0.5079317;
	results(3, 0) = 0.5140878;
	results(4, 0) = 0.5189752;
	
	DenseMatrix normalized = norm.normalizeMatrix(mat,poisson);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized(i,j), results(i,j), TOLERANCE);
	  }
	}
  
  
  	//Testing "normalizeRange" on some rows of mat (all rows should be the same here)
	DenseMatrix normalized2(num_rows, num_cols); // Matrix for results using "normalizeRange" on each num_rows
	
	RowMajorMatrixIterator<Matrix> row_it(&mat, 0);
	norm.normalizeRange(row_it->begin(), row_it->end(), normalized2, 0, poisson);
	EXPECT_NEAR(normalized2(0,0), results(0,0), TOLERANCE);
	
	RowMajorMatrixIterator<Matrix> row_it2(&mat, 2);
	norm.normalizeRange(row_it2->begin(), row_it2->end(), normalized2, 2, poisson);
	EXPECT_NEAR(normalized2(2,0), results(2,0), TOLERANCE);

}

/**
 * The results used in the test can be recreated by using the following R-code:
 * 
 * Gene1 <- c(1,2,3)
 * Gene2 <- c(4,5,6)
 * Gene3 <- c(7,8,9)
 * Gene4 <- c(10,11,12)
 * Gene5 <- c(13,14,15)
 * 
 * valueMatrix <- rbind(Gene1,Gene2,Gene3,Gene4,Gene5)
 * resultMatrix <- matrix(nrow=nrow(valueMatrix),ncol=ncol(valueMatrix))
 * 
 * for (i in 1:nrow(valueMatrix)) {
 *   for (j in 1:ncol(valueMatrix)) {
 *     result <- 0
 *     for (k in valueMatrix[i, ]) {
 *       result <- result + pnorm(((valueMatrix[i,j]-k)*4)/sd(valueMatrix[i, ]))
 *     }
 *     resultMatrix[i,j] <- result/ncol(valueMatrix)
 *   }
 * }
 * 
 * print(resultMatrix)
 * 
*/
TEST(MatrixNormalization, Gauss_simple_matrix)
{
	MatrixNormalization norm;
	GaussEstimator gauss;
  
	const unsigned int num_rows = 5;
	const unsigned int num_cols = 3;

	DenseMatrix mat(num_rows, num_cols);

	// Fill the matrix
	mat(0, 0) = 1;
	mat(0, 1) = 2;
	mat(0, 2) = 3;
	mat(1, 0) = 4;
	mat(1, 1) = 5;
	mat(1, 2) = 6;
	mat(2, 0) = 7;
	mat(2, 1) = 8;
	mat(2, 2) = 9;
	mat(3, 0) = 10;
	mat(3, 1) = 11;
	mat(3, 2) = 12;
	mat(4, 0) = 13;
	mat(4, 1) = 14;
	mat(4, 2) = 15;
	
	DenseMatrix results(num_rows, num_cols);
	
	results(0, 0) = 0.1666772;
	results(0, 1) = 0.5;
	results(0, 2) = 0.8333228;
	results(1, 0) = 0.1666772;
	results(1, 1) = 0.5;
	results(1, 2) = 0.8333228;
	results(2, 0) = 0.1666772;
	results(2, 1) = 0.5;
	results(2, 2) = 0.8333228;
	results(3, 0) = 0.1666772;
	results(3, 1) = 0.5;
	results(3, 2) = 0.8333228;
	results(4, 0) = 0.1666772;
	results(4, 1) = 0.5;
	results(4, 2) = 0.8333228;
		
	DenseMatrix normalized = norm.normalizeMatrix(mat,gauss);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized(i,j), results(i,j), TOLERANCE);
	  }
	}
	
	//Testing "normalizeRange" on some rows of mat
	DenseMatrix normalized2(num_rows, num_cols); // Matrix for results using "normalizeRange" on each num_rows
	
	RowMajorMatrixIterator<Matrix> row_it(&mat, 0);
	norm.normalizeRange(row_it->begin(), row_it->end(), normalized2, 0, gauss);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(0,j), results(0,j), TOLERANCE);
	}
	
	RowMajorMatrixIterator<Matrix> row_it2(&mat, 2);
	norm.normalizeRange(row_it2->begin(), row_it2->end(), normalized2, 2, gauss);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(2,j), results(2,j), TOLERANCE);
	}
}

/**
 * The results used in the test can be recreated by using the following R-code:
 * 
 * Gene1 <- c(0.6544,1.533,0.0,17.54)
 * Gene2 <- c(-25.5664,3.65,-5.323,11.223)
 * Gene3 <- c(15.6,16.5,17.99,5.3)
 * Gene4 <- c(-12.5,0.0,-87.4,-26.9876)
 * 
 * valueMatrix <- rbind(Gene1,Gene2,Gene3,Gene4)
 * resultMatrix <- matrix(nrow=nrow(valueMatrix),ncol=ncol(valueMatrix))
 * 
 * for (i in 1:nrow(valueMatrix)) {
 *   for (j in 1:ncol(valueMatrix)) {
 *     result <- 0
 *     for (k in valueMatrix[i, ]) {
 *       result <- result + pnorm(((valueMatrix[i,j]-k)*4)/sd(valueMatrix[i, ]))
 *     }
 *     resultMatrix[i,j] <- result/ncol(valueMatrix)
 *   }
 * }
 * 
 * print(resultMatrix)
 * 
*/
TEST(MatrixNormalization, Gauss_double_matrix)
{
	MatrixNormalization norm;
	GaussEstimator gauss;
  
	const unsigned int num_rows = 4;
	const unsigned int num_cols = 4;

	DenseMatrix mat(num_rows, num_cols);

	// Fill the matrix
	mat(0, 0) = 0.6544;
	mat(0, 1) = 1.533;
	mat(0, 2) = 0.0;
	mat(0, 3) = 17.54;
	mat(1, 0) = -25.5664;
	mat(1, 1) = 3.65;
	mat(1, 2) = -5.323;
	mat(1, 3) = 11.223;
	mat(2, 0) = 15.6;
	mat(2, 1) = 16.5;
	mat(2, 2) = 17.99;
	mat(2, 3) = 5.3;
	mat(3, 0) = -12.5;
	mat(3, 1) = 0.0;
	mat(3, 2) = -87.4;
	mat(3, 3) = -26.9876;

	
	DenseMatrix results(num_rows, num_cols);
	
	results(0, 0) = 0.3650719;
	results(0, 1) = 0.4820455;
	results(0, 2) = 0.2778826;
	results(0, 3) = 0.8750000;
	results(1, 0) = 0.1250000;
	results(1, 1) = 0.6290852;
	results(1, 2) = 0.3779869;
	results(1, 3) = 0.8679278;
	results(2, 0) = 0.4539877;
	results(2, 1) = 0.5961401;
	results(2, 2) = 0.8248722;
	results(2, 3) = 0.1250000;
	results(3, 0) = 0.6327640;
	results(3, 1) = 0.8497595;
	results(3, 2) = 0.1250000;
	results(3, 3) = 0.3924765;
	
	DenseMatrix normalized = norm.normalizeMatrix(mat,gauss);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized(i,j), results(i,j), TOLERANCE);
	  }
	}
	
	//Testing "normalizeRange" on some rows of mat
	DenseMatrix normalized2(num_rows, num_cols); // Matrix for results using "normalizeRange" on each num_rows
	
	RowMajorMatrixIterator<Matrix> row_it(&mat, 0);
	norm.normalizeRange(row_it->begin(), row_it->end(), normalized2, 0, gauss);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(0,j), results(0,j), TOLERANCE);
	}
	
	RowMajorMatrixIterator<Matrix> row_it2(&mat, 3);
	norm.normalizeRange(row_it2->begin(), row_it2->end(), normalized2, 3, gauss);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(3,j), results(3,j), TOLERANCE);
	}

  
}

TEST(MatrixNormalization, ZScore_simple_matrix)
{
	MatrixNormalization norm;
	ExclusiveZScore zscore;
  
	const unsigned int num_rows = 5;
	const unsigned int num_cols = 3;

	DenseMatrix mat(num_rows, num_cols);

	// Fill the matrix
	mat(0, 0) = 1;
	mat(0, 1) = 2;
	mat(0, 2) = 3;
	mat(1, 0) = 4;
	mat(1, 1) = 5;
	mat(1, 2) = 6;
	mat(2, 0) = 7;
	mat(2, 1) = 8;
	mat(2, 2) = 9;
	mat(3, 0) = 10;
	mat(3, 1) = 11;
	mat(3, 2) = 12;
	mat(4, 0) = 13;
	mat(4, 1) = 14;
	mat(4, 2) = 15;
	
	RowMajorMatrixIterator<Matrix> row_it(&mat, 0);
	DenseMatrix results(num_rows, num_cols);
	
	for(unsigned int i=0; i<num_rows;++i) {
	 for(unsigned int j=0; j<num_cols;++j) {
	   results.set(i,j,exclusiveZ(mat(i,j),row_it->begin(), row_it->end()));
	 }
	 ++row_it;
	}
	
	DenseMatrix normalized = norm.normalizeMatrix(mat,zscore);
	
	for(unsigned int i=0; i<num_rows;++i) {
	  for(unsigned int j=0; j<num_cols;++j) {
	    EXPECT_NEAR(normalized(i,j), results(i,j), TOLERANCE);
	  }
	}

	
	//Testing "normalizeRange" on some rows of mat
	DenseMatrix normalized2(num_rows, num_cols); // Matrix for results using "normalizeRange" on each num_rows
	
	RowMajorMatrixIterator<Matrix> row_it2(&mat, 0);
	norm.normalizeRange(row_it2->begin(), row_it2->end(), normalized2, 0, zscore);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(0,j), results(0,j), TOLERANCE);
	}
	
	RowMajorMatrixIterator<Matrix> row_it3(&mat, 2);
	norm.normalizeRange(row_it3->begin(), row_it3->end(), normalized2, 2, zscore);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(2,j), results(2,j), TOLERANCE);
	}
}

TEST(MatrixNormalization, ZScore_double_matrix)
{
	MatrixNormalization norm;
	ExclusiveZScore zscore;
  
	const unsigned int num_rows = 4;
	const unsigned int num_cols = 4;

	DenseMatrix mat(num_rows, num_cols);

	// Fill the matrix
	mat(0, 0) = 0.6544;
	mat(0, 1) = 1.533;
	mat(0, 2) = 0.0;
	mat(0, 3) = 17.54;
	mat(1, 0) = -25.5664;
	mat(1, 1) = 3.65;
	mat(1, 2) = -5.323;
	mat(1, 3) = 11.223;
	mat(2, 0) = 15.6;
	mat(2, 1) = 16.5;
	mat(2, 2) = 17.99;
	mat(2, 3) = 5.3;
	mat(3, 0) = -12.5;
	mat(3, 1) = 0.0;
	mat(3, 2) = -87.4;
	mat(3, 3) = -26.9876;
	
	RowMajorMatrixIterator<Matrix> row_it(&mat, 0);
	DenseMatrix results(num_rows, num_cols);
	
	for(unsigned int i=0; i<num_rows;++i) {
	 for(unsigned int j=0; j<num_cols;++j) {
	   results.set(i,j,exclusiveZ(mat(i,j),row_it->begin(), row_it->end()));
	 }
	 ++row_it;
	}
	
	DenseMatrix normalized = norm.normalizeMatrix(mat,zscore);
	
	for(unsigned int i=0; i<num_rows;++i) {
	  for(unsigned int j=0; j<num_cols;++j) {
	    EXPECT_NEAR(normalized(i,j), results(i,j), TOLERANCE);
	  }
	}

	//Testing "normalizeRange" on some rows of mat
	DenseMatrix normalized2(num_rows, num_cols); // Matrix for results using "normalizeRange" on each num_rows
	
	RowMajorMatrixIterator<Matrix> row_it2(&mat, 0);
	norm.normalizeRange(row_it2->begin(), row_it2->end(), normalized2, 0, zscore);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(0,j), results(0,j), TOLERANCE);
	}
	
	RowMajorMatrixIterator<Matrix> row_it3(&mat, 3);
	norm.normalizeRange(row_it3->begin(), row_it3->end(), normalized2, 3, zscore);
	
	for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_NEAR(normalized2(3,j), results(3,j), TOLERANCE);
	}

  
}