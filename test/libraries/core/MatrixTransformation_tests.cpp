#include <gtest/gtest.h>
#include <config.h>
#include <genetrail2/core/MatrixTransformation.h>
#include <genetrail2/core/DenseMatrix.h>

using namespace GeneTrail; 

void fillTestMatrix(Matrix& mat)
{
	mat(0,0) =  1.0;
	mat(0,1) =  2.0;
	mat(0,2) =  3.0;
	mat(1,0) =  4.0;
	mat(1,1) =  5.0;
	mat(1,2) =  6.0;
	mat(2,0) =  7.0;
	mat(2,1) =  8.0;
	mat(2,2) =  9.0;
	mat(3,0) = 10.0;
	mat(3,1) = 11.0;
	mat(3,2) = 12.0;
	mat(4,0) = 13.0;
	mat(4,1) = 14.0;
	mat(4,2) = 15.0;
}

void fillTestMatrixNegative(Matrix& mat)
{
	mat(0,0) =  -1.0;
	mat(0,1) =  2.0;
	mat(0,2) =  3.0;
	mat(1,0) =  -4.0;
	mat(1,1) =  5.0;
	mat(1,2) =  -6.0;
	mat(2,0) =  7.0;
	mat(2,1) =  8.0;
	mat(2,2) =  9.0;
	mat(3,0) = -10.0;
	mat(3,1) = 11.0;
	mat(3,2) = -12.0;
	mat(4,0) = -13.0;
	mat(4,1) = -14.0;
	mat(4,2) = 15.0;
}

TEST(MT, valuesToRanks1)
{
	const unsigned int num_rows = 5;
	const unsigned int num_cols = 3;

	DenseMatrix mat(num_rows, num_cols);

	// Fill the matrix
	mat(0, 0) = 180;
	mat(1, 0) = 7;
	mat(2, 0) = 43;
	mat(3, 0) = -43;
	mat(4, 0) = 0;
	
	mat(0, 1) = -6;
	mat(1, 1) = 77;
	mat(2, 1) = 99998;
	mat(3, 1) = 89;
	mat(4, 1) = 179;

	mat(0, 2) = 11;
	mat(1, 2) = 122;
	mat(2, 2) = 133;
	mat(3, 2) = 144;
	mat(4, 2) = 155;
	
	DenseMatrix results(num_rows, num_cols);
	
// 	results(0, 0) = 1;
// 	results(1, 0) = 3;
// 	results(2, 0) = 2;
// 	results(3, 0) = 5;
// 	results(4, 0) = 4;
// 	
// 	results(0, 1) = 3;
// 	results(1, 1) = 5;
// 	results(2, 1) = 4;
// 	results(3, 1) = 2;
// 	results(4, 1) = 1;
// 	results(0 ,2) = 5;
// 	results(1, 2) = 4;
// 	results(2, 2) = 3;
// 	results(3, 2) = 2;
// 	results(4, 2) = 1;
	
	results(0, 0) = 1;
	results(1, 0) = 3;
	results(2, 0) = 2;
	results(3, 0) = 5;
	results(4, 0) = 4;

	results(0, 1) = 5;
	results(1, 1) = 4;
	results(2, 1) = 1;
	results(3, 1) = 3;
	results(4, 1) = 2;
	
	results(0 ,2) = 5;
	results(1, 2) = 4;
	results(2, 2) = 3;
	results(3, 2) = 2;
	results(4, 2) = 1;
	
	valuesToRanks(mat);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_EQ(mat(i,j), results(i,j));
	  }
	}
}

TEST(MT, valuesToRanks2)
{
	const unsigned int num_rows = 4;
	const unsigned int num_cols = 4;

	DenseMatrix mat(num_rows, num_cols);

	// Fill the matrix
	mat(0, 0) = 18.5;
	mat(1, 0) = 7.009;
	mat(2, 0) = -43.7;
	mat(3, 0) = -43.6;
	
	mat(0, 1) = 0.0;
	mat(1, 1) = 12.2;
	mat(2, 1) = -17.7654;
	mat(3, 1) = -89.13;
	
	mat(0, 2) = 11111;
	mat(1, 2) = 2.1;
	mat(2, 2) = 0.0009;
	mat(3, 2) = -12.32;
	
	mat(0, 3) = 11.11;
	mat(1, 3) = 998.0;
	mat(2, 3) = 133.7;
	mat(3, 3) = 144.1;

	
	DenseMatrix results(num_rows, num_cols);
	
	results(0, 0) = 1;
	results(1, 0) = 2;
	results(2, 0) = 4;
	results(3, 0) = 3;
	
	results(0, 1) = 2;
	results(1, 1) = 1;
	results(2, 1) = 3;
	results(3, 1) = 4;
	
	results(0 ,2) = 1;
	results(1, 2) = 2;
	results(2, 2) = 3;
	results(3, 2) = 4;
	
	results(0 ,3) = 4;
	results(1, 3) = 1;
	results(2, 3) = 3;
	results(3, 3) = 2;
	
	valuesToRanks(mat);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_EQ(mat(i,j), results(i,j));
	  }
	}
}

TEST(MT, valuesToRanksEasy)
{
	const unsigned int num_rows = 2;
	const unsigned int num_cols = 2;

	DenseMatrix mat(num_rows, num_cols);

	// Fill the matrix
	mat(0, 0) = 18.5;
	mat(1, 0) = 7.009;
	mat(0, 1) = 0.0;
	mat(1, 1) = 12.2;
	
	DenseMatrix results(num_rows, num_cols);
	
	results(0, 0) = 1;
	results(1, 0) = 2;
	results(0, 1) = 2;
	results(1, 1) = 1;
	
	
	valuesToRanks(mat);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_EQ(mat(i,j), results(i,j));
	  }
	}

  
  
} 

TEST(MT, valuesToRanksInt)
{
  
	const unsigned int num_rows = 3;
	const unsigned int num_cols = 2;

	DenseMatrix mat(num_rows, num_cols);

	// Fill the matrix
	mat(0, 0) = 1;
	mat(1, 0) = 7;
	mat(2, 0) = 5;
	mat(0, 1) = 2;
	mat(1, 1) = 17;
	mat(2, 1) = 7;
	
	DenseMatrix results(num_rows, num_cols);
	
	results(0, 0) = 3;
	results(1, 0) = 1;
	results(2, 0) = 2;
	results(0, 1) = 3;
	results(1, 1) = 1;
	results(2, 1) = 2;
	
	
	valuesToRanks(mat);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_EQ(mat(i,j), results(i,j));
	  }
	}
}

TEST(MT, upweightEnds) 
{
	const unsigned int num_rows = 4;
	const unsigned int num_cols = 3;
	
	DenseMatrix mat(num_rows, num_cols);
	DenseMatrix mat2(num_rows, num_cols);
	DenseMatrix transformedMat(num_rows, num_cols);
	DenseMatrix transformedMat2(num_rows, num_cols);
	fillTestMatrix(mat);
	fillTestMatrix(transformedMat);
	fillTestMatrixNegative(mat2);
	fillTestMatrixNegative(transformedMat2);
	upweightEnds(transformedMat);
	upweightEnds(transformedMat2);
	
	for(size_t i=0; i<num_rows; ++i){
	  for(size_t j=0; j< num_cols; ++j){
	    mat.set(i,j,std::abs(num_rows/2 +0.5 -mat(i,j)));
	    mat2.set(i,j,std::abs(num_rows/2 + 0.5 -mat2(i,j)));
	  }
	}
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_EQ(transformedMat(i,j),mat(i,j));
	    EXPECT_EQ(transformedMat2(i,j),mat2(i,j));
	  }
	}
}

TEST(MT, upweightTail) 
{
	const unsigned int num_rows = 4;
	const unsigned int num_cols = 3;
	
	DenseMatrix mat(num_rows, num_cols);
	DenseMatrix mat2(num_rows, num_cols);
	DenseMatrix transformedMat(num_rows, num_cols);
	DenseMatrix transformedMat2(num_rows, num_cols);
	fillTestMatrix(mat);
	fillTestMatrix(transformedMat);
	fillTestMatrixNegative(mat2);
	fillTestMatrixNegative(transformedMat2);
	upweightTail(transformedMat);
	upweightTail(transformedMat2);
	
	for(size_t i=0; i<num_rows; ++i){
	  for(size_t j=0; j< num_cols; ++j){
	    mat.set(i,j,num_rows/2 + 0.5 -mat(i,j));
	    mat2.set(i,j,num_rows/2 + 0.5 -mat2(i,j));
	  }
	}
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_EQ(transformedMat(i,j),mat(i,j));
	    EXPECT_EQ(transformedMat2(i,j),mat2(i,j));
	  }
	}
}

TEST(MT, abs) 
{
	const unsigned int num_rows = 4;
	const unsigned int num_cols = 3;
	
	DenseMatrix mat(num_rows, num_cols);
	DenseMatrix mat2(num_rows, num_cols);
	DenseMatrix transformedMat(num_rows, num_cols);
	DenseMatrix transformedMat2(num_rows, num_cols);
	fillTestMatrix(mat);
	fillTestMatrix(transformedMat);
	fillTestMatrixNegative(mat2);
	fillTestMatrixNegative(transformedMat2);
	abs(transformedMat);
	abs(transformedMat2);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_EQ(transformedMat(i,j),std::abs(mat(i,j)));
	    EXPECT_EQ(transformedMat2(i,j),std::abs(mat2(i,j)));
	  }
	}
}

TEST(MT, sqrt) 
{
	const unsigned int num_rows = 4;
	const unsigned int num_cols = 3;
	
	DenseMatrix mat(num_rows, num_cols);
	DenseMatrix transformedMat(num_rows, num_cols);
	fillTestMatrix(mat);
	fillTestMatrix(transformedMat);
	sqrt(transformedMat);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_EQ(transformedMat(i,j),std::sqrt(mat(i,j)));
	  }
	}
}

TEST(MT, log) 
{
	const unsigned int num_rows = 4;
	const unsigned int num_cols = 3;
	
	DenseMatrix mat(num_rows, num_cols);
	DenseMatrix transformedMat(num_rows, num_cols);
	fillTestMatrix(mat);
	fillTestMatrix(transformedMat);
	log(transformedMat);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_EQ(transformedMat(i,j),std::log(mat(i,j)));
	  }
	}
}

TEST(MT, log2) 
{
	const unsigned int num_rows = 4;
	const unsigned int num_cols = 3;
	
	DenseMatrix mat(num_rows, num_cols);
	DenseMatrix transformedMat(num_rows, num_cols);
	fillTestMatrix(mat);
	fillTestMatrix(transformedMat);
	log2(transformedMat);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_EQ(transformedMat(i,j),std::log2(mat(i,j)));
	  }
	}
}

TEST(MT, log10) 
{
	const unsigned int num_rows = 4;
	const unsigned int num_cols = 3;
	
	DenseMatrix mat(num_rows, num_cols);
	DenseMatrix transformedMat(num_rows, num_cols);
	fillTestMatrix(mat);
	fillTestMatrix(transformedMat);
	log10(transformedMat);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_EQ(transformedMat(i,j),std::log10(mat(i,j)));
	  }
	}
}

TEST(MT, pow) 
{
	const unsigned int num_rows = 4;
	const unsigned int num_cols = 3;
	
	DenseMatrix mat(num_rows, num_cols);
	DenseMatrix mat2(num_rows, num_cols);
	DenseMatrix transformedMat(num_rows, num_cols);
	DenseMatrix transformedMat2(num_rows, num_cols);
	fillTestMatrix(mat);
	fillTestMatrix(transformedMat);
	fillTestMatrixNegative(mat2);
	fillTestMatrixNegative(transformedMat2);
	pow(transformedMat,7);
	pow(transformedMat2,7);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_EQ(transformedMat(i,j),std::pow(mat(i,j),7));
	    EXPECT_EQ(transformedMat2(i,j),std::pow(mat2(i,j),7));
	  }
	}
}

TEST(MT, pow2) 
{
const unsigned int num_rows = 4;
	const unsigned int num_cols = 3;
	
	DenseMatrix mat(num_rows, num_cols);
	DenseMatrix mat2(num_rows, num_cols);
	DenseMatrix transformedMat(num_rows, num_cols);
	DenseMatrix transformedMat2(num_rows, num_cols);
	fillTestMatrix(mat);
	fillTestMatrix(transformedMat);
	fillTestMatrixNegative(mat2);
	fillTestMatrixNegative(transformedMat2);
	pow2(transformedMat);
	pow2(transformedMat2);
	
	for(unsigned int i=0; i<num_rows; ++i) {
	  for(unsigned  int j=0; j<num_cols; ++j) {
	    EXPECT_EQ(transformedMat(i,j),std::pow(mat(i,j),2));
	    EXPECT_EQ(transformedMat2(i,j),std::pow(mat2(i,j),2));
	  }
	}
}