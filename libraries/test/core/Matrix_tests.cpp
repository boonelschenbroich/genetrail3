/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2013 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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

#include <gtest/gtest.h>

#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixSubset.h>
#include <genetrail2/core/SparseMatrix.h>
#include <config.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>

using namespace GeneTrail;

DenseMatrix global_matrix(0, 0);

template <class T>
Matrix* buildTestMatrix();

template <class T>
Matrix* buildTestMatrix(int, int);

template<>
Matrix* buildTestMatrix<DenseMatrix>()
{
	const unsigned int num_rows = 5;
	const unsigned int num_cols = 3;

	DenseMatrix* mat = new DenseMatrix(num_rows, num_cols);

	// Fill the matrix
	(*mat)(0,0) =  1.0;
	(*mat)(0,1) =  2.0;
	(*mat)(0,2) =  3.0;
	(*mat)(1,0) =  4.0;
	(*mat)(1,1) =  5.0;
	(*mat)(1,2) =  6.0;
	(*mat)(2,0) =  7.0;
	(*mat)(2,1) =  8.0;
	(*mat)(2,2) =  9.0;
	(*mat)(3,0) = 10.0;
	(*mat)(3,1) = 11.0;
	(*mat)(3,2) = 12.0;
	(*mat)(4,0) = 13.0;
	(*mat)(4,1) = 14.0;
	(*mat)(4,2) = 15.0;

	mat->setRowName(0, "row0");
	mat->setRowName(1, "row1");
	mat->setRowName(2, "row2");
	mat->setRowName(3, "row3");
	mat->setRowName(4, "row4");

	mat->setColName(0, "col0");
	mat->setColName(1, "col1");
	mat->setColName(2, "col2");

	return mat;
}

template<>
Matrix* buildTestMatrix<DenseMatrixSubset>()
{
	global_matrix = DenseMatrix(5,3);

	// Fill the matrix
	global_matrix(0,0) =  1.0;
	global_matrix(0,1) =  2.0;
	global_matrix(0,2) =  3.0;
	global_matrix(1,0) =  4.0;
	global_matrix(1,1) =  5.0;
	global_matrix(1,2) =  6.0;
	global_matrix(2,0) =  7.0;
	global_matrix(2,1) =  8.0;
	global_matrix(2,2) =  9.0;
	global_matrix(3,0) = 10.0;
	global_matrix(3,1) = 11.0;
	global_matrix(3,2) = 12.0;
	global_matrix(4,0) = 13.0;
	global_matrix(4,1) = 14.0;
	global_matrix(4,2) = 15.0;

	global_matrix.setRowName(0, "row0");
	global_matrix.setRowName(1, "row1");
	global_matrix.setRowName(2, "row2");
	global_matrix.setRowName(3, "row3");
	global_matrix.setRowName(4, "row4");

	global_matrix.setColName(0, "col0");
	global_matrix.setColName(1, "col1");
	global_matrix.setColName(2, "col2");

	DenseMatrixSubset::ISubset rows(5), cols(3);

	for(int k = 0; k < 5; ++k)
	{
		rows[k] = k;
	}

	for(int k = 0; k < 3; ++k)
	{
		cols[k] = k;
	}

	return new DenseMatrixSubset(&global_matrix, std::move(rows), std::move(cols));
}

template<>
Matrix* buildTestMatrix<SparseMatrix>()
{
	const unsigned int num_rows = 5;
	const unsigned int num_cols = 3;

	GeneTrail::SparseMatrix* mat = new GeneTrail::SparseMatrix(num_rows, num_cols);

	typedef Eigen::Triplet<GeneTrail::SparseMatrix::value_type> T;

	std::vector<T> tmp{
		T(0,0, 1.0),
		T(0,1, 2.0),
		T(0,2, 3.0),
		T(1,0, 4.0),
		T(1,1, 5.0),
		T(1,2, 6.0),
		T(2,0, 7.0),
		T(2,1, 8.0),
		T(2,2, 9.0),
		T(3,0,10.0),
		T(3,1,11.0),
		T(3,2,12.0),
		T(4,0,13.0),
		T(4,1,14.0),
		T(4,2,15.0)
	};

	// Fill the matrix
	mat->matrix().setFromTriplets(tmp.begin(), tmp.end());

	mat->setRowName(0, "row0");
	mat->setRowName(1, "row1");
	mat->setRowName(2, "row2");
	mat->setRowName(3, "row3");
	mat->setRowName(4, "row4");

	mat->setColName(0, "col0");
	mat->setColName(1, "col1");
	mat->setColName(2, "col2");

	return mat;
}

template<>
Matrix* buildTestMatrix<DenseMatrix>(int i, int j)
{
	return new DenseMatrix(i,j);
}

template<>
Matrix* buildTestMatrix<SparseMatrix>(int i, int j)
{
	return new SparseMatrix(i,j);
}

template<>
Matrix* buildTestMatrix<DenseMatrixSubset>(int i, int j)
{
	global_matrix = DenseMatrix(i, j);

	DenseMatrixSubset::ISubset rows(i), cols(j);

	for(int k = 0; k < i; ++k)
	{
		rows[k] = k;
	}

	for(int k = 0; k < j; ++k)
	{
		cols[k] = k;
	}

	return new DenseMatrixSubset(&global_matrix, std::move(rows), std::move(cols));
}

template <typename T>
class MatrixTest : public ::testing::Test
{
	public:
		MatrixTest()
			: matrix_name_file(TEST_DATA_PATH("matrix_names.txt")),
			  matrix_name_additional_column_file(TEST_DATA_PATH("matrix_names_additional_column.txt")),
			  matrix_noname_file(TEST_DATA_PATH("matrix_nonames.txt"))
		{
		}

	protected:
		const std::string matrix_name_file;
		const std::string matrix_name_additional_column_file;
		const std::string matrix_noname_file;
};

typedef ::testing::Types<DenseMatrix, DenseMatrixSubset, SparseMatrix> MyTypes;
TYPED_TEST_CASE(MatrixTest, MyTypes);

TYPED_TEST(MatrixTest, setRowName_string_)
{
	Matrix& mat = *buildTestMatrix<TypeParam>();

	ASSERT_TRUE(mat.hasRow("row1"));
	mat.setRowName("row1", "blablubb");
	ASSERT_EQ("blablubb", mat.rowName(1));
	mat.setRowName("blablubb", "blablubb2");
	ASSERT_EQ("blablubb2", mat.rowName(1));

	delete &mat;
}

TYPED_TEST(MatrixTest, setRowName_index_)
{
	Matrix& mat = *buildTestMatrix<TypeParam>(2,2);

	ASSERT_EQ("", mat.rowName(0));
	ASSERT_EQ("", mat.rowName(1));

	mat.setRowName(0, "row0");
	mat.setRowName(1, "row1");

	EXPECT_EQ("row0", mat.rowName(0));
	EXPECT_EQ("row1", mat.rowName(1));

	EXPECT_EQ(0, mat.rowIndex("row0"));
	EXPECT_EQ(1, mat.rowIndex("row1"));
}

TYPED_TEST(MatrixTest, setColName_string_)
{
	Matrix& mat = *buildTestMatrix<TypeParam>();

	ASSERT_TRUE(mat.hasCol("col1"));
	mat.setColName("col1", "blablubb");
	ASSERT_EQ("blablubb", mat.colName(1));
	mat.setColName("blablubb", "blablubb2");
	ASSERT_EQ("blablubb2", mat.colName(1));

	delete &mat;
}

TYPED_TEST(MatrixTest, setColName_index_)
{
	Matrix& mat = *buildTestMatrix<TypeParam>(2,2);

	ASSERT_EQ("", mat.colName(0));
	ASSERT_EQ("", mat.colName(1));

	mat.setColName(0, "col0");
	mat.setColName(1, "col1");

	EXPECT_EQ("col0", mat.colName(0));
	EXPECT_EQ("col1", mat.colName(1));

	EXPECT_EQ(0, mat.colIndex("col0"));
	EXPECT_EQ(1, mat.colIndex("col1"));

	delete &mat;
}

TYPED_TEST(MatrixTest, rowIndex)
{
	Matrix& mat = *buildTestMatrix<TypeParam>(2,2);

	ASSERT_EQ("", mat.rowName(0));
	ASSERT_EQ("", mat.rowName(1));

	mat.setRowName(0, "row0");
	mat.setRowName(1, "row1");

	EXPECT_EQ(0, mat.rowIndex("row0"));
	EXPECT_EQ(1, mat.rowIndex("row1"));

	mat.setRowName(0, "row1");
	mat.setRowName(1, "row0");

	EXPECT_EQ(0, mat.rowIndex("row1"));
	EXPECT_EQ(1, mat.rowIndex("row0"));

	delete &mat;
}

TYPED_TEST(MatrixTest, colIndex)
{
	Matrix& mat = *buildTestMatrix<TypeParam>(2,2);

	ASSERT_EQ("", mat.colName(0));
	ASSERT_EQ("", mat.colName(1));

	mat.setColName(0, "col0");
	mat.setColName(1, "col1");

	ASSERT_EQ("col0", mat.colName(0));
	ASSERT_EQ("col1", mat.colName(1));

	EXPECT_EQ(0, mat.colIndex("col0"));
	EXPECT_EQ(1, mat.colIndex("col1"));

	mat.setColName(0, "col1");
	mat.setColName(1, "col0");

	ASSERT_EQ("col1", mat.colName(0));
	ASSERT_EQ("col0", mat.colName(1));

	EXPECT_EQ(0, mat.colIndex("col1"));
	EXPECT_EQ(1, mat.colIndex("col0"));

	delete &mat;
}

TYPED_TEST(MatrixTest, hasRow)
{
	Matrix& mat = *buildTestMatrix<TypeParam>();

	EXPECT_TRUE(mat.hasRow("row0"));
	EXPECT_TRUE(mat.hasRow("row1"));
	EXPECT_TRUE(mat.hasRow("row2"));
	EXPECT_FALSE(mat.hasRow("blablubb"));

	delete &mat;
}

TYPED_TEST(MatrixTest, hasCol)
{
	Matrix& mat = *buildTestMatrix<TypeParam>();

	EXPECT_TRUE(mat.hasCol("col0"));
	EXPECT_TRUE(mat.hasCol("col1"));
	EXPECT_TRUE(mat.hasCol("col2"));
	EXPECT_FALSE(mat.hasCol("blablubb"));

	delete &mat;
}

TYPED_TEST(MatrixTest, cols)
{
	Matrix& mat = *buildTestMatrix<TypeParam>(7, 2);

	EXPECT_EQ(2, mat.cols());

	delete &mat;
}

TYPED_TEST(MatrixTest, rows)
{
	Matrix& mat = *buildTestMatrix<TypeParam>(7, 2);

	EXPECT_EQ(7, mat.rows());

	delete &mat;
}

TYPED_TEST(MatrixTest, transpose)
{
	Matrix& mat  = *buildTestMatrix<TypeParam>();
	Matrix& matt = *buildTestMatrix<TypeParam>();

	matt.transpose();

	ASSERT_EQ(mat.rows(), matt.cols());
	ASSERT_EQ(mat.cols(), matt.rows());

	for(unsigned int i = 0; i < mat.rows(); ++i) {
		for(unsigned int j = 0; j < mat.cols(); ++j) {
			EXPECT_EQ(mat(i, j), matt(j, i));
		}
	}

	for(unsigned int i = 0; i < mat.rows(); ++i) {
		EXPECT_EQ(mat.rowName(i), matt.colName(i));
	}

	for(unsigned int i = 0; i < mat.cols(); ++i) {
		EXPECT_EQ(mat.colName(i), matt.rowName(i));
	}

	delete &matt;
	delete &mat;
}

TYPED_TEST(MatrixTest, shuffleRows)
{
	Matrix& mat      = *buildTestMatrix<TypeParam>();
	Matrix& mat_orig = *buildTestMatrix<TypeParam>();

	ASSERT_EQ(5, mat.rows());
	ASSERT_EQ(3, mat.cols());

	std::vector<Matrix::index_type> perm = {
		2, 1, 3, 4, 0
	};

	mat.shuffleRows(perm);

	ASSERT_EQ(5, mat.rows());
	ASSERT_EQ(3, mat.cols());

	for(int i = 0; i < mat_orig.rows(); ++i) {
		EXPECT_EQ(mat_orig.rowName(perm[i]), mat.rowName(i));

		for(int j = 0; j < mat_orig.cols(); ++j) {
			EXPECT_EQ(mat_orig(perm[i], j), mat(i, j));
		}
	}

	delete &mat_orig;
	delete &mat;
}

TYPED_TEST(MatrixTest, shuffleCols)
{
	Matrix& mat      = *buildTestMatrix<TypeParam>();
	Matrix& mat_orig = *buildTestMatrix<TypeParam>();

	ASSERT_EQ(5, mat.rows());
	ASSERT_EQ(3, mat.cols());

	std::vector<Matrix::index_type> perm = {
		2, 0, 1
	};

	mat.shuffleCols(perm);

	ASSERT_EQ(5, mat.rows());
	ASSERT_EQ(3, mat.cols());

	for(int j = 0; j < mat_orig.cols(); ++j) {
		EXPECT_EQ(mat_orig.colName(perm[j]), mat.colName(j));

		for(int i = 0; i < mat_orig.rows(); ++i) {
			EXPECT_EQ(mat_orig(i, perm[j]), mat(i, j));
		}
	}

	delete &mat_orig;
	delete &mat;
}

TYPED_TEST(MatrixTest, removeRows)
{
	Matrix& matrix = *buildTestMatrix<TypeParam>();
	Matrix& matrix_copy = *buildTestMatrix<TypeParam>();

	ASSERT_EQ(5, matrix.rows());
	ASSERT_EQ(3, matrix.cols());

	ASSERT_EQ(0, matrix.rowIndex("row0"));
	ASSERT_EQ(1, matrix.rowIndex("row1"));
	ASSERT_EQ(2, matrix.rowIndex("row2"));
	ASSERT_EQ(3, matrix.rowIndex("row3"));
	ASSERT_EQ(4, matrix.rowIndex("row4"));

	std::vector<Matrix::index_type> indices = {1, 2, 4};
	matrix.removeRows(indices);

	EXPECT_EQ(2, matrix.rows());
	EXPECT_EQ(3, matrix.cols());

	for(int i = 0; i < matrix.cols(); ++i) {
		EXPECT_EQ(matrix_copy(0, i), matrix(0, i));
		EXPECT_EQ(matrix_copy(3, i), matrix(1, i));
	}

	EXPECT_EQ(0, matrix.rowIndex("row0"));
	EXPECT_EQ(1, matrix.rowIndex("row3"));
	EXPECT_EQ("row0", matrix.rowName(0));
	EXPECT_EQ("row3", matrix.rowName(1));

	delete &matrix_copy;
	delete &matrix;
}

TYPED_TEST(MatrixTest, removeCols)
{
	Matrix& matrix = *buildTestMatrix<TypeParam>();
	Matrix& matrix_copy = *buildTestMatrix<TypeParam>();

	matrix.transpose();
	matrix_copy.transpose();

	ASSERT_EQ(5, matrix.cols());
	ASSERT_EQ(3, matrix.rows());

	// The names should be row0 - row4 as we transposed the matrix
	ASSERT_EQ(0, matrix.colIndex("row0"));
	ASSERT_EQ(1, matrix.colIndex("row1"));
	ASSERT_EQ(2, matrix.colIndex("row2"));
	ASSERT_EQ(3, matrix.colIndex("row3"));
	ASSERT_EQ(4, matrix.colIndex("row4"));

	std::vector<Matrix::index_type> indices = {1, 2, 4};
	matrix.removeCols(indices);

	EXPECT_EQ(2, matrix.cols());
	EXPECT_EQ(3, matrix.rows());

	for(int i = 0; i < matrix.rows(); ++i) {
		EXPECT_EQ(matrix_copy(i, 0), matrix(i, 0));
		EXPECT_EQ(matrix_copy(i, 3), matrix(i, 1));
	}

	EXPECT_EQ(0, matrix.colIndex("row0"));
	EXPECT_EQ(1, matrix.colIndex("row3"));
	EXPECT_EQ("row0", matrix.colName(0));
	EXPECT_EQ("row3", matrix.colName(1));

	delete &matrix_copy;
	delete &matrix;
}

TYPED_TEST(MatrixTest, rowNames)
{
	Matrix& mat = *buildTestMatrix<TypeParam>();

	std::vector<std::string> row_names = mat.rowNames();

	ASSERT_EQ(mat.rows(), row_names.size());

	EXPECT_EQ("row0", row_names[0]);
	EXPECT_EQ("row1", row_names[1]);
	EXPECT_EQ("row2", row_names[2]);
	EXPECT_EQ("row3", row_names[3]);
	EXPECT_EQ("row4", row_names[4]);

	mat.setRowName(2, "blablubb");
	EXPECT_EQ("blablubb", mat.rowNames()[2]);

	delete &mat;
}

TYPED_TEST(MatrixTest, colNames)
{
	Matrix& mat = *buildTestMatrix<TypeParam>();

	std::vector<std::string> col_names = mat.colNames();

	ASSERT_EQ(mat.cols(), col_names.size());

	EXPECT_EQ("col0", col_names[0]);
	EXPECT_EQ("col1", col_names[1]);
	EXPECT_EQ("col2", col_names[2]);

	mat.setColName(2, "blablubb");
	EXPECT_EQ("blablubb", mat.colNames()[2]);

	delete &mat;
}

TYPED_TEST(MatrixTest, setRowNames)
{
	Matrix& matrix = *buildTestMatrix<TypeParam>(5,3);

	std::vector<std::string> row_names = {
		"row0", "row1", "row2", "row3", "row4"
	};

	ASSERT_EQ(5, matrix.rows());
	EXPECT_EQ("", matrix.rowName(0));
	EXPECT_EQ("", matrix.rowName(1));
	EXPECT_EQ("", matrix.rowName(2));
	EXPECT_EQ("", matrix.rowName(3));
	EXPECT_EQ("", matrix.rowName(4));

	matrix.setRowNames(row_names);

	EXPECT_EQ("row0", matrix.rowName(0));
	EXPECT_EQ("row1", matrix.rowName(1));
	EXPECT_EQ("row2", matrix.rowName(2));
	EXPECT_EQ("row3", matrix.rowName(3));
	EXPECT_EQ("row4", matrix.rowName(4));

	EXPECT_EQ(0, matrix.rowIndex("row0"));
	EXPECT_EQ(1, matrix.rowIndex("row1"));
	EXPECT_EQ(2, matrix.rowIndex("row2"));
	EXPECT_EQ(3, matrix.rowIndex("row3"));
	EXPECT_EQ(4, matrix.rowIndex("row4"));

	delete &matrix;
}

TYPED_TEST(MatrixTest, setColNames)
{
	Matrix& matrix = *buildTestMatrix<TypeParam>(3,5);

	std::vector<std::string> col_names = {
		"col0", "col1", "col2", "col3", "col4"
	};

	ASSERT_EQ(5, matrix.cols());
	EXPECT_EQ("", matrix.colName(0));
	EXPECT_EQ("", matrix.colName(1));
	EXPECT_EQ("", matrix.colName(2));
	EXPECT_EQ("", matrix.colName(3));
	EXPECT_EQ("", matrix.colName(4));

	matrix.setColNames(col_names);

	EXPECT_EQ("col0", matrix.colName(0));
	EXPECT_EQ("col1", matrix.colName(1));
	EXPECT_EQ("col2", matrix.colName(2));
	EXPECT_EQ("col3", matrix.colName(3));
	EXPECT_EQ("col4", matrix.colName(4));

	EXPECT_EQ(0, matrix.colIndex("col0"));
	EXPECT_EQ(1, matrix.colIndex("col1"));
	EXPECT_EQ(2, matrix.colIndex("col2"));
	EXPECT_EQ(3, matrix.colIndex("col3"));
	EXPECT_EQ(4, matrix.colIndex("col4"));

	delete &matrix;
}
