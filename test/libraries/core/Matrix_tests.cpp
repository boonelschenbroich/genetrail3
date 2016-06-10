/*
 * GeneTrail2 - An efficient library for interpreting genetic data
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

#include <genetrail2/core/Exception.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseRowSubset.h>
#include <genetrail2/core/DenseColumnSubset.h>
#include <genetrail2/core/SparseMatrix.h>
#include <genetrail2/core/compat.h>
#include <config.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>

#include <memory>
#include <numeric>

using namespace GeneTrail;

std::vector<DenseMatrix> global_matrices;

template <class T>
std::unique_ptr<Matrix> buildTestMatrix();

template <class T>
std::unique_ptr<Matrix> buildTestMatrix(int, int);

void fillTestMatrix(Matrix& mat)
{
	// Fill the matrix
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

	mat.setRowName(0, "row0");
	mat.setRowName(1, "row1");
	mat.setRowName(2, "row2");
	mat.setRowName(3, "row3");
	mat.setRowName(4, "row4");

	mat.setColName(0, "col0");
	mat.setColName(1, "col1");
	mat.setColName(2, "col2");
}

template<>
std::unique_ptr<Matrix> buildTestMatrix<DenseMatrix>()
{
	const unsigned int num_rows = 5;
	const unsigned int num_cols = 3;

	auto mat = std::make_unique<DenseMatrix>(num_rows, num_cols);

	fillTestMatrix(*mat.get());

	return mat;
}

template<>
std::unique_ptr<Matrix> buildTestMatrix<DenseColumnSubset>()
{
	global_matrices.emplace_back(5, 3);
	fillTestMatrix(global_matrices.back());

	DenseColumnSubset::ISubset cols(global_matrices.back().cols());
	std::iota(cols.begin(), cols.end(), 0u);

	return std::make_unique<DenseColumnSubset>(&global_matrices.back(), std::move(cols));
}

template<>
std::unique_ptr<Matrix> buildTestMatrix<DenseRowSubset>()
{
	global_matrices.emplace_back(5, 3);
	fillTestMatrix(global_matrices.back());

	DenseRowSubset::ISubset rows(global_matrices.back().rows());
	std::iota(rows.begin(), rows.end(), 0u);

	return std::make_unique<DenseRowSubset>(&global_matrices.back(), std::move(rows));
}

template<>
std::unique_ptr<Matrix> buildTestMatrix<SparseMatrix>()
{
	const unsigned int num_rows = 5;
	const unsigned int num_cols = 3;

	auto mat = std::make_unique<SparseMatrix>(num_rows, num_cols);

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
std::unique_ptr<Matrix> buildTestMatrix<DenseMatrix>(int i, int j)
{
	return std::make_unique<DenseMatrix>(i,j);
}

template<>
std::unique_ptr<Matrix> buildTestMatrix<SparseMatrix>(int i, int j)
{
	return std::make_unique<SparseMatrix>(i,j);
}

template<>
std::unique_ptr<Matrix> buildTestMatrix<DenseColumnSubset>(int i, int j)
{
	global_matrices.emplace_back(i, j);

	DenseColumnSubset::ISubset cols(j);
	std::iota(cols.begin(), cols.end(), 0u);

	return std::make_unique<DenseColumnSubset>(&global_matrices.back(), std::move(cols));
}

template<>
std::unique_ptr<Matrix> buildTestMatrix<DenseRowSubset>(int i, int j)
{
	global_matrices.emplace_back(i, j);

	DenseRowSubset::ISubset rows(i);
	std::iota(rows.begin(), rows.end(), 0u);

	return std::make_unique<DenseRowSubset>(&global_matrices.back(), std::move(rows));
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

		void TearDown() override
		{
			global_matrices.clear();
		}

	protected:
		const std::string matrix_name_file;
		const std::string matrix_name_additional_column_file;
		const std::string matrix_noname_file;
};

typedef ::testing::Types<DenseMatrix, DenseColumnSubset, DenseRowSubset, SparseMatrix> MyTypes;
TYPED_TEST_CASE(MatrixTest, MyTypes);

TYPED_TEST(MatrixTest, setRowName_string_)
{
	auto mat = buildTestMatrix<TypeParam>();

	ASSERT_TRUE(mat->hasRow("row1"));
	mat->setRowName("row1", "blablubb");
	ASSERT_EQ("blablubb", mat->rowName(1));
	mat->setRowName("blablubb", "blablubb2");
	ASSERT_EQ("blablubb2", mat->rowName(1));
}

TYPED_TEST(MatrixTest, setRowName_index_)
{
	auto mat = buildTestMatrix<TypeParam>(2,2);

	ASSERT_EQ("", mat->rowName(0));
	ASSERT_EQ("", mat->rowName(1));

	mat->setRowName(0, "row0");
	mat->setRowName(1, "row1");

	EXPECT_EQ("row0", mat->rowName(0));
	EXPECT_EQ("row1", mat->rowName(1));

	EXPECT_EQ(0u, mat->rowIndex("row0"));
	EXPECT_EQ(1u, mat->rowIndex("row1"));
}

TYPED_TEST(MatrixTest, setColName_string_)
{
	auto mat = buildTestMatrix<TypeParam>();

	ASSERT_TRUE(mat->hasCol("col1"));
	mat->setColName("col1", "blablubb");
	ASSERT_EQ("blablubb", mat->colName(1));
	mat->setColName("blablubb", "blablubb2");
	ASSERT_EQ("blablubb2", mat->colName(1));
}

TYPED_TEST(MatrixTest, setColName_index_)
{
	auto mat = buildTestMatrix<TypeParam>(2,2);

	ASSERT_EQ("", mat->colName(0));
	ASSERT_EQ("", mat->colName(1));

	mat->setColName(0, "col0");
	mat->setColName(1, "col1");

	EXPECT_EQ("col0", mat->colName(0));
	EXPECT_EQ("col1", mat->colName(1));

	EXPECT_EQ(0u, mat->colIndex("col0"));
	EXPECT_EQ(1u, mat->colIndex("col1"));
}

TYPED_TEST(MatrixTest, rowIndex)
{
	auto mat = buildTestMatrix<TypeParam>(2,2);

	ASSERT_EQ("", mat->rowName(0));
	ASSERT_EQ("", mat->rowName(1));

	mat->setRowName(0, "row0");
	mat->setRowName(1, "row1");

	EXPECT_EQ(0u, mat->rowIndex("row0"));
	EXPECT_EQ(1u, mat->rowIndex("row1"));

	mat->setRowName(0, "row1");
	mat->setRowName(1, "row0");

	EXPECT_EQ(0u, mat->rowIndex("row1"));
	EXPECT_EQ(1u, mat->rowIndex("row0"));
}

TYPED_TEST(MatrixTest, colIndex)
{
	auto mat = buildTestMatrix<TypeParam>(2,2);

	ASSERT_EQ("", mat->colName(0));
	ASSERT_EQ("", mat->colName(1));

	mat->setColName(0, "col0");
	mat->setColName(1, "col1");

	ASSERT_EQ("col0", mat->colName(0));
	ASSERT_EQ("col1", mat->colName(1));

	EXPECT_EQ(0u, mat->colIndex("col0"));
	EXPECT_EQ(1u, mat->colIndex("col1"));

	mat->setColName(0, "col1");
	mat->setColName(1, "col0");

	ASSERT_EQ("col1", mat->colName(0));
	ASSERT_EQ("col0", mat->colName(1));

	EXPECT_EQ(0u, mat->colIndex("col1"));
	EXPECT_EQ(1u, mat->colIndex("col0"));
}

TYPED_TEST(MatrixTest, hasRow)
{
	auto mat = buildTestMatrix<TypeParam>();

	EXPECT_TRUE(mat->hasRow("row0"));
	EXPECT_TRUE(mat->hasRow("row1"));
	EXPECT_TRUE(mat->hasRow("row2"));
	EXPECT_FALSE(mat->hasRow("blablubb"));
}

TYPED_TEST(MatrixTest, hasCol)
{
	auto mat = buildTestMatrix<TypeParam>();

	EXPECT_TRUE(mat->hasCol("col0"));
	EXPECT_TRUE(mat->hasCol("col1"));
	EXPECT_TRUE(mat->hasCol("col2"));
	EXPECT_FALSE(mat->hasCol("blablubb"));
}

TYPED_TEST(MatrixTest, cols)
{
	auto mat = buildTestMatrix<TypeParam>(7, 2);

	EXPECT_EQ(2u, mat->cols());
}

TYPED_TEST(MatrixTest, rows)
{
	auto mat = buildTestMatrix<TypeParam>(7, 2);

	EXPECT_EQ(7u, mat->rows());
}

TYPED_TEST(MatrixTest, transpose)
{
	auto mat  = buildTestMatrix<TypeParam>();
	auto matt = buildTestMatrix<TypeParam>();

	try {
		matt->transpose();
	} catch (NotImplemented&) {
		return;
	}

	ASSERT_EQ(mat->rows(), matt->cols());
	ASSERT_EQ(mat->cols(), matt->rows());

	for(unsigned int i = 0; i < mat->rows(); ++i) {
		for(unsigned int j = 0; j < mat->cols(); ++j) {
			EXPECT_EQ((*mat)(i, j), (*matt)(j, i));
		}
	}

	for(unsigned int i = 0; i < mat->rows(); ++i) {
		EXPECT_EQ(mat->rowName(i), matt->colName(i));
	}

	for(unsigned int i = 0; i < mat->cols(); ++i) {
		EXPECT_EQ(mat->colName(i), matt->rowName(i));
	}
}

TYPED_TEST(MatrixTest, shuffleRows)
{
	auto mat      = buildTestMatrix<TypeParam>();
	auto mat_orig = buildTestMatrix<TypeParam>();

	ASSERT_EQ(5u, mat->rows());
	ASSERT_EQ(3u, mat->cols());

	std::vector<Matrix::index_type> perm = {
		2, 1, 3, 4, 0
	};

	try {
		mat->shuffleRows(perm);
	} catch(NotImplemented&) {
		return;
	}

	ASSERT_EQ(5u, mat->rows());
	ASSERT_EQ(3u, mat->cols());

	for(unsigned int i = 0; i < mat_orig->rows(); ++i) {
		EXPECT_EQ(mat_orig->rowName(perm[i]), mat->rowName(i));

		for(unsigned int j = 0; j < mat_orig->cols(); ++j) {
			EXPECT_EQ((*mat_orig)(perm[i], j), (*mat)(i, j));
		}
	}
}

TYPED_TEST(MatrixTest, shuffleCols)
{
	auto mat      = buildTestMatrix<TypeParam>();
	auto mat_orig = buildTestMatrix<TypeParam>();

	ASSERT_EQ(5u, mat->rows());
	ASSERT_EQ(3u, mat->cols());

	std::vector<Matrix::index_type> perm = {
		2, 0, 1
	};

	try {
		mat->shuffleCols(perm);
	} catch(NotImplemented&) {
		return;
	}

	ASSERT_EQ(5u, mat->rows());
	ASSERT_EQ(3u, mat->cols());

	for(Matrix::index_type j = 0; j < mat_orig->cols(); ++j) {
		EXPECT_EQ(mat_orig->colName(perm[j]), mat->colName(j));

		for(Matrix::index_type i = 0; i < mat_orig->rows(); ++i) {
			EXPECT_EQ((*mat_orig)(i, perm[j]), (*mat)(i, j));
		}
	}
}

TYPED_TEST(MatrixTest, removeRows)
{
	auto matrix      = buildTestMatrix<TypeParam>();
	auto matrix_copy = buildTestMatrix<TypeParam>();

	ASSERT_EQ(5u, matrix->rows());
	ASSERT_EQ(3u, matrix->cols());

	ASSERT_EQ(0u, matrix->rowIndex("row0"));
	ASSERT_EQ(1u, matrix->rowIndex("row1"));
	ASSERT_EQ(2u, matrix->rowIndex("row2"));
	ASSERT_EQ(3u, matrix->rowIndex("row3"));
	ASSERT_EQ(4u, matrix->rowIndex("row4"));

	std::vector<Matrix::index_type> indices = {1, 2, 4};
	try {
		matrix->removeRows(indices);
	} catch(NotImplemented&) {
		return;
	}

	EXPECT_EQ(2u, matrix->rows());
	EXPECT_EQ(3u, matrix->cols());

	for(Matrix::index_type i = 0; i < matrix->cols(); ++i) {
		EXPECT_EQ((*matrix_copy)(0, i), (*matrix)(0, i));
		EXPECT_EQ((*matrix_copy)(3, i), (*matrix)(1, i));
	}

	EXPECT_EQ(0u, matrix->rowIndex("row0"));
	EXPECT_EQ(1u, matrix->rowIndex("row3"));
	EXPECT_EQ("row0", matrix->rowName(0));
	EXPECT_EQ("row3", matrix->rowName(1));
}

TYPED_TEST(MatrixTest, removeCols)
{
	auto matrix      = buildTestMatrix<TypeParam>(3, 5);
	auto matrix_copy = buildTestMatrix<TypeParam>(3, 5);

	matrix->setColName(0, "col0");
	matrix->setColName(1, "col1");
	matrix->setColName(2, "col2");
	matrix->setColName(3, "col3");
	matrix->setColName(4, "col4");

	// Fill the matrix with unique values
	for(Matrix::index_type i = 0; i < matrix->rows(); ++i) {
		for(Matrix::index_type j = 0; j < matrix->cols(); ++j) {
			(*matrix)(i,j) = (*matrix_copy)(i, j) = i * 5 + j;
		}
	}

	ASSERT_EQ(5u, matrix->cols());
	ASSERT_EQ(3u, matrix->rows());

	// The names should be row0 - row4 as we transposed the matrix
	ASSERT_EQ(0u, matrix->colIndex("col0"));
	ASSERT_EQ(1u, matrix->colIndex("col1"));
	ASSERT_EQ(2u, matrix->colIndex("col2"));
	ASSERT_EQ(3u, matrix->colIndex("col3"));
	ASSERT_EQ(4u, matrix->colIndex("col4"));

	std::vector<Matrix::index_type> indices = {1, 2, 4};
	try {
		matrix->removeCols(indices);
	} catch(NotImplemented&) {
		return;
	}

	EXPECT_EQ(2u, matrix->cols());
	EXPECT_EQ(3u, matrix->rows());

	for(Matrix::index_type i = 0; i < matrix->rows(); ++i) {
		EXPECT_EQ((*matrix_copy)(i, 0), (*matrix)(i, 0));
		EXPECT_EQ((*matrix_copy)(i, 3), (*matrix)(i, 1));
	}

	EXPECT_EQ(0u, matrix->colIndex("col0"));
	EXPECT_EQ(1u, matrix->colIndex("col3"));
	EXPECT_EQ("col0", matrix->colName(0));
	EXPECT_EQ("col3", matrix->colName(1));
}

TYPED_TEST(MatrixTest, rowNames)
{
	auto mat = buildTestMatrix<TypeParam>();

	std::vector<std::string> row_names = mat->rowNames();

	ASSERT_EQ(mat->rows(), row_names.size());

	EXPECT_EQ("row0", row_names[0]);
	EXPECT_EQ("row1", row_names[1]);
	EXPECT_EQ("row2", row_names[2]);
	EXPECT_EQ("row3", row_names[3]);
	EXPECT_EQ("row4", row_names[4]);

	mat->setRowName(2, "blablubb");
	EXPECT_EQ("blablubb", mat->rowNames()[2]);
}

TYPED_TEST(MatrixTest, colNames)
{
	auto mat = buildTestMatrix<TypeParam>();

	std::vector<std::string> col_names = mat->colNames();

	ASSERT_EQ(mat->cols(), col_names.size());

	EXPECT_EQ("col0", col_names[0]);
	EXPECT_EQ("col1", col_names[1]);
	EXPECT_EQ("col2", col_names[2]);

	mat->setColName(2, "blablubb");
	EXPECT_EQ("blablubb", mat->colNames()[2]);
}

TYPED_TEST(MatrixTest, setRowNames)
{
	auto matrix = buildTestMatrix<TypeParam>(5,3);

	std::vector<std::string> row_names = {
		"row0", "row1", "row2", "row3", "row4"
	};

	ASSERT_EQ(5u, matrix->rows());
	EXPECT_EQ("", matrix->rowName(0));
	EXPECT_EQ("", matrix->rowName(1));
	EXPECT_EQ("", matrix->rowName(2));
	EXPECT_EQ("", matrix->rowName(3));
	EXPECT_EQ("", matrix->rowName(4));

	matrix->setRowNames(row_names);

	EXPECT_EQ("row0", matrix->rowName(0));
	EXPECT_EQ("row1", matrix->rowName(1));
	EXPECT_EQ("row2", matrix->rowName(2));
	EXPECT_EQ("row3", matrix->rowName(3));
	EXPECT_EQ("row4", matrix->rowName(4));

	EXPECT_EQ(0u, matrix->rowIndex("row0"));
	EXPECT_EQ(1u, matrix->rowIndex("row1"));
	EXPECT_EQ(2u, matrix->rowIndex("row2"));
	EXPECT_EQ(3u, matrix->rowIndex("row3"));
	EXPECT_EQ(4u, matrix->rowIndex("row4"));
}

TYPED_TEST(MatrixTest, setColNames)
{
	auto matrix = buildTestMatrix<TypeParam>(3,5);

	std::vector<std::string> col_names = {
		"col0", "col1", "col2", "col3", "col4"
	};

	ASSERT_EQ(5u, matrix->cols());
	EXPECT_EQ("", matrix->colName(0));
	EXPECT_EQ("", matrix->colName(1));
	EXPECT_EQ("", matrix->colName(2));
	EXPECT_EQ("", matrix->colName(3));
	EXPECT_EQ("", matrix->colName(4));

	matrix->setColNames(col_names);

	EXPECT_EQ("col0", matrix->colName(0));
	EXPECT_EQ("col1", matrix->colName(1));
	EXPECT_EQ("col2", matrix->colName(2));
	EXPECT_EQ("col3", matrix->colName(3));
	EXPECT_EQ("col4", matrix->colName(4));

	EXPECT_EQ(0u, matrix->colIndex("col0"));
	EXPECT_EQ(1u, matrix->colIndex("col1"));
	EXPECT_EQ(2u, matrix->colIndex("col2"));
	EXPECT_EQ(3u, matrix->colIndex("col3"));
	EXPECT_EQ(4u, matrix->colIndex("col4"));
}
