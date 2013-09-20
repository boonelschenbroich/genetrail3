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

#include "../src/DenseMatrix.h"
#include <config.h>

using namespace GeneTrail;

class DenseMatrixTest : public ::testing::Test
{
	public:
		DenseMatrixTest()
			: matrix_name_file(TEST_DATA_PATH("matrix_names.txt")),
			  matrix_name_additional_column_file(TEST_DATA_PATH("matrix_names_additional_column.txt")),
			  matrix_noname_file(TEST_DATA_PATH("matrix_nonames.txt"))
		{
		}

	protected:
		void constructorCheck(const DenseMatrix& mat);
		std::vector<std::string> constructorRNames();
		std::vector<std::string> constructorCNames();

		DenseMatrix buildTestMatrix();

		const std::string matrix_name_file;
		const std::string matrix_name_additional_column_file;
		const std::string matrix_noname_file;
};

DenseMatrix DenseMatrixTest::buildTestMatrix()
{
	const unsigned int num_rows = 5;
	const unsigned int num_cols = 3;

	DenseMatrix mat(num_rows, num_cols);

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

	return mat;
}

TEST_F(DenseMatrixTest, Constructor)
{
	const unsigned int num_rows = 2;
	const unsigned int num_cols = 3;

	DenseMatrix mat(num_rows, num_cols);

	ASSERT_EQ(num_rows, mat.rows());
	ASSERT_EQ(num_cols, mat.cols());

	for(unsigned int i = 0; i < num_rows; ++i) {
		ASSERT_EQ("", mat.rowName(i));
	}

	for(unsigned int i = 0; i < num_cols; ++i) {
		ASSERT_EQ("", mat.colName(i));
	}
}

TEST_F(DenseMatrixTest, NameConstructor_ll)
{
	auto row_names = constructorRNames();
	auto col_names = constructorCNames();

	DenseMatrix result(row_names, col_names);

	ASSERT_EQ(row_names.size(), result.rows());
	ASSERT_EQ(col_names.size(), result.cols());

	constructorCheck(result);
}

TEST_F(DenseMatrixTest, NameConstructor_rl)
{
	auto row_names = constructorRNames();
	auto col_names = constructorCNames();

	std::string* tmp_ptr = &row_names[0];

	DenseMatrix result(std::move(row_names), col_names);

	ASSERT_EQ(4, result.rows());
	ASSERT_EQ(5, result.cols());

	EXPECT_EQ(tmp_ptr, &result.rowNames()[0]);

	constructorCheck(result);
}

TEST_F(DenseMatrixTest, NameConstructor_lr)
{
	auto row_names = constructorRNames();
	auto col_names = constructorCNames();

	std::string* tmp_ptr = &col_names[0];

	DenseMatrix result(row_names, std::move(col_names));

	ASSERT_EQ(4, result.rows());
	ASSERT_EQ(5, result.cols());

	EXPECT_EQ(tmp_ptr, &result.colNames()[0]);

	constructorCheck(result);
}

TEST_F(DenseMatrixTest, NameConstructor_rr)
{
	auto row_names = constructorRNames();
	auto col_names = constructorCNames();

	std::string* tmp_ptr_r = &row_names[0];
	std::string* tmp_ptr_c = &col_names[0];

	DenseMatrix result(std::move(row_names), std::move(col_names));

	ASSERT_EQ(4, result.rows());
	ASSERT_EQ(5, result.cols());

	EXPECT_EQ(tmp_ptr_r, &result.rowNames()[0]);
	EXPECT_EQ(tmp_ptr_c, &result.colNames()[0]);

	constructorCheck(result);
}

std::vector<std::string> DenseMatrixTest::constructorRNames()
{
	return std::vector<std::string> {
		"row1", "row2", "row3", "row4"
	};
}

std::vector<std::string> DenseMatrixTest::constructorCNames()
{
	return std::vector<std::string> {
		"col1", "col2", "col3", "col4", "col5"
	};
}

void DenseMatrixTest::constructorCheck(const DenseMatrix& result)
{
	EXPECT_EQ("row1", result.rowName(0));
	EXPECT_EQ("row2", result.rowName(1));
	EXPECT_EQ("row3", result.rowName(2));
	EXPECT_EQ("row4", result.rowName(3));

	EXPECT_EQ("col1", result.colName(0));
	EXPECT_EQ("col2", result.colName(1));
	EXPECT_EQ("col3", result.colName(2));
	EXPECT_EQ("col4", result.colName(3));
	EXPECT_EQ("col5", result.colName(4));
}

TEST_F(DenseMatrixTest, Move_Constructor)
{
	const unsigned int num_rows = 2;
	const unsigned int num_cols = 3;

	DenseMatrix mat(num_rows, num_cols);

	// Fill the matrix
	mat(0,0) = 1.0;
	mat(0,1) = 2.0;
	mat(0,2) = 3.0;
	mat(1,0) = 4.0;
	mat(1,1) = 5.0;
	mat(1,2) = 6.0;

	mat.setRowName(0, "row0");
	mat.setRowName(1, "row1");

	mat.setColName(0, "col0");
	mat.setColName(1, "col1");
	mat.setColName(2, "col2");

	ASSERT_EQ(num_rows, mat.rows());
	ASSERT_EQ(num_cols, mat.cols());

	// Store the pointers to the internal representations
	std::vector<double*> memptr(num_rows * num_cols);
	for(unsigned int i = 0; i < num_rows; ++i) {
		for(unsigned int j = 0; j < num_cols; ++j) {
			memptr[i*num_cols + j] = &mat(i,j);
		}
	}

	std::vector<const std::string*> row_name_ptr(num_rows);
	for(unsigned int i = 0; i < num_rows; ++i) {
		row_name_ptr[i] = &mat.rowName(i);
	}

	std::vector<const std::string*> col_name_ptr(num_cols);
	for(unsigned int i = 0; i < num_cols; ++i) {
		col_name_ptr[i] = &mat.colName(i);
	}

	// Move the matrix
	DenseMatrix mat2(std::move(mat));

	// Test that the matrix was moved, not copied
	for(unsigned int i = 0; i < num_rows; ++i) {
		ASSERT_EQ(*row_name_ptr[i], mat2.rowName(i));
		EXPECT_EQ(row_name_ptr[i], &mat2.rowName(i));
	}

	for(unsigned int i = 0; i < num_cols; ++i) {
		ASSERT_EQ(*col_name_ptr[i], mat2.colName(i));
		EXPECT_EQ(col_name_ptr[i], &mat2.colName(i));
	}

	for(unsigned int i = 0; i < num_rows; ++i) {
		for(unsigned int j = 0; j < num_cols; ++j) {
			ASSERT_EQ(*memptr[i*num_cols + j], mat2(i,j));
			EXPECT_EQ(memptr[i*num_cols + j], &mat2(i,j));
		}
	}
}

TEST_F(DenseMatrixTest, Move_Assignment)
{
	DenseMatrix mat = buildTestMatrix();
	
	const unsigned int num_rows = mat.rows();
	const unsigned int num_cols = mat.cols();
	
	// Store the pointers to the internal representations
	std::vector<double*> memptr(num_rows * num_cols);
	for(unsigned int i = 0; i < num_rows; ++i) {
		for(unsigned int j = 0; j < num_cols; ++j) {
			memptr[i*num_cols + j] = &mat(i,j);
		}
	}

	std::vector<const std::string*> row_name_ptr(num_rows);
	for(unsigned int i = 0; i < num_rows; ++i) {
		row_name_ptr[i] = &mat.rowName(i);
	}

	std::vector<const std::string*> col_name_ptr(num_cols);
	for(unsigned int i = 0; i < num_cols; ++i) {
		col_name_ptr[i] = &mat.colName(i);
	}

	// Move the matrix
	DenseMatrix mat2(10,10);

	ASSERT_EQ(10, mat2.rows());
	ASSERT_EQ(10, mat2.cols());

	mat2 = std::move(mat);

	ASSERT_EQ(num_rows, mat2.rows());
	ASSERT_EQ(num_cols, mat2.cols());

	// Test that the matrix was moved, not copied
	for(unsigned int i = 0; i < num_rows; ++i) {
		ASSERT_EQ(*row_name_ptr[i], mat2.rowName(i));
		EXPECT_EQ(row_name_ptr[i], &mat2.rowName(i));
	}

	for(unsigned int i = 0; i < num_cols; ++i) {
		ASSERT_EQ(*col_name_ptr[i], mat2.colName(i));
		EXPECT_EQ(col_name_ptr[i], &mat2.colName(i));
	}

	for(unsigned int i = 0; i < num_rows; ++i) {
		for(unsigned int j = 0; j < num_cols; ++j) {
			ASSERT_EQ(*memptr[i*num_cols + j], mat2(i,j));
			EXPECT_EQ(memptr[i*num_cols + j], &mat2(i,j));
		}
	}
}

TEST_F(DenseMatrixTest, setRowName_string_)
{
	DenseMatrix mat = buildTestMatrix();
	
	ASSERT_TRUE(mat.hasRow("row1"));
	mat.setRowName("row1", "blablubb");
	ASSERT_EQ("blablubb", mat.rowName(1));
	mat.setRowName("blablubb", "blablubb2");
	ASSERT_EQ("blablubb2", mat.rowName(1));
}

TEST_F(DenseMatrixTest, setRowName_index_)
{
	DenseMatrix mat(2,2);
	
	ASSERT_EQ("", mat.rowName(0));
	ASSERT_EQ("", mat.rowName(1));
	
	mat.setRowName(0, "row0");
	mat.setRowName(1, "row1");

	EXPECT_EQ("row0", mat.rowName(0));
	EXPECT_EQ("row1", mat.rowName(1));

	EXPECT_EQ(0, mat.rowIndex("row0"));
	EXPECT_EQ(1, mat.rowIndex("row1"));
}

TEST_F(DenseMatrixTest, setColName_string_)
{
	DenseMatrix mat = buildTestMatrix();
	
	ASSERT_TRUE(mat.hasCol("col1"));
	mat.setColName("col1", "blablubb");
	ASSERT_EQ("blablubb", mat.colName(1));
	mat.setColName("blablubb", "blablubb2");
	ASSERT_EQ("blablubb2", mat.colName(1));
}

TEST_F(DenseMatrixTest, setColName_index_)
{
	DenseMatrix mat(2,2);
	
	ASSERT_EQ("", mat.colName(0));
	ASSERT_EQ("", mat.colName(1));

	mat.setColName(0, "col0");
	mat.setColName(1, "col1");

	EXPECT_EQ("col0", mat.colName(0));
	EXPECT_EQ("col1", mat.colName(1));

	EXPECT_EQ(0, mat.colIndex("col0"));
	EXPECT_EQ(1, mat.colIndex("col1"));
}

TEST_F(DenseMatrixTest, rowIndex)
{
	DenseMatrix mat(2,2);
	
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
}

TEST_F(DenseMatrixTest, colIndex)
{
	DenseMatrix mat(2,2);
	
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
}

TEST_F(DenseMatrixTest, hasRow)
{
	DenseMatrix mat = buildTestMatrix();
	
	EXPECT_TRUE(mat.hasRow("row0"));
	EXPECT_TRUE(mat.hasRow("row1"));
	EXPECT_TRUE(mat.hasRow("row2"));
	EXPECT_FALSE(mat.hasRow("blablubb"));
}

TEST_F(DenseMatrixTest, hasCol)
{
	DenseMatrix mat = buildTestMatrix();
	
	EXPECT_TRUE(mat.hasCol("col0"));
	EXPECT_TRUE(mat.hasCol("col1"));
	EXPECT_TRUE(mat.hasCol("col2"));
	EXPECT_FALSE(mat.hasCol("blablubb"));
}

TEST_F(DenseMatrixTest, row)
{
	DenseMatrix mat = buildTestMatrix();

	for(unsigned int i = 0; i < mat.rows(); ++i) {
		EXPECT_TRUE(mat.row(i) == mat.matrix().row(i));
	}

	DenseMatrix::Vector tmp(mat.cols());
	
	for(unsigned int i = 0; i < mat.cols(); ++i) {
		tmp[i] = i * 0.1;
	}

	mat.row(1) = tmp.transpose();

	ASSERT_TRUE(mat.row(1) == mat.matrix().row(1));
	EXPECT_TRUE(mat.row(1) == tmp.transpose());
}

TEST_F(DenseMatrixTest, row_const)
{
	DenseMatrix mat = buildTestMatrix();

	for(unsigned int i = 0; i < mat.rows(); ++i) {
		EXPECT_TRUE(mat.row(i) == mat.matrix().row(i));
	}
}

TEST_F(DenseMatrixTest, col)
{
	DenseMatrix mat = buildTestMatrix();

	for(unsigned int i = 0; i < mat.cols(); ++i) {
		EXPECT_TRUE(mat.col(i) == mat.matrix().col(i));
	}

	DenseMatrix::Vector tmp(mat.rows());
	
	for(unsigned int i = 0; i < mat.rows(); ++i) {
		tmp[i] = i * 0.1;
	}

	mat.col(1) = tmp;

	ASSERT_TRUE(mat.col(1) == mat.matrix().col(1));
	EXPECT_TRUE(mat.col(1) == tmp);
}

TEST_F(DenseMatrixTest, col_const)
{
	DenseMatrix mat = buildTestMatrix();

	for(unsigned int i = 0; i < mat.cols(); ++i) {
		EXPECT_TRUE(mat.col(i) == mat.matrix().col(i));
	}
}

TEST_F(DenseMatrixTest, setRow_string)
{
	DenseMatrix mat(10, 5);
	
	mat.setRowName(3, "row3");

	DenseMatrix::Vector tmp(mat.cols());
	for(unsigned int i = 0; i < mat.cols(); ++i) {
		tmp[i] = i;
	}

	mat.setRow("row3",tmp);

	EXPECT_TRUE(mat.row(3) == tmp.transpose());
}

TEST_F(DenseMatrixTest, setRow_index)
{
	DenseMatrix mat(10, 5);
	
	DenseMatrix::Vector tmp(mat.cols());
	for(unsigned int i = 0; i < mat.cols(); ++i) {
		tmp[i] = i;
	}

	mat.setRow(3, tmp);

	EXPECT_TRUE(mat.row(3) == tmp.transpose());
}

TEST_F(DenseMatrixTest, setCol_string)
{
	DenseMatrix mat(10, 5);
	
	mat.setColName(3, "col3");

	DenseMatrix::Vector tmp(mat.rows());
	for(unsigned int i = 0; i < mat.rows(); ++i) {
		tmp[i] = i;
	}

	mat.setCol("col3",tmp);

	EXPECT_TRUE(mat.col(3) == tmp);
}

TEST_F(DenseMatrixTest, setCol_index)
{
	DenseMatrix mat(10, 5);
	
	DenseMatrix::Vector tmp(mat.rows());
	for(unsigned int i = 0; i < mat.rows(); ++i) {
		tmp[i] = i;
	}

	mat.setCol(3, tmp);

	EXPECT_TRUE(mat.col(3) == tmp);
}

TEST_F(DenseMatrixTest, cols)
{
	DenseMatrix mat(7,2);

	EXPECT_EQ(2, mat.cols());
}

TEST_F(DenseMatrixTest, rows)
{
	DenseMatrix mat(7,2);

	EXPECT_EQ(7, mat.rows());
}

TEST_F(DenseMatrixTest, matrix)
{
	DenseMatrix mat(7,2);

	EXPECT_EQ(mat.cols(), mat.matrix().cols());
	EXPECT_EQ(mat.rows(), mat.matrix().rows());

	mat(0,0) = 2323.2;

	EXPECT_EQ(mat(0,0), mat.matrix()(0,0));
}

TEST_F(DenseMatrixTest, matrix_const)
{
	DenseMatrix mat(7,2);
	mat(0,0) = 2323.2;
	
	const DenseMatrix& mat2 = mat;

	EXPECT_EQ(mat2.cols(), mat2.matrix().cols());
	EXPECT_EQ(mat2.rows(), mat2.matrix().rows());


	EXPECT_EQ(mat2(0,0), mat2.matrix()(0,0));
}

TEST_F(DenseMatrixTest, transpose)
{
	DenseMatrix mat  = buildTestMatrix();
	DenseMatrix matt = buildTestMatrix();

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
}

TEST_F(DenseMatrixTest, shuffleRows)
{
	// TODO implement
	FAIL();
}

TEST_F(DenseMatrixTest, shuffleCols)
{
	// TODO implement
	FAIL();
}

TEST_F(DenseMatrixTest, removeRows)
{
	DenseMatrix matrix      = buildTestMatrix();
	DenseMatrix matrix_copy = matrix;

	ASSERT_EQ(5, matrix.rows());
	ASSERT_EQ(3, matrix.cols());

	ASSERT_EQ(0, matrix.rowIndex("row0"));
	ASSERT_EQ(1, matrix.rowIndex("row1"));
	ASSERT_EQ(2, matrix.rowIndex("row2"));
	ASSERT_EQ(3, matrix.rowIndex("row3"));
	ASSERT_EQ(4, matrix.rowIndex("row4"));

	std::vector<DenseMatrix::index_type> indices = {1, 2, 4};
	matrix.removeRows(indices);

	EXPECT_EQ(2, matrix.rows());
	EXPECT_EQ(3, matrix.cols());

	EXPECT_TRUE(matrix_copy.row(0) == matrix.row(0));
	EXPECT_TRUE(matrix_copy.row(3) == matrix.row(1));
	EXPECT_EQ(0, matrix.rowIndex("row0"));
	EXPECT_EQ(1, matrix.rowIndex("row3"));
	EXPECT_EQ("row0", matrix.rowName(0));
	EXPECT_EQ("row3", matrix.rowName(1));
}

TEST_F(DenseMatrixTest, removeCols)
{
	DenseMatrix matrix      = buildTestMatrix();
	matrix.transpose();
	DenseMatrix matrix_copy = matrix;

	ASSERT_EQ(5, matrix.cols());
	ASSERT_EQ(3, matrix.rows());

	// The names should be row0 - row4 as we transposed the matrix
	ASSERT_EQ(0, matrix.colIndex("row0"));
	ASSERT_EQ(1, matrix.colIndex("row1"));
	ASSERT_EQ(2, matrix.colIndex("row2"));
	ASSERT_EQ(3, matrix.colIndex("row3"));
	ASSERT_EQ(4, matrix.colIndex("row4"));

	std::vector<DenseMatrix::index_type> indices = {1, 2, 4};
	matrix.removeCols(indices);

	EXPECT_EQ(2, matrix.cols());
	EXPECT_EQ(3, matrix.rows());

	EXPECT_TRUE(matrix_copy.col(0) == matrix.col(0));
	EXPECT_TRUE(matrix_copy.col(3) == matrix.col(1));
	EXPECT_EQ(0, matrix.colIndex("row0"));
	EXPECT_EQ(1, matrix.colIndex("row3"));
	EXPECT_EQ("row0", matrix.colName(0));
	EXPECT_EQ("row3", matrix.colName(1));
}

TEST_F(DenseMatrixTest, rowNames)
{
	DenseMatrix mat = buildTestMatrix();

	std::vector<std::string> row_names = mat.rowNames();

	ASSERT_EQ(mat.rows(), row_names.size());

	EXPECT_EQ("row0", row_names[0]);
	EXPECT_EQ("row1", row_names[1]);
	EXPECT_EQ("row2", row_names[2]);
	EXPECT_EQ("row3", row_names[3]);
	EXPECT_EQ("row4", row_names[4]);

	mat.setRowName(2, "blablubb");
	EXPECT_EQ("blablubb", mat.rowNames()[2]);
}

TEST_F(DenseMatrixTest, colNames)
{
	DenseMatrix mat = buildTestMatrix();

	std::vector<std::string> col_names = mat.colNames();

	ASSERT_EQ(mat.cols(), col_names.size());

	EXPECT_EQ("col0", col_names[0]);
	EXPECT_EQ("col1", col_names[1]);
	EXPECT_EQ("col2", col_names[2]);

	mat.setColName(2, "blablubb");
	EXPECT_EQ("blablubb", mat.colNames()[2]);
}

TEST_F(DenseMatrixTest, setRowNames)
{
	DenseMatrix matrix(5,3);
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
}

TEST_F(DenseMatrixTest, setColNames)
{
	DenseMatrix matrix(3,5);
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
}
