/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2016 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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

#include <genetrail2/core/DenseColumnSubset.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixIterator.h>
#include <genetrail2/core/DenseRowSubset.h>

using namespace GeneTrail;

class DenseMatrixIteratorTest : public ::testing::Test
{
};

DenseMatrix makeMatrix(size_t rows, size_t cols) {
	DenseMatrix mat(rows, cols);

	for(DenseMatrix::index_type i = 0; i < mat.rows(); ++i) {
		for(DenseMatrix::index_type j = 0; j < mat.cols(); ++j) {
			mat(i,j) = i*mat.cols() + j;
		}
	}

	return mat;
}

TEST_F(DenseMatrixIteratorTest, compile_column_iterator)
{
	DenseMatrix mat(10,20);
	make_column_iterator(&mat, 0);
}

TEST_F(DenseMatrixIteratorTest, compile_row_iterator)
{
	auto mat = makeMatrix(10, 20);
	make_row_iterator(&mat, 0);
}

TEST_F(DenseMatrixIteratorTest, compile_column_iterator_subset)
{
	auto mat = makeMatrix(10, 20);
	DenseColumnSubset subs(&mat, {0,4,6});
	make_column_iterator(&subs, 0);
}

TEST_F(DenseMatrixIteratorTest, compile_row_iterator_subset)
{
	auto mat = makeMatrix(10, 20);
	DenseRowSubset subs(&mat, {0,4,6});
	make_row_iterator(&subs, 0);
}

TEST_F(DenseMatrixIteratorTest, basic_column_iteration)
{
	auto mat = makeMatrix(10, 20);
	auto first = make_column_iterator(&mat, 0);
	auto last  = make_column_iterator(&mat, mat.cols());

	for(DenseMatrix::index_type j = 0; j < mat.cols(); ++j, ++first) {
		EXPECT_EQ(mat.col(j), *first);
	}

	EXPECT_TRUE(first == last);
}

TEST_F(DenseMatrixIteratorTest, basic_row_iteration)
{
	auto mat = makeMatrix(10, 20);
	auto first = make_row_iterator(&mat, 0);
	auto last  = make_row_iterator(&mat, mat.rows());

	for(DenseMatrix::index_type i = 0; i < mat.rows(); ++i, ++first) {
		EXPECT_EQ(mat.row(i), *first);
	}

	EXPECT_TRUE(first == last);
}

TEST_F(DenseMatrixIteratorTest, basic_column_iteration_subset)
{
	auto mat = makeMatrix(10, 20);
	DenseColumnSubset subs(&mat, {0,4,6});
	auto first = make_column_iterator(&subs, 0);
	auto last  = make_column_iterator(&subs, subs.cols());

	for(DenseMatrix::index_type j = 0; j < subs.cols(); ++j, ++first) {
		EXPECT_EQ(subs.col(j), *first);
	}

	EXPECT_TRUE(first == last);
}

TEST_F(DenseMatrixIteratorTest, basic_row_iteration_subset)
{
	auto mat = makeMatrix(10, 20);
	DenseRowSubset subs(&mat, {0,4,6});
	auto first = make_row_iterator(&subs, 0);
	auto last  = make_row_iterator(&subs, subs.rows());

	for(DenseMatrix::index_type i = 0; i < subs.rows(); ++i, ++first) {
		EXPECT_EQ(subs.row(i), *first);
	}

	EXPECT_TRUE(first == last);
}

TEST_F(DenseMatrixIteratorTest, index_column_iteration)
{
	auto mat = makeMatrix(10, 20);
	auto it = make_column_iterator(&mat, 0);

	for(DenseMatrix::index_type j = 0; j < mat.cols(); ++j) {
		EXPECT_EQ(mat.col(j), it[j]);
	}
}

TEST_F(DenseMatrixIteratorTest, index_row_iteration)
{
	auto mat = makeMatrix(10, 20);
	auto it = make_row_iterator(&mat, 0);

	for(DenseMatrix::index_type i = 0; i < mat.rows(); ++i) {
		EXPECT_EQ(mat.row(i), it[i]);
	}
}
