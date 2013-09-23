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
#include "../src/DenseMatrixReader.h"
#include <config.h>

#include <fstream>

using namespace GeneTrail;

class DenseMatrixReaderTest : public ::testing::Test
{
	public:
		DenseMatrixReaderTest()
			: matrix45_rm_(TEST_DATA_PATH("binary_matrix4x5_rm.bmat")),
			  matrix_name_file(TEST_DATA_PATH("matrix_names.txt")),
			  matrix_name_additional_column_file(TEST_DATA_PATH("matrix_names_additional_column.txt")),
			  matrix_noname_file(TEST_DATA_PATH("matrix_nonames.txt"))
		{
		}

	protected:
		const std::string matrix45_rm_;
		const std::string matrix_name_file;
		const std::string matrix_name_additional_column_file;
		const std::string matrix_noname_file;
};


TEST_F(DenseMatrixReaderTest, binaryRead)
{
	DenseMatrixReader reader;

	std::ifstream strm(matrix45_rm_, std::ios::binary);

	ASSERT_TRUE(strm);

	DenseMatrix result = reader.read(strm);

	ASSERT_EQ(4, result.rows());
	ASSERT_EQ(5, result.cols());

	EXPECT_EQ("row1", result.rowName(0));
	EXPECT_EQ("row2", result.rowName(1));
	EXPECT_EQ("row3", result.rowName(2));
	EXPECT_EQ("row4", result.rowName(3));

	EXPECT_EQ("col1", result.colName(0));
	EXPECT_EQ("col2", result.colName(1));
	EXPECT_EQ("col3", result.colName(2));
	EXPECT_EQ("col4", result.colName(3));
	EXPECT_EQ("col5", result.colName(4));

	EXPECT_EQ( 1.0, result(0,0));
	EXPECT_EQ( 2.0, result(0,1));
	EXPECT_EQ( 3.0, result(0,2));
	EXPECT_EQ( 4.0, result(0,3));
	EXPECT_EQ( 5.0, result(0,4));
	EXPECT_EQ( 6.0, result(1,0));
	EXPECT_EQ( 7.0, result(1,1));
	EXPECT_EQ( 8.0, result(1,2));
	EXPECT_EQ( 9.0, result(1,3));
	EXPECT_EQ(10.0, result(1,4));
	EXPECT_EQ(11.0, result(2,0));
	EXPECT_EQ(12.0, result(2,1));
	EXPECT_EQ(13.0, result(2,2));
	EXPECT_EQ(14.0, result(2,3));
	EXPECT_EQ(15.0, result(2,4));
	EXPECT_EQ(16.0, result(3,0));
	EXPECT_EQ(17.0, result(3,1));
	EXPECT_EQ(18.0, result(3,2));
	EXPECT_EQ(19.0, result(3,3));
	EXPECT_EQ(20.0, result(3,4));
}

// Test that calls to the constructor that contains invalid
// data is handled properly
TEST_F(DenseMatrixReaderTest, read_invalid_stream)
{
	std::ifstream strm("does_not_exist.txt");

	ASSERT_FALSE(strm.good());

	DenseMatrixReader reader;
	DenseMatrix mat = reader.read(strm);
	//TODO fix this
//	EXPECT_FALSE(reader.readMatrix(strm));

	EXPECT_EQ(-1, mat.rowIndex("blablubb"));
//	EXPECT_EQ("", mat.rowName(0));
	EXPECT_EQ(DenseMatrix::Matrix(), mat.matrix());
}

TEST_F(DenseMatrixReaderTest, read_valid_names)
{
	std::ifstream strm(matrix_name_file);

	ASSERT_TRUE(strm);

	DenseMatrixReader reader;
	DenseMatrix matrix = reader.read(strm, DenseMatrixReader::READ_ROW_NAMES);

	ASSERT_EQ(5, matrix.rows());
	ASSERT_EQ(3, matrix.cols());

	EXPECT_EQ(0, matrix.rowIndex("blablubb1"));
	EXPECT_EQ(1, matrix.rowIndex("blablubb2"));
	EXPECT_EQ(2, matrix.rowIndex("blablubb3"));
	EXPECT_EQ(3, matrix.rowIndex("blablubb4"));
	EXPECT_EQ(4, matrix.rowIndex("blablubb5"));

	EXPECT_EQ(-1, matrix.rowIndex("hossa"));
	// Make sure that we did not accidentally introduce an index
	EXPECT_EQ(-1, matrix.rowIndex("hossa"));

	EXPECT_EQ("blablubb1", matrix.rowName(0));
	EXPECT_EQ("blablubb2", matrix.rowName(1));
	EXPECT_EQ("blablubb3", matrix.rowName(2));
	EXPECT_EQ("blablubb4", matrix.rowName(3));
	EXPECT_EQ("blablubb5", matrix.rowName(4));

	EXPECT_EQ( 1.0, matrix(0, 0));
	EXPECT_EQ( 2.0, matrix(0, 1));
	EXPECT_EQ( 3.0, matrix(0, 2));
	EXPECT_EQ( 4.0, matrix(1, 0));
	EXPECT_EQ( 5.0, matrix(1, 1));
	EXPECT_EQ( 6.0, matrix(1, 2));
	EXPECT_EQ( 7.0, matrix(2, 0));
	EXPECT_EQ( 8.0, matrix(2, 1));
	EXPECT_EQ( 9.0, matrix(2, 2));
	EXPECT_EQ(10.0, matrix(3, 0));
	EXPECT_EQ(11.0, matrix(3, 1));
	EXPECT_EQ(12.0, matrix(3, 2));
	EXPECT_EQ(13.0, matrix(4, 0));
	EXPECT_EQ(14.0, matrix(4, 1));
	EXPECT_EQ(15.0, matrix(4, 2));
}

TEST_F(DenseMatrixReaderTest, read_valid_names_transposed)
{
	std::ifstream strm(matrix_name_file);

	ASSERT_TRUE(strm);

	DenseMatrixReader reader;
	DenseMatrix matrix = reader.read(strm, DenseMatrixReader::READ_ROW_NAMES | DenseMatrixReader::TRANSPOSE);

	ASSERT_EQ(3, matrix.rows());
	ASSERT_EQ(5, matrix.cols());

	EXPECT_EQ(0, matrix.colIndex("blablubb1"));
	EXPECT_EQ(1, matrix.colIndex("blablubb2"));
	EXPECT_EQ(2, matrix.colIndex("blablubb3"));
	EXPECT_EQ(3, matrix.colIndex("blablubb4"));
	EXPECT_EQ(4, matrix.colIndex("blablubb5"));

	EXPECT_EQ(-1, matrix.colIndex("hossa"));
	// Make sure that we did not accidentally introduce an index
	EXPECT_EQ(-1, matrix.colIndex("hossa"));

	EXPECT_EQ("blablubb1", matrix.colName(0));
	EXPECT_EQ("blablubb2", matrix.colName(1));
	EXPECT_EQ("blablubb3", matrix.colName(2));
	EXPECT_EQ("blablubb4", matrix.colName(3));
	EXPECT_EQ("blablubb5", matrix.colName(4));

	EXPECT_EQ( 1.0, matrix(0, 0));
	EXPECT_EQ( 2.0, matrix(1, 0));
	EXPECT_EQ( 3.0, matrix(2, 0));
	EXPECT_EQ( 4.0, matrix(0, 1));
	EXPECT_EQ( 5.0, matrix(1, 1));
	EXPECT_EQ( 6.0, matrix(2, 1));
	EXPECT_EQ( 7.0, matrix(0, 2));
	EXPECT_EQ( 8.0, matrix(1, 2));
	EXPECT_EQ( 9.0, matrix(2, 2));
	EXPECT_EQ(10.0, matrix(0, 3));
	EXPECT_EQ(11.0, matrix(1, 3));
	EXPECT_EQ(12.0, matrix(2, 3));
	EXPECT_EQ(13.0, matrix(0, 4));
	EXPECT_EQ(14.0, matrix(1, 4));
	EXPECT_EQ(15.0, matrix(2, 4));
}

TEST_F(DenseMatrixReaderTest, read_valid_names_additional_column)
{
	std::ifstream strm(matrix_name_additional_column_file);

	ASSERT_TRUE(strm);

	DenseMatrixReader reader;
	DenseMatrix matrix = reader.read(strm, DenseMatrixReader::READ_COL_NAMES | DenseMatrixReader::READ_ROW_NAMES | DenseMatrixReader::ADDITIONAL_COL_NAME);

	ASSERT_EQ(5, matrix.rows());
	ASSERT_EQ(3, matrix.cols());

	EXPECT_EQ(0, matrix.rowIndex("blablubb1"));
	EXPECT_EQ(1, matrix.rowIndex("blablubb2"));
	EXPECT_EQ(2, matrix.rowIndex("blablubb3"));
	EXPECT_EQ(3, matrix.rowIndex("blablubb4"));
	EXPECT_EQ(4, matrix.rowIndex("blablubb5"));

	EXPECT_EQ(-1, matrix.rowIndex("hossa"));
	// Make sure that we did not accidentally introduce an index
	EXPECT_EQ(-1, matrix.rowIndex("hossa"));

	EXPECT_EQ("blablubb1", matrix.rowName(0));
	EXPECT_EQ("blablubb2", matrix.rowName(1));
	EXPECT_EQ("blablubb3", matrix.rowName(2));
	EXPECT_EQ("blablubb4", matrix.rowName(3));
	EXPECT_EQ("blablubb5", matrix.rowName(4));

	EXPECT_EQ("A", matrix.colName(0));
	EXPECT_EQ("B", matrix.colName(1));
	EXPECT_EQ("C", matrix.colName(2));

	EXPECT_EQ(0, matrix.colIndex("A"));
	EXPECT_EQ(1, matrix.colIndex("B"));
	EXPECT_EQ(2, matrix.colIndex("C"));

	EXPECT_EQ( 1.0, matrix(0, 0));
	EXPECT_EQ( 2.0, matrix(0, 1));
	EXPECT_EQ( 3.0, matrix(0, 2));
	EXPECT_EQ( 4.0, matrix(1, 0));
	EXPECT_EQ( 5.0, matrix(1, 1));
	EXPECT_EQ( 6.0, matrix(1, 2));
	EXPECT_EQ( 7.0, matrix(2, 0));
	EXPECT_EQ( 8.0, matrix(2, 1));
	EXPECT_EQ( 9.0, matrix(2, 2));
	EXPECT_EQ(10.0, matrix(3, 0));
	EXPECT_EQ(11.0, matrix(3, 1));
	EXPECT_EQ(12.0, matrix(3, 2));
	EXPECT_EQ(13.0, matrix(4, 0));
	EXPECT_EQ(14.0, matrix(4, 1));
	EXPECT_EQ(15.0, matrix(4, 2));
}

TEST_F(DenseMatrixReaderTest, read_valid_nonames)
{
	std::ifstream strm(matrix_noname_file);

	ASSERT_TRUE(strm);

	DenseMatrixReader reader;
	DenseMatrix matrix = reader.read(strm, DenseMatrixReader::NO_OPTIONS);

	ASSERT_EQ(3, matrix.rows());
	ASSERT_EQ(4, matrix.cols());

	EXPECT_EQ(-1, matrix.rowIndex("hossa"));
	// Make sure that we did not accidentally introduce an index
	EXPECT_EQ(-1, matrix.rowIndex("hossa"));

	EXPECT_EQ("", matrix.rowName(0));
	EXPECT_EQ("", matrix.rowName(1));
	EXPECT_EQ("", matrix.rowName(2));

	EXPECT_EQ( 1.0, matrix(0, 0));
	EXPECT_EQ( 2.0, matrix(0, 1));
	EXPECT_EQ( 3.0, matrix(0, 2));
	EXPECT_EQ( 4.0, matrix(0, 3));
	EXPECT_EQ( 5.0, matrix(1, 0));
	EXPECT_EQ( 6.0, matrix(1, 1));
	EXPECT_EQ( 7.0, matrix(1, 2));
	EXPECT_EQ( 8.0, matrix(1, 3));
	EXPECT_EQ( 9.0, matrix(2, 0));
	EXPECT_EQ(10.0, matrix(2, 1));
	EXPECT_EQ(11.0, matrix(2, 2));
	EXPECT_EQ(12.0, matrix(2, 3));
}
