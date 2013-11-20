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
#include "../src/DenseMatrixWriter.h"
#include "../src/DenseMatrixReader.h"
#include <config.h>

#include <fstream>
#include <cstdlib>
#include <ctime>

#include <boost/filesystem.hpp>

using namespace GeneTrail;
namespace fs = boost::filesystem;

class DenseMatrixWriterTest : public ::testing::Test
{
	public:
		DenseMatrixWriterTest()
			: matrix45_rm_(TEST_DATA_PATH("binary_matrix4x5_rm.bmat")),
			  matrix45_cm_(TEST_DATA_PATH("binary_matrix4x5_cm.bmat")),
			  matrix45_ascii_(TEST_DATA_PATH("ascii_matrix4x5.mat")),
			  temp_file_name_(fs::unique_path().native())
		{
			srand(time(nullptr));
		}

		virtual void TearDown() override {
			// Delete the temporary file
			fs::remove(temp_file_name_);
		};

	protected:
		DenseMatrix buildKnownMatrix();
		DenseMatrix buildRandomMatrix();

		const std::string matrix45_rm_;
		const std::string matrix45_cm_;
		const std::string matrix45_ascii_;
		const std::string temp_file_name_;
};

DenseMatrix DenseMatrixWriterTest::buildKnownMatrix()
{
	DenseMatrix result(
		std::vector<std::string>{
			"row1",
			"row2",
			"row3",
			"row4"
		},
		std::vector<std::string>{
			"col1",
			"col2",
			"col3",
			"col4",
			"col5"
	});

	for(unsigned int j = 0; j < result.cols(); ++j) {
		for(unsigned int i = 0; i < result.rows(); ++i) {
			result(i,j) = (double)i * result.cols() + j;
		}
	}

	return result;
}

DenseMatrix DenseMatrixWriterTest::buildRandomMatrix()
{
	std::vector<std::string> row_names(100);
	std::vector<std::string> col_names(200);

	// Create random row names
	for(auto& s : row_names) {
		int length = 1 + rand() % 100;
		for(int j = 0; j < length; ++j) {
			s += (char)(65 + rand() % 26);
		}
	}

	// Create random column names
	for(auto& s : col_names) {
		int length = 1 + rand() % 100;
		for(int j = 0; j < length; ++j) {
			s += (char)(65 + rand() % 26);
		}
	}

	DenseMatrix out(row_names, col_names);
	out.matrix().setRandom();

	return out;
}

TEST_F(DenseMatrixWriterTest, writeText_known)
{
	DenseMatrix result = buildKnownMatrix();

	ASSERT_EQ(4, result.rows());
	ASSERT_EQ(5, result.cols());

	std::ifstream istrm(matrix45_ascii_);
	ASSERT_TRUE(istrm);

	// Read the correct file from file
	char in_buffer[1024];
	istrm.read(in_buffer, 1024);
	uint64_t bytes_read = istrm.gcount();

	std::ostringstream ostrm;
	ASSERT_TRUE(ostrm);

	// Write the test to the buffer
	DenseMatrixWriter writer;
	writer.writeText(ostrm, result);
	std::string tmp = ostrm.str();

	ASSERT_EQ(bytes_read, tmp.length());
	ASSERT_TRUE(memcmp(in_buffer, &tmp[0], bytes_read) == 0);
}

TEST_F(DenseMatrixWriterTest, binaryReadWrite_random)
{
	DenseMatrix out = buildRandomMatrix();

	std::ofstream ostrm(temp_file_name_, std::ios::binary);
	ASSERT_TRUE(ostrm);
	DenseMatrixWriter writer;
	writer.writeBinary(ostrm, out);
	ostrm.close();

	std::ifstream istrm(temp_file_name_, std::ios::binary);
	ASSERT_TRUE(istrm);
	DenseMatrixReader reader;
	DenseMatrix in = reader.read(istrm);
	istrm.close();

	ASSERT_EQ(out.rows(), in.rows());
	ASSERT_EQ(out.cols(), in.cols());

	for(unsigned int i = 0; i < in.rows(); ++i) {
		EXPECT_EQ(in.rowName(i), out.rowName(i));
	}

	for(unsigned int j = 0; j < in.cols(); ++j) {
		EXPECT_EQ(in.colName(j), out.colName(j));
	}

	for(unsigned int j = 0; j < in.cols(); ++j) {
		for(unsigned int i = 0; i < in.rows(); ++i) {
			EXPECT_EQ(in(i,j), out(i,j));
		}
	}
}

TEST_F(DenseMatrixWriterTest, binaryWrite_known)
{
	DenseMatrix result = buildKnownMatrix();

	std::ifstream istrm(matrix45_cm_, std::ios::binary);
	ASSERT_TRUE(istrm);

	// Read the correct file from file
	char in_buffer[1024];
	istrm.read(in_buffer, 1024);
	uint64_t bytes_read = istrm.gcount();

	std::ostringstream ostrm;
	ASSERT_TRUE(ostrm);

	// Write the test to the buffer
	DenseMatrixWriter writer;
	writer.writeBinary(ostrm, result);
	std::string tmp = ostrm.str();

	ASSERT_EQ(bytes_read, tmp.length());
	ASSERT_TRUE(memcmp(in_buffer, &tmp[0], bytes_read) == 0);
}

TEST_F(DenseMatrixWriterTest, textReadWrite_random)
{
	DenseMatrix out = buildRandomMatrix();

	std::ofstream ostrm(temp_file_name_);
	ASSERT_TRUE(ostrm);
	DenseMatrixWriter writer;
	writer.writeText(ostrm, out);
	ostrm.close();

	std::ifstream istrm(temp_file_name_);
	ASSERT_TRUE(istrm);
	DenseMatrixReader reader;
	DenseMatrix in = reader.read(istrm);
	istrm.close();

	ASSERT_EQ(out.rows(), in.rows());
	ASSERT_EQ(out.cols(), in.cols());

	for(unsigned int i = 0; i < in.rows(); ++i) {
		EXPECT_EQ(in.rowName(i), out.rowName(i));
	}

	for(unsigned int j = 0; j < in.cols(); ++j) {
		EXPECT_EQ(in.colName(j), out.colName(j));
	}

	for(unsigned int j = 0; j < in.cols(); ++j) {
		for(unsigned int i = 0; i < in.rows(); ++i) {
			EXPECT_NEAR(in(i,j), out(i,j), 0.00001);
		}
	}
}