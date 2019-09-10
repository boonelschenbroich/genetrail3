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

#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/Exception.h>
#include <genetrail2/core/NAHandler.h>
#include <config.h>

#include <fstream>

using namespace GeneTrail;


TEST(NAHandlerTest, noNANs){
	DenseMatrixReader reader;
	std::ifstream strm(TEST_DATA_PATH("naHandlerMatrix1.txt"), std::ios::binary);
	ASSERT_TRUE(strm.good());
	DenseMatrix m = reader.read(strm);
	
	NAHandler h;
	h.handle(m, NAStrategyRemove());
	
	ASSERT_EQ(3, m.rows());
	ASSERT_EQ(4, m.cols());

	EXPECT_EQ("G1", m.rowName(0));
	EXPECT_EQ("G2", m.rowName(1));
	EXPECT_EQ("G3", m.rowName(2));

	EXPECT_EQ("Hi", m.colName(0));
	EXPECT_EQ("Hi2", m.colName(1));
	EXPECT_EQ("Hi3", m.colName(2));
	EXPECT_EQ("Hi4", m.colName(3));

	EXPECT_EQ( 1.0, m(0,0));
	EXPECT_EQ( 2.0, m(0,1));
	EXPECT_EQ( 4.0, m(0,2));
	EXPECT_EQ( 3.0, m(0,3));
	EXPECT_EQ( 8.0, m(1,0));
	EXPECT_EQ( 7.0, m(1,1));
	EXPECT_EQ( 9.0, m(1,2));
	EXPECT_EQ( 4.0, m(1,3));
	EXPECT_EQ(10.0, m(2,0));
	EXPECT_EQ( 2.0, m(2,1));
	EXPECT_EQ(-3.0, m(2,2));
	EXPECT_EQ( 0.0, m(2,3));
}

TEST(NAHandlerTest, remove){
	DenseMatrixReader reader;
	std::ifstream strm(TEST_DATA_PATH("naHandlerMatrix2.txt"), std::ios::binary);
	ASSERT_TRUE(strm.good());
	DenseMatrix m = reader.read(strm);
	
	NAHandler h;
	h.handle(m, NAStrategyRemove());
	
	ASSERT_EQ(2, m.rows());
	ASSERT_EQ(4, m.cols());

	EXPECT_EQ("G1", m.rowName(0));
	EXPECT_EQ("G3", m.rowName(2));

	EXPECT_EQ("Hi", m.colName(0));
	EXPECT_EQ("Hi2", m.colName(1));
	EXPECT_EQ("Hi3", m.colName(2));
	EXPECT_EQ("Hi4", m.colName(3));

	EXPECT_EQ( 1.0, m(0,0));
	EXPECT_EQ( 2.0, m(0,1));
	EXPECT_EQ( 4.0, m(0,2));
	EXPECT_EQ( 3.0, m(0,3));
	EXPECT_EQ(10.0, m(1,0));
	EXPECT_EQ( 2.0, m(1,1));
	EXPECT_EQ(-3.0, m(1,2));
	EXPECT_EQ( 0.0, m(1,3));
}

TEST(NAHandlerTest, zero){
	DenseMatrixReader reader;
	std::ifstream strm(TEST_DATA_PATH("naHandlerMatrix3.txt"), std::ios::binary);
	ASSERT_TRUE(strm.good());
	DenseMatrix m = reader.read(strm);
	
	NAHandler h;
	h.handle(m, NAStrategyZero());
	
	ASSERT_EQ(3, m.rows());
	ASSERT_EQ(4, m.cols());

	EXPECT_EQ("G1", m.rowName(0));
	EXPECT_EQ("G2", m.rowName(1));
	EXPECT_EQ("G3", m.rowName(2));

	EXPECT_EQ("Hi", m.colName(0));
	EXPECT_EQ("Hi2", m.colName(1));
	EXPECT_EQ("Hi3", m.colName(2));
	EXPECT_EQ("Hi4", m.colName(3));

	EXPECT_EQ( 0.0, m(0,0));
	EXPECT_EQ( 0.0, m(0,1));
	EXPECT_EQ( 4.0, m(0,2));
	EXPECT_EQ( 3.0, m(0,3));
	EXPECT_EQ( 8.0, m(1,0));
	EXPECT_EQ( 0.0, m(1,1));
	EXPECT_EQ( 9.0, m(1,2));
	EXPECT_EQ( 4.0, m(1,3));
	EXPECT_EQ( 0.0, m(2,0));
	EXPECT_EQ( 2.0, m(2,1));
	EXPECT_EQ( 0.0, m(2,2));
	EXPECT_EQ( 0.0, m(2,3));
}

TEST(NAHandlerTest, mean){
	DenseMatrixReader reader;
	std::ifstream strm(TEST_DATA_PATH("naHandlerMatrix3.txt"), std::ios::binary);
	ASSERT_TRUE(strm.good());
	DenseMatrix m = reader.read(strm);
	
	NAHandler h;
	h.handle(m, NAStrategyMean());
	
	ASSERT_EQ(3, m.rows());
	ASSERT_EQ(4, m.cols());

	EXPECT_EQ("G1", m.rowName(0));
	EXPECT_EQ("G2", m.rowName(1));
	EXPECT_EQ("G3", m.rowName(2));

	EXPECT_EQ("Hi", m.colName(0));
	EXPECT_EQ("Hi2", m.colName(1));
	EXPECT_EQ("Hi3", m.colName(2));
	EXPECT_EQ("Hi4", m.colName(3));

	EXPECT_EQ( 3.5, m(0,0));
	EXPECT_EQ( 3.5, m(0,1));
	EXPECT_EQ( 4.0, m(0,2));
	EXPECT_EQ( 3.0, m(0,3));
	EXPECT_EQ( 8.0, m(1,0));
	EXPECT_EQ( 7.0, m(1,1));
	EXPECT_EQ( 9.0, m(1,2));
	EXPECT_EQ( 4.0, m(1,3));
	EXPECT_EQ( 1.0, m(2,0));
	EXPECT_EQ( 2.0, m(2,1));
	EXPECT_EQ( 1.0, m(2,2));
	EXPECT_EQ( 0.0, m(2,3));
}

TEST(NAHandlerTest, median){
	DenseMatrixReader reader;
	std::ifstream strm(TEST_DATA_PATH("naHandlerMatrix3.txt"), std::ios::binary);
	ASSERT_TRUE(strm.good());
	DenseMatrix m = reader.read(strm);
	
	NAHandler h;
	h.handle(m, NAStrategyMedian());
	
	ASSERT_EQ(3, m.rows());
	ASSERT_EQ(4, m.cols());

	EXPECT_EQ("G1", m.rowName(0));
	EXPECT_EQ("G2", m.rowName(1));
	EXPECT_EQ("G3", m.rowName(2));

	EXPECT_EQ("Hi", m.colName(0));
	EXPECT_EQ("Hi2", m.colName(1));
	EXPECT_EQ("Hi3", m.colName(2));
	EXPECT_EQ("Hi4", m.colName(3));

	EXPECT_EQ( 3.5, m(0,0));
	EXPECT_EQ( 3.5, m(0,1));
	EXPECT_EQ( 4.0, m(0,2));
	EXPECT_EQ( 3.0, m(0,3));
	EXPECT_EQ( 8.0, m(1,0));
	EXPECT_EQ( 8.0, m(1,1));
	EXPECT_EQ( 9.0, m(1,2));
	EXPECT_EQ( 4.0, m(1,3));
	EXPECT_EQ( 1.0, m(2,0));
	EXPECT_EQ( 2.0, m(2,1));
	EXPECT_EQ( 1.0, m(2,2));
	EXPECT_EQ( 0.0, m(2,3));
}


