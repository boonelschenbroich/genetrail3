/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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

#include <genetrail2/core/Category.h>
#include <genetrail2/core/Scores.h>

using namespace GeneTrail;

class ScoresTest : public ::testing::Test
{
};

TEST_F(ScoresTest, Constructor)
{
	Scores scores;

	EXPECT_EQ(size_t(0), scores.size());
}

TEST_F(ScoresTest, Iteration)
{
	Scores scores;

	scores.emplace_back("A", 1.2);
	scores.emplace_back("D", 1.1);
	scores.emplace_back("G", 1.0);
	scores.emplace_back("B", -1.2);
	scores.emplace_back("C", 1.0);

	ASSERT_EQ(size_t(5), scores.size());

	auto begin = scores.begin();

	EXPECT_EQ("A", begin->name());
	EXPECT_EQ(1.2, begin->score());
	++begin;

	EXPECT_EQ("D", begin->name());
	EXPECT_EQ(1.1, begin->score());
	++begin;

	EXPECT_EQ("G", begin->name());
	EXPECT_EQ(1.0, begin->score());
	++begin;

	EXPECT_EQ("B", begin->name());
	EXPECT_EQ(-1.2, begin->score());
	++begin;

	EXPECT_EQ("C", begin->name());
	EXPECT_EQ(1.0, begin->score());
	++begin;

	EXPECT_EQ(scores.end(), begin);
}

TEST_F(ScoresTest, NameIteration)
{
	Scores scores;

	scores.emplace_back("A", 1.2);
	scores.emplace_back("D", 1.1);
	scores.emplace_back("G", 1.0);
	scores.emplace_back("B", -1.2);
	scores.emplace_back("C", 1.0);

	ASSERT_EQ(size_t(5), scores.size());

	auto begin = scores.names().begin();

	EXPECT_EQ("A", *begin);
	++begin;

	EXPECT_EQ("D", *begin);
	++begin;

	EXPECT_EQ("G", *begin);
	++begin;

	EXPECT_EQ("B", *begin);
	++begin;

	EXPECT_EQ("C", *begin);
	++begin;

	EXPECT_EQ(scores.names().end(), begin);
}

TEST_F(ScoresTest, ScoreIteration)
{
	Scores scores;

	scores.emplace_back("A", 1.2);
	scores.emplace_back("D", 1.1);
	scores.emplace_back("G", 1.0);
	scores.emplace_back("B", -1.2);
	scores.emplace_back("C", 1.0);

	ASSERT_EQ(size_t(5), scores.size());

	auto begin = scores.scores().begin();

	EXPECT_EQ(1.2, *begin);
	++begin;

	EXPECT_EQ(1.1, *begin);
	++begin;

	EXPECT_EQ(1.0, *begin);
	++begin;

	EXPECT_EQ(-1.2, *begin);
	++begin;

	EXPECT_EQ(1.0, *begin);
	++begin;

	EXPECT_EQ(scores.scores().end(), begin);
}

TEST_F(ScoresTest, sortByNameInvariant)
{
	Scores scores;
	EXPECT_TRUE(scores.isSortedByName());

	scores.emplace_back("AA", 1.0);
	ASSERT_TRUE(scores.isSortedByName());
	scores.emplace_back("B", 1.0);
	ASSERT_TRUE(scores.isSortedByName());
	scores.emplace_back("C", 1.0);
	ASSERT_TRUE(scores.isSortedByName());
	scores.emplace_back("E", 1.0);
	ASSERT_TRUE(scores.isSortedByName());
	scores.emplace_back("AB", 1.0);
	ASSERT_FALSE(scores.isSortedByName());

	scores.sortByName();
	ASSERT_TRUE(scores.isSortedByName());
}

TEST_F(ScoresTest, subsetSorted)
{
	Category c("");

	c.insert("A");
	c.insert("G");
	c.insert("H");

	Scores scores;
	scores.emplace_back("A", 1.1);
	scores.emplace_back("B", 1.1);
	scores.emplace_back("G", 1.1);
	EXPECT_TRUE(scores.isSortedByName());

	Scores subset = scores.subset(c);

	EXPECT_EQ(size_t(2), subset.size());
	EXPECT_TRUE(subset.contains("A"));
	EXPECT_TRUE(subset.contains("G"));
	EXPECT_FALSE(subset.contains("H"));
	EXPECT_FALSE(subset.contains("B"));
}

TEST_F(ScoresTest, subsetUnsorted)
{
	Category c("");

	c.insert("P");
	c.insert("A");
	c.insert("G");
	c.insert("H");

	Scores scores;
	scores.emplace_back("F", 1.1);
	scores.emplace_back("A", 1.1);
	scores.emplace_back("P", 1.1);
	scores.emplace_back("B", 1.1);
	scores.emplace_back("G", 1.1);
	EXPECT_FALSE(scores.isSortedByName());

	Scores subset = scores.subset(c);

	EXPECT_EQ(size_t(3), subset.size());
	EXPECT_TRUE(subset.contains("A"));
	EXPECT_TRUE(subset.contains("G"));
	EXPECT_TRUE(subset.contains("P"));
	EXPECT_FALSE(subset.contains("H"));
	EXPECT_FALSE(subset.contains("B"));
}
