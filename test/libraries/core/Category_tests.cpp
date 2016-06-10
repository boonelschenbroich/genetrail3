#include <gtest/gtest.h>

#include <genetrail2/core/Category.h>
#include <genetrail2/core/EntityDatabase.h>

#include <string>
#include <vector>

using namespace GeneTrail;

TEST(Category, construct_string_iterator)
{
	EntityDatabase db;
	std::vector<std::string> names = {"A", "B", "C"};

	Category cat(&db, names.begin(), names.end());

	EXPECT_EQ(names.size(), cat.size());
	EXPECT_TRUE(cat.contains("A"));
	EXPECT_TRUE(cat.contains("B"));
	EXPECT_TRUE(cat.contains("C"));
}

TEST(Category, construct_string_iterator_duplicates)
{
	EntityDatabase db;
	std::vector<std::string> names = {"A", "B", "B", "C", "C"};

	Category cat(&db, names.begin(), names.end());

	EXPECT_EQ(3u, cat.size());
	EXPECT_TRUE(cat.contains("A"));
	EXPECT_TRUE(cat.contains("B"));
	EXPECT_TRUE(cat.contains("C"));
}

TEST(Category, construct_int_iterator)
{
	EntityDatabase db;
	std::vector<size_t> indices {
		db("A"), db("B"), db("C"), db("D")
	};

	Category cat(&db, indices.begin(), indices.end());

	EXPECT_EQ(indices.size(), cat.size());
	EXPECT_TRUE(cat.contains("A"));
	EXPECT_TRUE(cat.contains("B"));
	EXPECT_TRUE(cat.contains("C"));
	EXPECT_TRUE(cat.contains("D"));

	EXPECT_EQ("", cat.name());
	EXPECT_EQ("", cat.reference());
}

TEST(Category, construct_name)
{
	EntityDatabase db;
	Category cat(&db, "MyCategory");

	EXPECT_EQ(0u, cat.size());
	EXPECT_FALSE(cat.contains("Something"));

	EXPECT_EQ("MyCategory", cat.name());
}

TEST(Category, construct_entity_db)
{
	EntityDatabase db;
	Category cat(&db);

	EXPECT_EQ(0u, cat.size());
	EXPECT_FALSE(cat.contains("Something"));

	EXPECT_EQ("", cat.name());
	EXPECT_EQ("", cat.reference());
}

TEST(Category, insert)
{
	EntityDatabase db;
	Category cat(&db);

	EXPECT_EQ(0u, cat.size());
	EXPECT_FALSE(cat.contains("Something"));
	EXPECT_FALSE(cat.contains("SomethingElse"));
	cat.insert("Something");
	EXPECT_EQ(1u, cat.size());
	cat.insert("SomethingElse");
	EXPECT_EQ(2u, cat.size());

	EXPECT_TRUE(cat.contains("Something"));
	EXPECT_TRUE(cat.contains("SomethingElse"));
}

TEST(Category, operator_equals_1)
{
	EntityDatabase db;

	std::vector<std::string> names = {"X", "A", "haha", "Bla"};

	Category cat_a(&db, names.begin(), names.end());
	Category cat_b(&db, names.begin(), names.end());

	EXPECT_TRUE(cat_a == cat_b);

	cat_b.insert("Z");

	EXPECT_FALSE(cat_a == cat_b);
}

TEST(Category, operator_equals_2)
{
	EntityDatabase db;

	std::vector<std::string> names = {"X", "A", "haha", "Bla"};

	Category cat_a(&db, names.begin(), names.end());
	Category cat_b(&db, names.begin(), names.end());

	EXPECT_TRUE(cat_a == cat_b);

	cat_a.insert("Y");
	cat_b.insert("Z");

	EXPECT_FALSE(cat_a == cat_b);
}

TEST(Category, intersect)
{
	EntityDatabase db;

	std::vector<std::string> names = {"X", "A", "haha", "Bla"};

	Category cat_a(&db, names.begin(), names.end());
	Category cat_b(&db, names.begin(), names.end());

	cat_a.insert("Y");
	cat_b.insert("Z");

	auto result  = Category::intersect("new", cat_a, cat_b);
	auto result2 = Category::intersect("new", cat_b, cat_a);

	EXPECT_EQ(result, result2);

	EXPECT_EQ(4u, result.size());
	EXPECT_TRUE(result.contains("A"));
	EXPECT_TRUE(result.contains("Bla"));
	EXPECT_TRUE(result.contains("haha"));
	EXPECT_TRUE(result.contains("X"));
}

TEST(Category, intersect_throw)
{
	EntityDatabase db_a;
	EntityDatabase db_b;

	Category cat_a(&db_a);
	Category cat_b(&db_b);

	EXPECT_THROW(Category::intersect("new", cat_a, cat_b), std::invalid_argument);
}

TEST(Category, combine)
{
	EntityDatabase db;

	std::vector<std::string> names = {"X", "A", "haha", "Bla"};

	Category cat_a(&db, names.begin(), names.end());
	Category cat_b(&db, names.begin(), names.end());

	cat_a.insert("Y");
	cat_b.insert("Z");

	auto result  = Category::combine("new", cat_a, cat_b);
	auto result2 = Category::combine("new", cat_b, cat_a);

	EXPECT_EQ(result, result2);

	EXPECT_EQ(6u, result.size());
	EXPECT_TRUE(result.contains("A"));
	EXPECT_TRUE(result.contains("Bla"));
	EXPECT_TRUE(result.contains("haha"));
	EXPECT_TRUE(result.contains("X"));
	EXPECT_TRUE(result.contains("Y"));
	EXPECT_TRUE(result.contains("Z"));
}

TEST(Category, combine_throw)
{
	EntityDatabase db_a;
	EntityDatabase db_b;

	Category cat_a(&db_a);
	Category cat_b(&db_b);

	EXPECT_THROW(Category::combine("new", cat_a, cat_b), std::invalid_argument);
}

TEST(Category, operator_less)
{
	EntityDatabase db;

	std::vector<std::string> names = {"X", "A", "haha", "Bla"};

	Category cat_a(&db, names.begin(), names.end());
	Category cat_b(&db, names.begin(), names.end());

	EXPECT_FALSE(cat_a < cat_b);

	cat_b.insert("Z");

	EXPECT_TRUE(cat_a < cat_b);
}
