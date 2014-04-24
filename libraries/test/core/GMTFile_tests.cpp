#include <gtest/gtest.h>

#include <genetrail2/core/Category.h>
#include <genetrail2/core/GMTFile.h>

#include <config.h>

using namespace GeneTrail;

TEST(GMTFile, read)
{
	GMTFile f(TEST_DATA_PATH("categories.gmt"));

	ASSERT_TRUE(f);

	std::vector<Category> categories;

	while(f) {
		categories.push_back(f.read());
	}

	EXPECT_EQ(5, categories.size());

	EXPECT_EQ("CatA", categories[0].name());
	EXPECT_EQ("CatB", categories[1].name());
	EXPECT_EQ("Cat C", categories[2].name());
	EXPECT_EQ("Cat D", categories[3].name());
	EXPECT_EQ("Crazy Cat E", categories[4].name());

	EXPECT_EQ("http://somthing.de", categories[0].reference());
	EXPECT_EQ("Strange but valid", categories[1].reference());
	EXPECT_EQ("https://aaa.xxx", categories[2].reference());
	EXPECT_EQ("", categories[3].reference());
	EXPECT_EQ("\" \"", categories[4].reference());

	EXPECT_EQ(3, categories[0].size());
	EXPECT_FALSE(categories[0].empty());
	EXPECT_TRUE(categories[0].contains("A"));
	EXPECT_TRUE(categories[0].contains("B"));
	EXPECT_TRUE(categories[0].contains("Bla Bla"));

	EXPECT_EQ(2, categories[1].size());
	EXPECT_FALSE(categories[1].empty());
	EXPECT_TRUE(categories[1].contains("123"));
	EXPECT_TRUE(categories[1].contains("323"));

	EXPECT_EQ(2, categories[2].size());
	EXPECT_FALSE(categories[2].empty());
	EXPECT_TRUE(categories[2].contains("A"));
	EXPECT_TRUE(categories[2].contains("B"));

	EXPECT_EQ(2, categories[3].size());
	EXPECT_FALSE(categories[3].empty());
	EXPECT_TRUE(categories[3].contains("A"));
	EXPECT_TRUE(categories[3].contains("B"));

	EXPECT_EQ(0, categories[4].size());
	EXPECT_TRUE(categories[4].empty());
}

