#include <gtest/gtest.h>

#include <genetrail2/core/Category.h>
#include <genetrail2/core/CategoryDatabase.h>
#include <genetrail2/core/JsonCategoryFile.h>

#include <config.h>

#include <boost/filesystem.hpp>

using namespace GeneTrail;
namespace fs = boost::filesystem;

class JsonCategoryFileTest : public ::testing::Test
{
  public:
	JsonCategoryFileTest()
	    : tmp_file_name_("/tmp/" + fs::unique_path().native())
	{
	}

	void TearDown() override { fs::remove(tmp_file_name_); }

  protected:
	std::string tmp_file_name_;
};

TEST_F(JsonCategoryFileTest, read)
{
	auto db = std::make_shared<EntityDatabase>();
	JsonCategoryFile f(db, TEST_DATA_PATH("Categories.json"));

	ASSERT_TRUE(f);

	auto database = f.read();

	EXPECT_FALSE(f);

	EXPECT_EQ("Test Database", database.name());
	EXPECT_EQ("Jon Doe", database.editor().name);
	EXPECT_EQ("jon@doe.org", database.editor().email);
	EXPECT_EQ("http://doe.org", database.sourceUrl());
	EXPECT_EQ("GeneSymbol", database.identifier());

	ASSERT_EQ(2, database.size());

	EXPECT_EQ("Test Category 1", database[0].name());
	EXPECT_EQ("http://doe.org/1", database[0].reference());
	EXPECT_EQ(4, database[0].size());
	EXPECT_TRUE(database[0].contains("TP53"));
	EXPECT_TRUE(database[0].contains("MDM2"));
	EXPECT_TRUE(database[0].contains("TRIM71"));
	EXPECT_TRUE(database[0].contains("TCF3"));

	EXPECT_EQ("Test Category 2", database[1].name());
	EXPECT_EQ("http://doe.org/2", database[1].reference());
	EXPECT_EQ(2, database[1].size());
	EXPECT_TRUE(database[1].contains("ABC"));
	EXPECT_TRUE(database[1].contains("OR1021"));
}

TEST_F(JsonCategoryFileTest, readMetadata)
{
	auto db = std::make_shared<EntityDatabase>();

	JsonCategoryFile f(db, TEST_DATA_PATH("CategoryMetadata.json"));

	ASSERT_TRUE(f);

	auto database = f.read();

	EXPECT_FALSE(f);

	EXPECT_EQ("Metadata", database.name());
	EXPECT_EQ("Jane Doe", database.editor().name);
	EXPECT_EQ("jane@doe.net", database.editor().email);
	EXPECT_EQ("http://doe.net", database.sourceUrl());
	EXPECT_EQ("EntrezGene", database.identifier());

	const Metadata& md = database.metadata();

	EXPECT_EQ(3, md.size());
	EXPECT_TRUE(md.has("purpose"));
	EXPECT_TRUE(md.has("stuff"));
	EXPECT_TRUE(md.has("objTest"));

	EXPECT_EQ("Unit tests", get<std::string>(md.get("purpose")));

	const auto& array = get<Metadata::Array>(md.get("stuff"));
	EXPECT_EQ(7, array.size());
	EXPECT_EQ(1, get<int64_t>(array[0]));
	EXPECT_EQ(2, get<int64_t>(array[1]));
	EXPECT_EQ(3, get<int64_t>(array[2]));
	EXPECT_EQ(4, get<int64_t>(array[3]));
	EXPECT_EQ(true, get<bool>(array[4]));
	EXPECT_EQ("", get<std::string>(array[5]));
	EXPECT_EQ(nullptr, get<std::nullptr_t>(array[6]));

	const auto& obj = get<Metadata::Object>(md.get("objTest"));
	EXPECT_EQ(3, obj.size());
	EXPECT_NE(obj.end(), obj.find("asdasd"));
	EXPECT_EQ(1.0, get<double>(obj.find("asdasd")->second));
	EXPECT_NE(obj.end(), obj.find("blabla"));
	EXPECT_EQ(true, get<bool>(obj.find("blabla")->second));
	EXPECT_NE(obj.end(), obj.find("nested"));

	const auto& nested = get<Metadata::Array>(obj.find("nested")->second);
	EXPECT_EQ(3, nested.size());
	EXPECT_EQ(1, get<int64_t>(nested[0]));
	EXPECT_EQ(2, get<int64_t>(nested[1]));
	EXPECT_EQ(3, get<int64_t>(nested[2]));

	ASSERT_EQ(1, database.size());
	const auto& cat = database[0];

	EXPECT_EQ("CatA", cat.name());
	EXPECT_EQ("http://doe.net/CatA", cat.reference());
	EXPECT_EQ(3, cat.size());
	EXPECT_EQ(1, cat.metadata().size());
	EXPECT_TRUE(cat.metadata().has("concentration"));
	EXPECT_EQ(0.4, get<double>(cat.metadata().get("concentration")));
}

TEST_F(JsonCategoryFileTest, readWrite)
{
	auto db = std::make_shared<EntityDatabase>();
	JsonCategoryFile in(db, TEST_DATA_PATH("Categories.json"));
	ASSERT_TRUE(in);

	auto database = in.read();

	EXPECT_FALSE(in);

	{
		JsonCategoryFile out(db, tmp_file_name_, FileOpenMode::WRITE);
		ASSERT_TRUE(out.write(database));
	}

	JsonCategoryFile in2(db, tmp_file_name_);
	ASSERT_TRUE(in2);

	auto database2 = in2.read();

	EXPECT_FALSE(in2);

	EXPECT_EQ(database.name(), database2.name());
	EXPECT_EQ(database.identifier(), database2.identifier());
	EXPECT_EQ(database.sourceUrl(), database2.sourceUrl());
	EXPECT_EQ(database.editor().name, database2.editor().name);
	EXPECT_EQ(database.editor().email, database2.editor().email);
	EXPECT_EQ(database.size(), database2.size());

	EXPECT_EQ(database[0].size(), database2[0].size());
	EXPECT_EQ(database[0].name(), database2[0].name());
	EXPECT_EQ(database[0].reference(), database2[0].reference());

	EXPECT_EQ(database[1].size(), database2[1].size());
	EXPECT_EQ(database[1].name(), database2[1].name());
	EXPECT_EQ(database[1].reference(), database2[1].reference());
}

TEST_F(JsonCategoryFileTest, readWriteMetadata)
{
	auto db = std::make_shared<EntityDatabase>();
	JsonCategoryFile in(db, TEST_DATA_PATH("CategoryMetadata.json"));

	ASSERT_TRUE(in);

	auto database = in.read();

	EXPECT_FALSE(in);

	{
		JsonCategoryFile out(db, tmp_file_name_, FileOpenMode::WRITE);
		ASSERT_TRUE(out.write(database));
	}

	JsonCategoryFile in2(db, tmp_file_name_);
	ASSERT_TRUE(in2);

	auto database2 = in2.read();

	EXPECT_FALSE(in2);

	EXPECT_EQ(database.name(), database2.name());
	EXPECT_EQ(database.identifier(), database2.identifier());
	EXPECT_EQ(database.sourceUrl(), database2.sourceUrl());
	EXPECT_EQ(database.editor().name, database2.editor().name);
	EXPECT_EQ(database.editor().email, database2.editor().email);
	EXPECT_EQ(database.size(), database2.size());

	const auto& md1 = database.metadata();
	const auto& md2 = database2.metadata();

	EXPECT_EQ(md1.size(), md2.size());

	EXPECT_TRUE(md2.has("purpose"));
	EXPECT_TRUE(md2.has("stuff"));
	EXPECT_TRUE(md2.has("objTest"));

	EXPECT_EQ("Unit tests", get<std::string>(md2.get("purpose")));

	const auto& array = get<Metadata::Array>(md2.get("stuff"));
	EXPECT_EQ(7, array.size());

	const auto& obj = get<Metadata::Object>(md2.get("objTest"));
	EXPECT_EQ(3, obj.size());
}
