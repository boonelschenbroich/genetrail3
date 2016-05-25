/*
 * GeneTrail2 - An efficient library for interpreting genetic data
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

#include <genetrail2/core/Metadata.h>

using namespace GeneTrail;

TEST(MetadataTest, addNull)
{
	Metadata metadata;

	metadata.add("Value", nullptr);

	std::nullptr_t tmp;
	ASSERT_NO_THROW(tmp = get<std::nullptr_t>(metadata.get("Value")));
	EXPECT_EQ(nullptr, tmp);
}

TEST(MetadataTest, addString)
{
	Metadata metadata;

	// Test that we can add std::strings
	metadata.add("Value", std::string("blubb"));
	std::string tmp;
	ASSERT_NO_THROW(tmp = get<std::string>(metadata.get("Value")));
	EXPECT_EQ("blubb", tmp);

	// Test that the conversion from const char* is working
	metadata.add("Value2", "bla");
	ASSERT_NO_THROW(tmp = get<std::string>(metadata.get("Value2")));
	EXPECT_EQ("bla", tmp);
}

TEST(MetadataTest, addInt)
{
	Metadata metadata;

	// Check that the overloads are working
	metadata.add("Value1", static_cast<uint8_t>(128));
	metadata.add("Value2", static_cast<uint16_t>(34343));
	metadata.add("Value3", static_cast<uint32_t>(2315813));
	metadata.add("Value4", static_cast<uint64_t>(112232312335813));
	metadata.add("Value", 11235813);

	int64_t tmp;
	ASSERT_NO_THROW(tmp = get<int64_t>(metadata.get("Value")));
	EXPECT_EQ(11235813, tmp);
}

TEST(MetadataTest, addBool)
{
	Metadata metadata;

	metadata.add("Value", true);
	bool tmp;
	ASSERT_NO_THROW(tmp = get<bool>(metadata.get("Value")));
	EXPECT_TRUE(tmp);

	metadata.add("Value2", false);
	ASSERT_NO_THROW(tmp = get<bool>(metadata.get("Value2")));
	EXPECT_FALSE(tmp);
}

TEST(MetadataTest, addDouble)
{
	Metadata metadata;

	double tmp;
	metadata.add("Value", 0.1);
	ASSERT_NO_THROW(tmp = get<double>(metadata.get("Value")));
	EXPECT_EQ(0.1, tmp);

	metadata.add("Value2", -1.2f);
	ASSERT_NO_THROW(tmp = get<double>(metadata.get("Value2")));
	EXPECT_NEAR(-1.2, tmp, 1e-7);
}

TEST(MetadataTest, addArray)
{
	Metadata metadata;

	Metadata::Array values;
	values.emplace_back("Hallo");
	values.emplace_back(1.0);
	values.emplace_back(true);

	metadata.add("Value", values);

	Metadata::Array tmp;
	ASSERT_NO_THROW(tmp = get<Metadata::Array>(metadata.get("Value")));

	ASSERT_EQ(3, tmp.size());

	std::string tmp_string;
	ASSERT_NO_THROW(tmp_string = get<std::string>(tmp[0]));
	EXPECT_EQ("Hallo", tmp_string);

	double tmp_double;
	ASSERT_NO_THROW(tmp_double = get<double>(tmp[1]));
	EXPECT_EQ(1.0, tmp_double);

	bool tmp_bool;
	ASSERT_NO_THROW(tmp_bool = get<bool>(tmp[2]));
	EXPECT_EQ(true, tmp_bool);
}

TEST(MetadataTest, addObject)
{
	Metadata metadata;

	Metadata::Object values;
	values.emplace("Hallo", Metadata::Value("Hiho"));
	values.emplace("Bla", Metadata::Value(1.0));
	values.emplace("Blubber", Metadata::Value(true));

	metadata.add("Value", values);

	Metadata::Object tmp;
	ASSERT_NO_THROW(tmp = get<Metadata::Object>(metadata.get("Value")));

	ASSERT_EQ(3, tmp.size());

	{
	std::string tmp_string;
	auto res = tmp.find("Hallo");
	ASSERT_NE(tmp.end(), res);
	ASSERT_NO_THROW(tmp_string = get<std::string>(res->second));
	EXPECT_EQ("Hiho", tmp_string);
	}

	{
	double tmp_double;
	auto res = tmp.find("Bla");
	ASSERT_NE(tmp.end(), res);
	ASSERT_NO_THROW(tmp_double = get<double>(res->second));
	EXPECT_EQ(1.0, tmp_double);
	}

	{
	bool tmp_bool;
	auto res = tmp.find("Blubber");
	ASSERT_NE(tmp.end(), res);
	ASSERT_NO_THROW(tmp_bool = get<bool>(res->second));
	EXPECT_EQ(true, tmp_bool);
	}
}

TEST(MetadataTest, remove)
{
	Metadata metadata;

	metadata.add("Value1", 1.0);
	metadata.add("Value2", true);
	metadata.add("Value3", nullptr);
	metadata.add("Value4", "Hiho");

	EXPECT_EQ(4, metadata.size());

	EXPECT_TRUE(metadata.has("Value1"));
	EXPECT_TRUE(metadata.has("Value2"));
	EXPECT_TRUE(metadata.has("Value3"));
	EXPECT_TRUE(metadata.has("Value4"));

	EXPECT_FALSE(metadata.remove("adsdf"));

	EXPECT_EQ(4, metadata.size());

	EXPECT_TRUE(metadata.remove("Value1"));
	EXPECT_EQ(3, metadata.size());
	EXPECT_FALSE(metadata.has("Value1"));
	EXPECT_TRUE(metadata.has("Value2"));
	EXPECT_TRUE(metadata.has("Value2"));
	EXPECT_TRUE(metadata.has("Value3"));
	EXPECT_TRUE(metadata.has("Value4"));
}

TEST(MetadataTest, iterationEmpty)
{
	const Metadata metadata;

	EXPECT_EQ(0, metadata.size());
	EXPECT_TRUE(metadata.empty());

	EXPECT_TRUE(metadata.begin() == metadata.end());
}

TEST(MetadataTest, iterationValues)
{
	Metadata metadata;

	metadata.add("Value1", "bla");
	metadata.add("Value2", 1);
	metadata.add("Value3", false);

	EXPECT_EQ(3, metadata.size());
	EXPECT_FALSE(metadata.empty());

	EXPECT_EQ(3, std::distance(metadata.begin(), metadata.end()));
}

TEST(MetadataTest, empty)
{
	Metadata metadata;

	EXPECT_TRUE(metadata.empty());
	metadata.add("asdsd", "");
	EXPECT_FALSE(metadata.empty());
	metadata.remove("asdsd");
	EXPECT_TRUE(metadata.empty());
}

TEST(MetadataTest, indexOperator)
{
	Metadata metadata;
	EXPECT_TRUE(metadata.empty());
	metadata["Bla"] = "Test";
	EXPECT_EQ(1, metadata.size());
	metadata["Hiho"] = 1;
	EXPECT_EQ(2, metadata.size());
	metadata["Hiho"] = true;
	EXPECT_EQ(2, metadata.size());
}
