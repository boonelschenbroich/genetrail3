/*
 * GeneTrail2 - An efficent library for interpreting genetic data
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
#include "JsonCategoryFile.h"

#include "Exception.h"

#include <rapidjson/document.h>
#include <rapidjson/writer.h>

#include <cassert>
#include <cstring>

namespace GeneTrail
{
// This code is taken from the RapidJSON documentation
class IStreamWrapper
{
  public:
	typedef char Ch;
	IStreamWrapper(std::istream& is) : is_(is) {}
	Ch Peek() const
	{
		int c = is_.peek();
		return c == std::char_traits<char>::eof() ? '\0' : (Ch)c;
	}
	Ch Take()
	{
		int c = is_.get();
		return c == std::char_traits<char>::eof() ? '\0' : (Ch)c;
	}
	size_t Tell() const { return (size_t)is_.tellg(); }
	Ch* PutBegin()
	{
		assert(false);
		return 0;
	}
	void Put(Ch) { assert(false); }
	void Flush() { assert(false); }
	size_t PutEnd(Ch*)
	{
		assert(false);
		return 0;
	}

  private:
	IStreamWrapper(const IStreamWrapper&);
	IStreamWrapper& operator=(const IStreamWrapper&);
	std::istream& is_;
};

class OStreamWrapper
{
  public:
	typedef char Ch;
	OStreamWrapper(std::ostream& os) : os_(os) {}
	Ch Peek() const
	{
		assert(false);
		return '\0';
	}
	Ch Take()
	{
		assert(false);
		return '\0';
	}
	size_t Tell() const { return 0; }
	Ch* PutBegin()
	{
		assert(false);
		return 0;
	}
	void Put(Ch c) { os_.put(c); }
	void Flush() { os_.flush(); }
	size_t PutEnd(Ch*)
	{
		assert(false);
		return 0;
	}

  private:
	OStreamWrapper(const OStreamWrapper&);
	OStreamWrapper& operator=(const OStreamWrapper&);
	std::ostream& os_;
};


using Writer = rapidjson::Writer<OStreamWrapper>;

static const char* getString(const rapidjson::Value& member, const char* memberName)
{
	if(!member.IsString()) {
		throw IOError(std::string("Category field '") + memberName +
		              "' is not a string.");
	}

	return member.GetString();
}

static Metadata::Value toMetadata(const rapidjson::Value& value)
{
	if(value.IsNull()) {
		return Metadata::Value(nullptr);
	} else if(value.IsString()) {
		return Metadata::Value(value.GetString());
	} else if(value.IsDouble()) {
		return Metadata::Value(value.GetDouble());
	} else if(value.IsBool()) {
		return Metadata::Value(value.GetBool());
	} else if(value.IsInt()) {
		return Metadata::Value(value.GetInt());
	} else if(value.IsInt64()) {
		return Metadata::Value(value.GetInt64());
	} else if(value.IsUint()) {
		return Metadata::Value(value.GetUint());
	} else if(value.IsUint64()) {
		return Metadata::Value(value.GetUint64());
	} else if(value.IsArray()) {
		Metadata::Array values;
		values.reserve(value.Size());
		for(auto it = value.Begin(); it != value.End(); ++it) {
			values.emplace_back(toMetadata(*it));
		}

		return Metadata::Value(std::move(values));
	} else if(value.IsObject()) {
		Metadata::Object values;
		values.reserve(value.MemberCount());
		for(auto it = value.MemberBegin(); it != value.MemberEnd(); ++it) {
			values.emplace(it->name.GetString(), toMetadata(it->value));
		}

		return Metadata::Value(std::move(values));
	} else {
		throw NotImplemented(__FILE__, __LINE__, "Unhandled case for JSON Value.");
	}
}

static void readCategoryFields(const rapidjson::Value& category, Category& cat)
{
	if(!category.HasMember("name")) {
		throw IOError("Found category without name.");
	}

	if(!category.HasMember("url")) {
		throw IOError("Found category without URL.");
	}

	for(auto it = category.MemberBegin(); it != category.MemberEnd(); ++it) {
		const char* name = it->name.GetString();
		if(std::strcmp(name, "name") == 0) {
			cat.setName(getString(it->value, "name"));
		} else if(std::strcmp(name, "url") == 0) {
			cat.setReference(getString(it->value, "url"));
		} else if(std::strcmp(name, "members") != 0) {
			cat.metadata().add(name, toMetadata(it->value));
		}
	}
}

static void readMembers(const rapidjson::Value& members, Category& cat)
{
	for(auto it = members.Begin(); it != members.End(); ++it) {
		if(!it->IsString()) {
			throw IOError("Found a non-string element im members array.");
		}

		cat.insert(it->GetString());
	}
}

static void readCategory(const rapidjson::Value& category, Category& cat)
{
	if(!category.HasMember("members") || !category["members"].IsArray()) {
		throw IOError("Category does not have a member array.");
	}

	readCategoryFields(category, cat);
	readMembers(category["members"], cat);
}

static void readCategories(const rapidjson::Value& categories,
                           CategoryDatabase& database)
{
	database.reserve(categories.Size());

	for(rapidjson::SizeType i = 0; i < categories.Size(); ++i) {
		if(!categories[i].IsObject()) {
			throw IOError("Entry " + std::to_string(i) + " is not an object.");
		}

		Category& cat = database.addCategory();
		readCategory(categories[i], cat);
	}
}


static void readDatabaseFields(const rapidjson::Document& document,
                               CategoryDatabase& database)
{
	if(!document.HasMember("name")) {
		throw IOError("Database does not have a name.");
	}

	if(!document.HasMember("editor") || !document["editor"].IsObject()) {
		throw IOError("Database does not have an editor.");
	}

	if(!document.HasMember("creationDate")) {
		throw IOError("Database does not have a creation date.");
	}

	if(!document.HasMember("sourceUrl")) {
		throw IOError("Database does not have a source url.");
	}

	if(!document.HasMember("identifier")) {
		throw IOError("Database does not have a identifier.");
	}

	for(auto it = document.MemberBegin(); it != document.MemberEnd(); ++it) {
		const auto* name = it->name.GetString();

		if(std::strcmp(name, "name") == 0) {
			database.setName(getString(it->value, "name"));
		} else if(std::strcmp(name, "editor") == 0) {
			database.editor().name = getString(it->value["name"], "editor.name");
			database.editor().email = getString(it->value["email"], "editor.email");
		} else if(std::strcmp(name, "creationDate") == 0) {
			database.setCreationDate(getString(it->value, "creationDate"));
		} else if(std::strcmp(name, "sourceUrl") == 0) {
			database.setSourceUrl(getString(it->value, "sourceUrl"));
		} else if(std::strcmp(name, "identifier") == 0) {
			database.setIdentifier(getString(it->value, "identifier"));
		} else if(std::strcmp(name, "categories") != 0) {
			database.metadata().add(name, toMetadata(it->value));
		}
	}
}

static void buildDatabase(const rapidjson::Document& document,
                          CategoryDatabase& database)
{
	if(!document.HasMember("categories") || !document["categories"].IsArray()) {
		throw IOError("Provided category database does not contain a list of "
		              "categories.");
	}

	readDatabaseFields(document, database);
	readCategories(document["categories"], database);
}

static void writeString(Writer& writer, const std::string& key)
{
	writer.String(key.c_str(), key.size());
}

class MetadataWriter: public boost::static_visitor<void>
{
public:
	MetadataWriter(Writer& writer)
		: writer_(writer){
	}

	void operator()(const std::string& value) const {
		writer_.String(value.c_str(), value.length());
	}

	void operator()(nullptr_t) const {
		writer_.Null();
	}

	void operator()(bool value) const {
		writer_.Bool(value);
	}

	void operator()(int64_t value) const {
		writer_.Int64(value);
	}

	void operator()(double value) const {
		writer_.Double(value);
	}

	void operator()(const Metadata::Array& value) const {
		writer_.StartArray();
		for(const auto& elem : value) {
			boost::apply_visitor(MetadataWriter(writer_), *elem);
		}
		writer_.EndArray();
	}

	void operator()(const Metadata::Object& value) const {
		writer_.StartObject();
		for(const auto& elem : value) {
			writer_.Key(elem.first.c_str(), elem.first.length());
			boost::apply_visitor(MetadataWriter(writer_), *elem.second);
		}
		writer_.EndObject();
	}

private:
	Writer& writer_;
};

static void writeFreeFields(Writer& writer, const Metadata& metadata)
{
	for(const auto& m : metadata) {
		writer.Key(m.first.c_str(), m.first.length());
		boost::apply_visitor(MetadataWriter(writer), *m.second);
	}
}

static void writeCategory(Writer& writer, const Category& cat)
{
	writer.StartObject();
	writer.Key("name", 4);
	writeString(writer, cat.name());
	writer.Key("url", 3);
	writeString(writer, cat.reference());

	writer.Key("members", 7);
	writer.StartArray();
	for(const auto& e : cat.names()) {
		writeString(writer, e);
	}
	writer.EndArray();

	writeFreeFields(writer, cat.metadata());

	writer.EndObject();
}

static void writeDatabaseFields(Writer& writer, const CategoryDatabase& db)
{
	writer.Key("name", 4);
	writeString(writer, db.name());

	writer.Key("identifier", 10);
	writeString(writer, db.identifier());

	writer.Key("sourceUrl", 9);
	writeString(writer, db.sourceUrl());

	writer.Key("creationDate", 12);
	writeString(writer, db.creationDate());

	writer.Key("editor", 6);
	writer.StartObject();
	writer.Key("name", 4);
	writeString(writer, db.editor().name);
	writer.Key("email", 5);
	writeString(writer, db.editor().email);
	writer.EndObject();

	writer.Key("categories", 10);
	writer.StartArray();
	for(const auto& cat : db) {
		writeCategory(writer, cat);
	}
	writer.EndArray();

	writeFreeFields(writer, db.metadata());
}

JsonCategoryFile::JsonCategoryFile(const std::shared_ptr<EntityDatabase>& db,
                                   const std::string& path, FileOpenMode mode)
    : CategoryDatabaseFile(path, mode), entity_database_(db)
{
}

bool JsonCategoryFile::write(const CategoryDatabase& db)
{
	if(!isValid_() || !isWriting()) {
		throw IOError("The file object is not in a valid state to write.");
	}

	OStreamWrapper strm(*out_strm_);
	Writer writer(strm);

	writer.StartObject();
	writeDatabaseFields(writer, db);
	writer.EndObject();

	return true;
}

CategoryDatabase JsonCategoryFile::read()
{
	CategoryDatabase database(entity_database_);

	rapidjson::Document document;
	IStreamWrapper strm(*in_strm_);
	document.ParseStream(strm);

	if(document.GetParseError() != rapidjson::kParseErrorNone) {
		throw IOError("Invalid Json file!");
	}

	buildDatabase(document, database);

	return database;
}
}
