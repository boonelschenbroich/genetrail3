#include <genetrail2/core/EntityDatabase.h>
#include <genetrail2/core/File.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/JsonCategoryFile.h>

#include <boost/program_options.hpp>

#include <iostream>

using namespace GeneTrail;
namespace bpo = boost::program_options;

void setParameters(const bpo::variables_map& vm, CategoryDatabase& db)
{
	if(!vm["editor-name"].empty()) {
		db.editor().name = vm["editor-name"].as<std::string>();
	}

	if(!vm["editor-email"].empty()) {
		db.editor().email = vm["editor-email"].as<std::string>();
	}

	if(!vm["creation-date"].empty()) {
		db.setCreationDate(vm["creation-date"].as<std::string>());
	}

	if(!vm["source-url"].empty()) {
		db.setSourceUrl(vm["source-url"].as<std::string>());
	}

	if(!vm["name"].empty()) {
		db.setName(vm["name"].as<std::string>());
	}

	if(!vm["identifier"].empty()) {
		db.setIdentifier(vm["identifier"].as<std::string>());
	}
}

template <typename CatDBFileA, typename CatDBFileB>
void readWriteCategoryDatabase(const std::string& input,
                               const std::string& output,
                               const bpo::variables_map& vm)
{
	auto db = std::make_shared<EntityDatabase>();
	CatDBFileA in(db, input);
	CatDBFileB out(db, output, FileOpenMode::WRITE);
	auto categories = in.read();
	setParameters(vm, categories);
	out.write(categories);
}

int main(int argc, char* argv[])
{
	std::string input, output;

	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()
		("input,i",         bpo::value(&input)->required(), "The input file that should be converted.")
		("output,o",        bpo::value(&output)->required(), "The output file that should be written.")
		("editor-name,e",   bpo::value<std::string>(), "Set the editor name.")
		("editor-email,m",  bpo::value<std::string>(), "Set the editor email.")
		("creation-date,c", bpo::value<std::string>(), "Set the creation date.")
		("source-url,u",    bpo::value<std::string>(), "Set the source URL.")
		("name,n",          bpo::value<std::string>(), "Set the database name.")
		("identifier,d",    bpo::value<std::string>(), "Set the identifier.")
	;

	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),
		           vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cout << e.what() << std::endl;
		return -1;
	}

	if(input.substr(input.length() - 5) == ".json") {
		readWriteCategoryDatabase<JsonCategoryFile, GMTFile>(input, output, vm);
	} else {
		readWriteCategoryDatabase<GMTFile, JsonCategoryFile>(input, output, vm);
	}

	return 0;
}
