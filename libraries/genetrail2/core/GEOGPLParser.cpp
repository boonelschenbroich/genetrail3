/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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
#include "GEOGPLParser.h"

using namespace GeneTrail;

GPL_Parser::GPL_Parser() : GEO()
{
}

// default destructor
GPL_Parser::~GPL_Parser()
{
}

void GPL_Parser::downloadGPLFile(const std::string& filename,
                                 const std::string& geo_dir)
{
	int r = system(
	    ("wget ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/annotation/platforms/" +
	     filename + " -P " + geo_dir).c_str());
	if (r != 0)
	{
		std::cerr << "Could not download GPL file." << std::endl;
	}
}

std::map<std::string, std::string>
GPL_Parser::readGPLFile(const std::string& filename,
                        const std::string& mappingColumn)
{
	bool within_platform_table = false;
	uint index_col = 999;

	std::map<std::string, std::string> cloneid2otherid_;

	std::ifstream file(filename.c_str(),
	                   std::ios_base::in | std::ios_base::binary);

	if(!file) {
		std::cerr << "GSE_Parser::readPlatformFile cannot open file: "
		          << filename << std::endl;
		return cloneid2otherid_;
	}

	boost::iostreams::filtering_istream input;
	if(filename.find(".gz") != std::string::npos) {
		input.push(boost::iostreams::gzip_decompressor());
		input.push(file);
	} else {
		input.push(file);
	}

	for(std::string line; std::getline(input, line);) {
		if(line != "") {
			if(line.find("platform_table_begin") != std::string::npos) {
				within_platform_table = true;
				// retrieve the next line containing the headers of the columns
				getline(input, line);
				std::vector<std::string> entries;
				boost::split(entries, line, boost::is_any_of("\t"));
				// find the "mapping" column
				for(uint i = 0; i < entries.size(); ++i) {
					if(entries[i] == mappingColumn) {

						std::cout << i << "\t" << entries[i] << std::endl;
						index_col = i;
					}
				}
				// check if we have found the preferred mapping column
				if(index_col == 999) {
					std::cerr << "GSE_Parser::readPlatformFile cannot find the "
					             "mapping header: " << mappingColumn
					          << std::endl;
					return cloneid2otherid_;
				}
				continue;
			}
			if(within_platform_table) {
				std::vector<std::string> entries;
				boost::split(entries, line, boost::is_any_of("\t"));
				if(entries.size() >= index_col) {
					if((boost::trim_copy(entries[0]) != "") &&
					   (boost::trim_copy(entries[index_col]) != "")) {
						if(cloneid2otherid_.find(entries[0]) ==
						   cloneid2otherid_.end()) {
							cloneid2otherid_[entries[0]] = entries[index_col];
						} else {
							std::cerr << "warning: "
							             "GSE_Parser::readPlatformFile -> same "
							             "cloneid mapped to different other "
							             "ids!" << std::endl;
						}
					}
				}
			}

			if(line.find("platform_table_end") != std::string::npos) {
				break;
			}
		}
	}
	input.auto_close();
	return cloneid2otherid_;
}

std::map<std::string, std::string>
GPL_Parser::annotate(const std::string& geo_dir,
					 const std::string& platform,
                     const std::string& mappingColumn)
{
	std::string filename = platform + ".annot.gz";
	std::string path = geo_dir + "/" + filename;
	std::ifstream file(path.c_str(), std::ios_base::in | std::ios_base::binary);

	if(!file) {
		downloadGPLFile(filename, geo_dir);
		std::cout << "Downloading: " << filename << " to " << geo_dir
		          << std::endl;
	}

	boost::iostreams::filtering_istream input;
	if(filename.find(".gz") != std::string::npos) {
		input.push(boost::iostreams::gzip_decompressor());
		input.push(file);
	} else {
		input.push(file);
	}

	if(input.empty()) {
		downloadGPLFile(filename, geo_dir);
		std::cout << "Downloading: " << filename << " to " << geo_dir
		          << std::endl;
	}

	input.auto_close();

	return readGPLFile(path, mappingColumn);
}
