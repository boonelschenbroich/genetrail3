/*
 * GeneTrail2 - An efficent library for interpreting genetic data
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
#include "GEOGDSParser.h"

using namespace GeneTrail;

GEOGDSParser::GEOGDSParser()
{
}

GEOGDSParser::~GEOGDSParser()
{
}

/**
 * This function parses a GDS file and saves all columns.
 * @param filename The name of the GDS file.
 */
GEOMap GEOGDSParser::readGDSFile(const std::string& filename)
{
	GEOMap result;
	std::map<std::string, std::vector<double>> clone2expression_values;
	unsigned int sample_count = 0;

	std::ifstream file(filename.c_str(),
	                   std::ios_base::in | std::ios_base::binary);
	if(!file) {
		std::cerr << "ERROR: Cannot open GDS file : " << filename << std::endl;
		return result;
	}

	std::cout << "INFO: Parsing - " << filename << std::endl;

	boost::iostreams::filtering_istream gds_stream;
	if(filename.find(".gz") != std::string::npos) {
		gds_stream.push(boost::iostreams::gzip_decompressor());
		gds_stream.push(file);
	} else {
		gds_stream.push(file);
	}

	// Flag to indicate if the program is in the header section
	bool header_already_done = false;

	// Start reading file line by line
	for(std::string line; std::getline(gds_stream, line);) {
		// The line containing the words ID_REF IDENTIFIER is the end of the header
		if(!header_already_done) {
			// Parse needed information from the header setion
			if(line.find("dataset_platform =") != std::string::npos) {
				std::string s = line;
				result.platform = s.substr(s.find("GPL"));
				std::cout << "INFO: Platform: " << result.platform << std::endl;
			} else if(line.find("dataset_sample_count") != std::string::npos) {
				std::vector<std::string> strs;
				std::string s = line;
				boost::split(strs, s, boost::is_any_of("="));
				sample_count = atoi(strs[1].c_str());
				std::cout << "INFO: Number of samples: " << sample_count << std::endl;
			} else if(line.find("DATASET") != std::string::npos){
				std::vector<std::string> strs;
				boost::split(strs, line, boost::is_any_of("="));
				boost::trim(strs[1]);
				result.dataset = strs[1];
			} else if(line.find("ID_REF\tIDENTIFIER") != std::string::npos) {
				header_already_done = true;
				// Save the contained GSMs
				std::vector<std::string> gsms;
				boost::split(gsms, line, boost::is_any_of(" \t"));
				// Delete the ID_REF IDENTIFIER tags
				gsms.erase(gsms.begin(),gsms.begin()+2);
				result.sampleNames = gsms;
			}
			continue;
		}

		// The line containing the text "_table_end" indicates the file end
		else if(line.find("_table_end") != std::string::npos)
			break;

		// A vector to save the entries of read line
		std::vector<std::string> tabs;

		boost::algorithm::split(tabs, line, [](char c) { return c == '\t'; });

		for(unsigned int i = 0; i < sample_count; i++) {
			double expression_value;
			if(tabs.at(i + 2) == "null")
				expression_value = std::numeric_limits<double>::quiet_NaN();
			else
				expression_value = boost::lexical_cast<double>(tabs.at(i + 2));

			// The probe id is the first entry in 'tabs'
			std::string id = tabs.at(0);

			// cloneid2otherid_.insert(std::pair<std::string,std::string>(tabs.at
			// ( 0 ),tabs.at ( 1 )));

			// find out whether this id has been read before (from files read earlier in this analysis)
			std::map<std::string, std::vector<double>>::iterator id_iterator = clone2expression_values.find(id);

			// if id is being read for first time, make a vector, save the just
			// read expression value into the vector
			// and save the combination of the id and the vector into the right
			// map
			// if not, save the just read value into the existing vector
			if(id_iterator == clone2expression_values.end()) {
				std::vector<double> values;
				values.push_back(expression_value);
				clone2expression_values.insert(make_pair(id, values));
			} else {
				id_iterator->second.push_back(expression_value);
			}
		}
	}

	gds_stream.auto_close();
	result.gene2exprs = clone2expression_values;
	return result;
}

