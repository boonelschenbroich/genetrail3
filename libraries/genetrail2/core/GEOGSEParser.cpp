#include <boost/lexical_cast.hpp>
#include <limits>

#include "GEOGSEParser.h"

using namespace GeneTrail;

// default ctor
GEOGSEParser::GEOGSEParser() : GEO()
{
}

// default destructor
GEOGSEParser::~GEOGSEParser()
{
}

GEOMap
GEOGSEParser::readGSEFile(const std::string& filename)
{
	GEOMap result;
	std::map<std::string, std::vector<double>> clone2expression_values;
	std::string gsm = "";
	bool within_sample_table = false;

	std::ifstream file(filename.c_str(),
	                   std::ios_base::in | std::ios_base::binary);
	if(!file) {
		std::cerr << "ERROR: Cannot open file: " << filename << std::endl;
		return result;
	}

	boost::iostreams::filtering_istream input;
	if(filename.find(".gz") != std::string::npos) {
		input.push(boost::iostreams::gzip_decompressor());
		input.push(file);
	} else {
		input.push(file);
	}

	int value_idx = -1;
	for(std::string line; std::getline(input, line);) {
		if(line != "") {
			if(boost::starts_with(line, "!Series_platform_id")) {
				result.platform = line.substr(line.find("GPL"));
				std::cout << "INFO: Platform: " << result.platform << std::endl;
			}
			// start extraction for found sample
			if(boost::starts_with(line, "^SAMPLE")) {
				if(line.find("GSM") != std::string::npos) {
					// extract GSM number
					gsm = line.substr(line.find("GSM"));
					result.sampleNames.push_back(gsm);
					std::cout << "INFO: Parsing - " << gsm << std::endl;
					continue;
				}
			}

			if(boost::starts_with(line, "!sample_table_begin")) {
				within_sample_table = true;
				continue;
			}

			if(within_sample_table && (gsm != "")) {

				if(line[0] == '#') {
					continue;
				}
				// skip header line
				if(boost::starts_with(line, "ID_REF")) {
					std::vector<std::string> entries;
					boost::split(entries, line, boost::is_any_of("\t"),
					             boost::token_compress_on);

					for(value_idx = 0; (size_t)value_idx < entries.size() &&
					                       (entries[value_idx] != "VALUE");
					    ++value_idx) {
					}

					continue;
				}

				if(value_idx == -1) {
					continue;
				}

				typedef boost::split_iterator<std::string::iterator>
				string_split_iterator;

				std::string probeset;
				std::string value;

				int i = 0;
				for(string_split_iterator
				        it = boost::make_split_iterator(
				            line, boost::first_finder("\t", boost::is_equal()));
				    (i <= value_idx) && (it != string_split_iterator());
				    ++it, ++i) {
					if(i == 0) {
						probeset = boost::copy_range<std::string>(*it);
					} else if(i == value_idx) {
						value = boost::copy_range<std::string>(*it);
					}
				}

				if(i > value_idx) {
					double expression_value;
					try
					{
						expression_value = boost::lexical_cast<double>(value);
					}
					catch(boost::bad_lexical_cast& e)
					{
						expression_value =
						    std::numeric_limits<double>::quiet_NaN();
					}
					// new clone id
					if(clone2expression_values.find(probeset) ==
					   clone2expression_values.end()) {
						if(boost::trim_copy(probeset) != "") {
							std::vector<double> tmp_vec;
							tmp_vec.push_back(expression_value);
							clone2expression_values[probeset] = tmp_vec;
						}
					} else // clone id is known... just append the expression
					       // value
					{
						if(boost::trim_copy(probeset) != "") {
							clone2expression_values[probeset]
							    .push_back(expression_value);
						}
					}
				}
			}

			if(boost::starts_with(line, "!sample_table_end")) {
				within_sample_table = false;
				gsm = "";
				value_idx = -1;
			}
		}
	}

	input.auto_close();
	result.gene2exprs = clone2expression_values;
	return result;
}

