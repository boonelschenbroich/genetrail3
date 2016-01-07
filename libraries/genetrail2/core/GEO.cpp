/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *               2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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
#include "GEO.h"

#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/iter_find.hpp>

using namespace GeneTrail;

GEO::GEO() {}

GEO::~GEO() {}

std::vector<double>
GEO::mergeDuplicatedVectors(std::vector<std::vector<double>> matrix,
                            std::string method)
{
	std::vector<double> merged;
	for(unsigned int i = 0; i < matrix[0].size(); ++i) {
		std::vector<double> tmp;
		for(unsigned int j = 0; j < matrix.size(); ++j) {
			tmp.push_back(matrix[j][i]);
		}
		double merged_value = apply(method, tmp);
		merged.push_back(merged_value);
	}
	return merged;
}

std::map<std::string, std::vector<double>> GEO::mapAndRemoveDuplicates(
    std::map<std::string, std::vector<double>> gene2expression_values_,
    std::map<std::string, std::string> cloneid2otherid_, std::string method)
{
	std::map<std::string, std::vector<std::vector<double>>> matrix;
	for(std::map<std::string, std::vector<double>>::iterator it =
	        gene2expression_values_.begin();
	    it != gene2expression_values_.end(); ++it) {
		if(matrix.find(cloneid2otherid_[it->first]) == matrix.end()) {
			std::vector<std::vector<double>> new_tmp;
			new_tmp.push_back(it->second);
			matrix[cloneid2otherid_[it->first]] = new_tmp;
		} else {
			matrix[cloneid2otherid_[it->first]].push_back(it->second);
		}
	}
	std::map<std::string, std::vector<double>> return_map;
	for(std::map<std::string, std::vector<std::vector<double>>>::iterator it =
	        matrix.begin();
	    it != matrix.end(); ++it) {
		if(it->second.size() == 1) {
			return_map[it->first] = it->second[0];
		} else {
			return_map[it->first] = mergeDuplicatedVectors(it->second, method);
		}
	}
	return return_map;
}

double GEO::apply(std::string method, std::vector<double> values)
{
	if(method == "mean") {
		return statistic::mean<double>(values.begin(), values.end());
	} else if(method == "median") {
		return statistic::median<double>(values.begin(), values.end());
	} else if(method == "max") {
		return statistic::max<double>(values.begin(), values.end());
	} else if(method == "min") {
		return statistic::min<double>(values.begin(), values.end());
	}
	return 0.0;
}

void GEO::writeGEOMap(const std::string& filename, const GEOMap& map)
{
	std::ofstream out(filename.c_str());
	bool first = true;
	for(auto id : map.sampleNames) {
		if(first) {
			out << id;
			first = false;
			continue;
		}
		out << "\t" << id;
	}
	out << std::endl;

	for(auto it = map.gene2exprs.begin(); it != map.gene2exprs.end(); it++) {
		if(it->first != "") {
			std::vector<std::string> strs;
			boost::algorithm::iter_split(strs, it->first, boost::first_finder("///"));
			out << strs[0];
			for(unsigned int i = 0; i < it->second.size(); ++i) {
				out << "\t" << it->second[i];
			}
			out << std::endl;
		}
	}
	out.close();
}
