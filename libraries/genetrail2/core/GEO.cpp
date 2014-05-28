#include "GEO.h"

using namespace GeneTrail;

GEO::GEO()
{
}

GEO::~GEO()
{
}

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
	if( method == "mean" ){
		return statistic::mean<double, _viter>(values.begin(), values.end());
	}
	else if( method == "mean" ){
		return statistic::median<double, _viter>(values.begin(), values.end());
	}
	else if( method == "max" ){
		return statistic::max<double, _viter>(values.begin(), values.end());
	}
	else if( method == "min" ){
		return statistic::min<double, _viter>(values.begin(), values.end());
	}
	return 0.0;
}

void GEO::writeGEOMap(const std::string& filename, const GEOMap& map)
{
	std::ofstream out(filename.c_str());
	bool first = true;
	for(auto id : map.sampleNames)
	{
		if(first){
			out << id;
			first = false;
		}
		out << "\t" << id;
	}
	out << std::endl;

	for(auto it = map.gene2exprs.begin(); it != map.gene2exprs.end(); it++) {
		if(it->first != "") {
			std::vector<std::string> strs;
			boost::algorithm::split_regex(strs, it->first, boost::regex("///"));
			out << strs[0];
			for(unsigned int i = 0; i < it->second.size(); ++i) {
				out << "\t" << it->second[i];
			}
			out << std::endl;
		}
	}
	out.close();
}
