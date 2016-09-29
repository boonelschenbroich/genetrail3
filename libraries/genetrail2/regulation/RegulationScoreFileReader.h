/**
* GeneTrail2 - An efficent library for interpreting genetic data
* Copyright (C) 2016 Tim Kehl tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_REGULATION_SCORE_FILE_READER_H
#define GT2_REGULATION_SCORE_FILE_READER_H

#include "Regulator.h"

#include <genetrail2/core/macros.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>

namespace GeneTrail
{
template <typename ValueType> class GT2_EXPORT RegulationScoreFileReader
{
  public:
	using value_type = ValueType;

	RegulationScoreFileReader(std::string file, size_t permutations, bool useAbsoluteValues)
	    : file_(file), permutations_(permutations), useAbsoluteValues_(useAbsoluteValues)
	{
	}

	RegulationScoreFileReader(RegulationScoreFileReader&&) = default;
	RegulationScoreFileReader& operator=(RegulationScoreFileReader&&) = default;
	RegulationScoreFileReader(const RegulationScoreFileReader&) = default;

	template <typename Function>
	std::vector<Regulator<value_type>> parse(Function aggregator)
	{
		read_();
		build_regulators_(aggregator);
		return regulators_;
	}

	std::vector<Regulator<value_type>> getRegulators() { return regulators_; }

	std::vector<value_type> getScores() { return scores_; }

  private:
	void read_()
	{
		std::ifstream input(file_);
		if(!input) {
			throw GeneTrail::IOError("File (" + file_ +
			                         ") is not open for reading");
		}

		for(std::string line; getline(input, line);) {
			std::vector<std::string> sline(2);
			boost::split(sline, line, boost::is_any_of(" \t"));
			if(sline.size() == 3) {
				boost::trim(sline[0]);
				boost::trim(sline[2]);
				value_type score = boost::lexical_cast<value_type>(sline[2]);
				if(regulator_2_scores_.find(sline[0]) ==
				   regulator_2_scores_.end()) {
					std::vector<value_type> tmp;
					regulator_2_scores_[sline[0]] = tmp;
				}
				if(useAbsoluteValues_){
					score = std::abs(score);
				}
				regulator_2_scores_[sline[0]].push_back(score);
				scores_.push_back(score);
			} else {
				throw GeneTrail::IOError("Wrong file format.");
			}
		}
	}

	template <typename Function> void build_regulators_(Function aggregator)
	{
		for(auto it = regulator_2_scores_.begin();
		    it != regulator_2_scores_.end(); ++it) {
			value_type score = aggregator(it->second.begin(), it->second.end());
			Regulator<value_type> r(it->first, it->second.size(), score,
			                        permutations_);
			regulators_.emplace_back(r);
		}
	}

	std::string file_;
	size_t permutations_;
	bool useAbsoluteValues_;
	std::map<std::string, std::vector<value_type>> regulator_2_scores_;
	std::vector<Regulator<value_type>> regulators_;
	std::vector<value_type> scores_;
};
}

#endif // GT2_REGULATION_SCORE_FILE_READER_H
