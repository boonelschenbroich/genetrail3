/**
* GeneTrail2 - An efficent library for interpreting genetic data
* Copyright (C) 2017 Tim Kehl tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_RTI_NETWORK_H
#define GT2_RTI_NETWORK_H

#include <genetrail2/core/macros.h>

#include "RegulatorEffectResult.h"
#include "RegulationFile.h"

#include <vector>
#include <set>
#include <iostream>
#include <fstream>

namespace GeneTrail
{
template <typename NameDatabase, typename ValueType> class GT2_EXPORT RTINetwork
{
  public:
	using value_type = ValueType;
	using Regulation = std::tuple<size_t, size_t, value_type>;

	RTINetwork(NameDatabase& name_database, std::vector<RegulatorEffectResult>& results, RegulationFile<value_type>& regulationFile)
	    : name_database_(name_database),
          results_(results),
          regulationFile_(regulationFile)
	{}

    void build_network() {
        for(RegulatorEffectResult rer : results_) {
            if (rer.corrected_p_value < 0.05) {
                for(auto& r : regulationFile_.regulator2regulations(name_database_(rer.name))) {
                    network_.emplace_back(r);
                    genes_.emplace(std::get<0>(r));
                    regulators_.emplace(std::get<0>(r));
                    genes_.emplace(std::get<1>(r));
                }
            }
        }
    }

    void save(const std::string& directory) {
        save_sif(directory);
        save_ids(directory, "names");
        save_ids(directory, "Gene-Symbol");
        save_types(directory);
        save_scores(directory, "");
        save_scores(directory, "_orig");
    }

  private:

    void save_sif(const std::string& directory) {
        std::ofstream out;
        out.open (directory + "/network.k1.sif");
        for(Regulation& r : network_) {
            std::string reg = (results_[std::get<1>(r)].mean_correlation > 0) ? "ACTIVATION" : "INHIBITION";
            out << name_database_(std::get<0>(r)) << "\t"
                << reg << "\t"
                << name_database_(std::get<1>(r)) << std::endl;
        }
        out.close();
    }

    void save_ids(const std::string& directory, const std::string& id) {
        std::ofstream out;
        out.open (directory + "/network." + id + ".na");
        out << "names (class=String)" << std::endl;
        for(auto r : genes_) {
            out << name_database_(r) << "\t=\t" << name_database_(r) << std::endl;
        }
        out.close();
    }

    void save_types(const std::string& directory) {
        std::ofstream out;
        out.open (directory + "/network.types.na");
        out << "types (class=String)" << std::endl;
        for(auto r : genes_) {
            std::string type = (regulators_.find(r) != regulators_.end()) ? "Regulator" : "Target";
            out << name_database_(r) << "\t=\t" << type << std::endl;
        }
        out.close();
    }

    void save_scores(const std::string& directory, const std::string& orig) {
        std::ofstream out;
        out.open (directory + "/network.score" + orig + ".na");
        out << "Scores (class=Double)" << std::endl;
        for(auto r : genes_) {
            double score = (regulators_.find(r) != regulators_.end()) ? results_[r].score : 0.0;
            score = (results_[r].mean_correlation < 0) ? score : score;
            out << name_database_(r) << "\t=\t" << score << std::endl;
        }
        out.close();
    }

	NameDatabase& name_database_;
    std::vector<RegulatorEffectResult>& results_;
    RegulationFile<value_type>& regulationFile_;
    std::vector<Regulation> network_;
    std::set<size_t> genes_;
    std::set<size_t> regulators_;
};
}

#endif // GT2_RTI_NETWORK_H
