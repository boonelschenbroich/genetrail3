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

#ifndef GT2_REGULATION_REGULATOR_GENE_ASSOCIATION_ENRICHMENT_ANALYSIS_H
#define GT2_REGULATION_REGULATOR_GENE_ASSOCIATION_ENRICHMENT_ANALYSIS_H

#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/Exception.h>
#include <genetrail2/core/macros.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/Statistic.h>

#include <algorithm>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <fstream>

#include <boost/numeric/conversion/cast.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>

#include "RegulatorEffectResult.h"
#include "RegulationFile.h"
//#include "RegulationBootstrapper.h"

namespace GeneTrail
{
template <typename Bootstrapper, typename NameDatabase, typename ValueType> class GT2_EXPORT RegulatorGeneAssociationEnrichmentAnalysis
{
  public:
	using value_type = ValueType;
	using Regulation = std::tuple<size_t, size_t, value_type>;

	RegulatorGeneAssociationEnrichmentAnalysis(std::vector<size_t>& sorted_targets,
	                            RegulationFile<value_type>& regulationFile,
								Bootstrapper bootstraper,
                                NameDatabase name_database,
	                            bool use_absolute_value, size_t runs)
	    : sorted_targets_(sorted_targets),
	      regulationFile_(regulationFile),
		  bootstrapper_(bootstraper),
	      name_database_(name_database),
	      use_absolute_value_(use_absolute_value),
	      runs_(runs)
	{
		max_number_of_regulators = 0;
		biggest_regulator_idx = 0;
	}

	
	/**
	 * Performs the entire RegulatorGeneAssociationEnrichmentAnalysis.
	 *
	 * @param algorithm RegulatorEnrichmentAlgorithm that should be performed.
	 * @return A vector of RegulatorEffectResults.
	 */
	template <typename Algorithm, typename RegulatorImpactScore>
	std::vector<RegulatorEffectResult>
	run(Algorithm algorithm, RegulatorImpactScore impactScore)
	{
		// Perform algorithm without bootstrapping
		perform_bootstrapping_run_(impactScore);
		create_regulator_lists_();
		compute_scores_(algorithm);
		extract_correlations_();

		// Perform runs_ bootstrapping runs;
		for(size_t i = 0; i < runs_; ++i) {
			std::cout << "INFO: Performing bootstrapping run " << i + 1
			          << std::endl;
			bootstrapper_.create_bootstrap_sample();
			perform_bootstrapping_run_(impactScore);
			create_regulator_lists_();
			compute_scores_(algorithm);
		}

		// Compute p-value
		compute_pvalues_(algorithm);

		return results_;
	}

  private:
	/**
	 * Creates a sorted list of regulators REA algorithm.
	 */
	template <typename RegulatorImpactScore>
	void perform_bootstrapping_run_(RegulatorImpactScore impactScore)
	{
		std::vector<size_t> tf_list_tmp;
		size_t size_of_tf_list = 0;
		for(size_t targetname : sorted_targets_) {
			
			//std::cout << targetname << std::endl;
			// Check if gene is targetted by any Regulator
			if(!regulationFile_.checkTarget(targetname)) {
				continue;
			}

			auto& targets = regulationFile_.target2regulations(targetname);

			// Get the number of regulators
			size_of_tf_list += targets.size();

			max_number_of_regulators =
			    std::max(max_number_of_regulators, targets.size());

			// Perform single bootstrapping run
			bootstrapper_.perform_bootstrapping_run(targets,
			                                        use_absolute_value_, impactScore);
		}

		tf_list_tmp.reserve(size_of_tf_list);

		for(size_t i = 0; i < max_number_of_regulators; ++i) {
			for(size_t targetname : sorted_targets_) {
				// Check if gene is targetted by any Regulator
				if(!regulationFile_.checkTarget(targetname)) {
					continue;
				}

				std::vector<Regulation>& regulators = regulationFile_.target2regulations(targetname);

				if(regulators.size() <= i) {
					continue;
				}

				size_t regulator_i = std::get<0>(regulators[i]);

				biggest_regulator_idx =
				    std::max(regulator_i, biggest_regulator_idx);

				tf_list_tmp.emplace_back(regulator_i);
			}
		}

		tf_list = tf_list_tmp;
	}

	/**
	 * Save for each regulator the indices of the occurrences in the list.
	 */
	void create_regulator_lists_()
	{
		std::vector<std::vector<size_t>> regulator_inidces_tmp(biggest_regulator_idx + 1);
		for(size_t i = 0; i < tf_list.size(); ++i) {
			regulator_inidces_tmp[tf_list[i]].emplace_back(i);
		}
		
		regulator_inidces_ = regulator_inidces_tmp;
	}

	/**
	 * Perform the REA algorithm for each regulator.
	 *
	 * @param algorithm The algorithm to be performed (ks-test, wrs-test)
	 */
	template <typename Algorithm> void compute_scores_(Algorithm algorithm)
	{
		bool empty = results_.size() == 0;
		if(empty) {
			results_.resize(biggest_regulator_idx + 1);
		}
		for(size_t i = 0; i < regulator_inidces_.size(); ++i) {
			if(regulator_inidces_[i].size() == 0){
				continue;
			}
			auto score = algorithm.compute_score(tf_list.size(),
			                                     regulator_inidces_[i].begin(),
			                                     regulator_inidces_[i].end());

			if(empty) {
				results_[i].name = name_database_(i);
				results_[i].hits = regulator_inidces_[i].size();
				results_[i].scores.reserve(runs_);
			}

			results_[i].addScore(score);
		}
	}

	/**
	 * Extracts the mean correlations so we judge what effect the regulator has.
	 */
	void extract_correlations_()
	{
		regulator2correlations_.resize(biggest_regulator_idx + 1);
		for(size_t targetname : sorted_targets_) {
			// Check if gene is targetted by any Regulator
			if(!regulationFile_.checkTarget(targetname)) {
				continue;
			}

			for(const auto& reg : regulationFile_.target2regulations(targetname)) {
				if(std::get<2>(reg) <= 1.0 && std::get<2>(reg) >= -1.0 &&
				   std::get<2>(reg) != 0.0) {
					regulator2correlations_[std::get<0>(reg)].emplace_back(
					    std::get<2>(reg));
				}
			}
		}

		for(size_t i = 0; i < regulator2correlations_.size(); ++i) {
			results_[i].mean_correlation =
			    statistic::mean<value_type>(regulator2correlations_[i].begin(),
			                                regulator2correlations_[i].end());
		}
	}

	/**
	 * Computes a p-value for the middle element of all the bootstrap samples.
	 */
	template <typename Algorithm> void compute_pvalues_(Algorithm algorithm)
	{
		for(auto& result : results_) {
			if(result.name == "") {
				continue;
			}
			
			std::cout << "INFO: Computing p-values for category '"
			          << result.name << "'" << std::endl;
			result.score = statistic::middle<value_type>(result.scores.begin(),
			                                             result.scores.end());
			algorithm.compute_p_value(result, tf_list.size());
			
			result.normalize_scores(algorithm, tf_list.size());
		}
	}

	std::vector<size_t> sorted_targets_;
	RegulationFile<value_type>& regulationFile_;
	Bootstrapper bootstrapper_;
	NameDatabase name_database_;

	bool use_absolute_value_;
	size_t runs_;

	size_t max_number_of_regulators;
	size_t biggest_regulator_idx;

	std::vector<std::vector<value_type>> regulator2correlations_;
	std::vector<std::vector<size_t>> regulator_inidces_;
	std::vector<RegulatorEffectResult> results_;
	std::vector<size_t> tf_list;
};
}

#endif // GT2_REGULATION_REGULATOR_GENE_ASSOCIATION_ENRICHMENT_ANALYSIS_H
