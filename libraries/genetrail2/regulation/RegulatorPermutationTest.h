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

#ifndef GT2_REGULATOR_PERMUTATION_TEST_H
#define GT2_REGULATOR_PERMUTATION_TEST_H

#include "Regulator.h"

#include <genetrail2/core/macros.h>

#include <vector>
#include <algorithm>
#include <random>
#include <iostream>

namespace GeneTrail
{
template <typename ValueType> class GT2_EXPORT RegulatorPermutationTest
{
  public:
	using value_type = ValueType;

	RegulatorPermutationTest(const std::vector<Regulator<value_type>>& regulators,
	                         const std::vector<value_type>& scores,
	                         size_t permutations, uint64_t randomSeed)
	    : regulators_(regulators),
	      scores_(scores),
	      permutations_(permutations),
	      twister_(randomSeed)
	{
		sort(
		    regulators_.begin(), regulators_.end(),
		    [](const Regulator<value_type>& l, const Regulator<value_type>& r) {
			    return l.numberOfTargets() < r.numberOfTargets();
			});
	}

	RegulatorPermutationTest(RegulatorPermutationTest&&) = default;
	RegulatorPermutationTest& operator=(RegulatorPermutationTest&&) = default;
	RegulatorPermutationTest(const RegulatorPermutationTest&) = default;

	template <typename Function>
	std::vector<Regulator<value_type>> execute(Function aggregator)
	{
		performAllPermutations_(aggregator);
		return regulators_;
	}

  private:
	void shuffle_(size_t n)
	{
		for(size_t i = 0; i < n; ++i) {
			std::swap(scores_[i],
			          scores_[i + twister_() % (scores_.size() - i)]);
		}
	}

	template <typename Function>
	void performAllPermutations_(Function aggregator)
	{
		for(size_t i = 0; i < permutations_; ++i) {
			std::cout << "INFO: Permutation " << i << "\n";
			performSinglePermutation_(aggregator);
		}
	}

	template <typename Function>
	void performSinglePermutation_(Function aggregator)
	{
		size_t lastNumberOfTargets = 0;
		value_type lastScore = 0.0;
		for(size_t i=0; i<regulators_.size(); ++i) {
			if(lastNumberOfTargets < regulators_[i].numberOfTargets()) {
				shuffle_(regulators_[i].numberOfTargets());
				lastScore = aggregator(scores_.begin(),
				                       scores_.begin() + regulators_[i].numberOfTargets());
			}
			regulators_[i].addScore(lastScore);
		}
	}

	std::vector<Regulator<value_type>> regulators_;
	std::vector<value_type> scores_;
	size_t permutations_;
	std::mt19937_64 twister_;
};
}

#endif // GT2_REGULATOR_PERMUTATION_TEST_H
