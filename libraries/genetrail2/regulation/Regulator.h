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

#ifndef GT2_REGULATOR_H
#define GT2_REGULATOR_H

#include <genetrail2/core/Exception.h>
#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/macros.h>

#include <vector>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>

namespace GeneTrail
{
template <typename ValueType> class GT2_EXPORT Regulator
{
  public:
	using value_type = ValueType;

	Regulator(const std::string& name, size_t numberOfTargets,
	          const value_type& original_score, size_t numberOfPermutations)
	    : name_(name),
	      numberOfTargets_(numberOfTargets),
	      original_score_(original_score)
	{
		permuted_scores_.reserve(numberOfPermutations);
	}

	Regulator(Regulator&&) = default;
	Regulator& operator=(Regulator&&) = default;
	Regulator(const Regulator&) = default;

	void addScore(const value_type& score)
	{
		permuted_scores_.emplace_back(score);
	}

	size_t numberOfTargets() const { return numberOfTargets_; }

	std::string name() const { return name_; }

	value_type original_score() const { return original_score_; }

	boost::math::normal distribution() const
	{
		value_type mean = statistic::mean<value_type>(permuted_scores_.begin(),
		                                              permuted_scores_.end());
		value_type sd = statistic::sd<value_type>(permuted_scores_.begin(),
		                                          permuted_scores_.end());
		boost::math::normal dist(mean, sd);
		return dist;
	}

	value_type estimateTwoSidedPValue() const
	{
		return 2.0 * boost::math::cdf(boost::math::complement(
		                 distribution(), fabs(original_score_)));
	}

	value_type estimateLowerTailedPValue() const
	{
		return boost::math::cdf(distribution(), original_score_);
	}

	value_type estimateUpperTailedPValue() const
	{
		return boost::math::cdf(
		    boost::math::complement(distribution(), original_score_));
	}

  private:
	std::string name_;
	size_t numberOfTargets_;
	value_type original_score_;
	std::vector<value_type> permuted_scores_;
};
}

#endif // GT2_REGULATOR_H
