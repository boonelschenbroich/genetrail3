/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2016 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_REGULATION_REGULATORY_IMPACT_SCORES_H
#define GT2_REGULATION_REGULATORY_IMPACT_SCORES_H

#include <boost/numeric/conversion/cast.hpp>

#include <genetrail2/core/Statistic.h>

namespace GeneTrail
{

struct GT2_EXPORT RIF1
{
	template <typename Iterator>
	double compute(Iterator reference_regulator_it,
	               Iterator reference_target_it, Iterator sample_regulator_it,
	               Iterator sample_target_it)
	{
		size_t ref_size = std::distance(reference_target_it->begin(),
		                                reference_target_it->end());
		size_t sample_size =
		    std::distance(sample_target_it->begin(), sample_target_it->end());
		double ref_target_sum = std::accumulate(
		    reference_target_it->begin(), reference_target_it->end(), 0.0);
		double sample_target_sum = std::accumulate(
		    sample_target_it->begin(), sample_target_it->end(), 0.0);
		double a = (ref_target_sum + sample_target_sum) /
		           boost::numeric_cast<double>(ref_size + sample_size);
		double d =
		    (ref_target_sum / boost::numeric_cast<double>(ref_size)) -
		    (sample_target_sum / boost::numeric_cast<double>(sample_size));
		double r1 = statistic::pearson_correlation<double>(
		    reference_regulator_it->begin(), reference_regulator_it->end(),
		    sample_regulator_it->begin(), sample_regulator_it->end());
		double r2 = statistic::pearson_correlation<double>(
		    reference_target_it->begin(), reference_target_it->end(),
		    sample_target_it->begin(), sample_target_it->end());
		return a * d * (r1 + r2);
	}
};

struct GT2_EXPORT RIF2
{
	template <typename Iterator>
	double compute(Iterator reference_regulator_it, Iterator reference_target_it,
	            Iterator sample_regulator_it, Iterator sample_target_it)
	{
		size_t ref_size = std::distance(reference_target_it->begin(),
		                                reference_target_it->end());
		size_t sample_size =
		    std::distance(sample_target_it->begin(), sample_target_it->end());
		double ref_target_sum = std::accumulate(
		    reference_target_it->begin(), reference_target_it->end(), 0.0);
		double sample_target_sum = std::accumulate(
		    sample_target_it->begin(), sample_target_it->end(), 0.0);
		double ex1 = (ref_target_sum / boost::numeric_cast<double>(ref_size));
		double ex2 =
		    (sample_target_sum / boost::numeric_cast<double>(sample_size));
		double r1 = statistic::pearson_correlation<double>(
		    reference_regulator_it->begin(), reference_regulator_it->end(),
		    sample_regulator_it->begin(), sample_regulator_it->end());
		double r2 = statistic::pearson_correlation<double>(
		    reference_target_it->begin(), reference_target_it->end(),
		    sample_target_it->begin(), sample_target_it->end());
		return pow((ex1 * r1), 2) - pow((ex2 * r2), 2);
	}
};
}

#endif // GT2_REGULATOR_EFFECT_ANALYSIS_H
