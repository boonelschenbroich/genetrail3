/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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

#include "SetLevelStatistics.h"

namespace GeneTrail
{
	StatisticsEnrichment::StatisticsEnrichment(const Statistics& test,
	                                           const Scores& scores)
	    : test_(test), scores_(scores.db()), cached_(0.0)
	{
		setInputScores(scores);
	}

	StatisticsEnrichment::StatisticsEnrichment(const Statistics& test,
	                                           const Scores& scores,
	                                           double cached)
	    : test_(test), scores_(scores), cached_(cached)
	{
		scores_.sortByIndex();
	}

	void StatisticsEnrichment::setInputScores(const Scores& scores)
	{
		scores_ = scores;
		scores_.sortByIndex();
		cached_ = cacheStatistic_();
	}

	std::tuple<double, double>
	StatisticsEnrichment::computeScore(const Category& c) const
	{
		const auto intersection = scores_.subset(c);
		auto score =
		    test_(intersection.scores().begin(), intersection.scores().end());

		return std::make_tuple(score, getExpectedValue_(c));
	}

	double StatisticsEnrichment::cacheStatistic_() const
	{
		return test_(scores_.scores().begin(), scores_.scores().end());
	}

	double StatisticsEnrichment::getExpectedValue_(const Category&) const
	{
		return cached_;
	}

	double SumEnrichment::cacheStatistic_() const
	{
		return statistic::mean<double>(scores_.scores().begin(),
		                               scores_.scores().end());
	}
}
