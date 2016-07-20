/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2015 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_REGULATION_REGULATOR_ENRICHMENT_RESULT_H
#define GT2_REGULATION_REGULATOR_ENRICHMENT_RESULT_H

#include <string>
#include <vector>

#include <genetrail2/core/ConfidenceInterval.h>
#include <genetrail2/core/Statistic.h>

namespace GeneTrail
{
struct RegulatorEnrichmentResult
{
	std::string name = "";
	size_t hits = 0;
	std::vector<double> scores;
	double score = 0.0;
	double mean_correlation = 0.0;
	double p_value = 1.0;
	double corrected_p_value = 1.0;
	std::string targets;

	void addScore(double score) { scores.emplace_back(score); }

	virtual std::string header() const
	{
		std::string header = "#";
		header += "Name\t";
		header += "Hits\t";
		header += "Score\t";
		header += "CI\t";
		header += "Sd\t";
		header += "MAD\t";
		header += "P-value\t";
		header += "Mean(correlation)\n";
		return header;
	}

	template <typename Algorithm>
	void normalize_scores(Algorithm algorithm, size_t n)
	{
		score = algorithm.normalize_score(n, hits, score);
		for(size_t i = 0; i < scores.size(); ++i) {
			scores[i] = algorithm.normalize_score(n, hits, scores[i]);
		}
	}

	void serialize(std::ostream& strm, double alpha,
	               const std::string& ci_method)
	{
		auto sd = statistic::sd<double>(scores.begin(), scores.end());
		auto mad = statistic::mean_absolute_deviation<double>(scores.begin(),
		                                                      scores.end());
		auto ci =
		    confidence_interval<double, std::vector<double>::iterator>::compute(
		        ci_method, scores.begin(), scores.end(), score, alpha);
		strm << name << '\t' << hits << '\t' << score << '\t' << "["
		     << std::get<0>(ci) << ',' << std::get<1>(ci) << "]\t" << sd << '\t'
		     << mad << '\t' << corrected_p_value << '\t' << mean_correlation
		     << '\n';
	}

	template <typename Writer>
	void serializeJSON(Writer& writer, double alpha,
	                   const std::string& ci_method)
	{
		auto sd = statistic::sd<double>(scores.begin(), scores.end());
		auto mad = statistic::mean_absolute_deviation<double>(scores.begin(),
		                                                      scores.end());
		auto ci =
		    confidence_interval<double, std::vector<double>::iterator>::compute(
		        ci_method, scores.begin(), scores.end(), score, alpha);
		writer.StartObject();

		writer.String("regulator");
		writer.String(name.c_str());

		writer.String("score");
		writer.Double(score);

		writer.String("hits");
		writer.Double(hits);

		writer.String("confidenceInterval");
		writer.StartArray();
		writer.Double(std::get<0>(ci));
		writer.Double(std::get<1>(ci));
		writer.EndArray();

		writer.String("stDev");
		writer.Double(sd);

		writer.String("mad");
		writer.Double(mad);

		writer.String("pValue");
		writer.Double(p_value);

		writer.String("correctedPValue");
		writer.Double(corrected_p_value);

		writer.String("meanCorrelation");
		writer.Double(mean_correlation);

		writer.EndObject();
	}
};
}

#endif // GT2_REGULATION_REGULATOR_ENRICHMENT_RESULT_H
