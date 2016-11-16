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

#ifndef GT2_REGULATION_REGULATOR_EFFECT_RESULT_H
#define GT2_REGULATION_REGULATOR_EFFECT_RESULT_H

#include <string>
#include <vector>
#include <tuple>

#include <genetrail2/core/ConfidenceInterval.h>
#include <genetrail2/core/Statistic.h>
#include <genetrail2/core/macros.h>

namespace GeneTrail
{
struct GT2_EXPORT RegulatorEffectResult
{
	std::string name = "";
	size_t hits = 0;
	double expected_hits = 0;
	size_t rank = 0;
	std::tuple<double, double> ci;
	double sd = 0.0;
	double mad = 0.0;
	std::vector<double> scores;
	double score = 0.0;
	double mean_correlation = 0.0;
	double p_value = 1.0;
	double corrected_p_value = 1.0;
	bool skip = true;

	void addScore(double score) { scores.emplace_back(score); }

	virtual std::string header() const
	{
		std::string header = "#";
		header += "Name\t";
		header += "Rank\t";
		header += "Hits\t";
		header += "Score\t";
		if(!skip) {
			header += "CI\t";
			header += "Sd\t";
			header += "MAD\t";
		}
		header += "P-value\t";
		if(!skip) {
			header += "Mean(correlation)\n";
		}
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

	void serialize(std::ostream& strm)
	{
		strm << name << '\t' << rank << '\t' << hits << '\t' << score << '\t';
		if(!skip) {
			strm << "[" << std::get<0>(ci) << ',' << std::get<1>(ci) << "]\t"
			     << sd << '\t' << mad;
		}
		strm << '\t' << corrected_p_value;
		if(!skip) {
			strm << '\t' << mean_correlation;
		}
		strm << '\n';
	}

	void calculate_bootstrap_parameters(double alpha,
	                                    const std::string& ci_method)
	{
		sd = statistic::sd<double>(scores.begin(), scores.end());
		mad = statistic::mean_absolute_deviation<double>(scores.begin(),
		                                                 scores.end());
		ci =
		    confidence_interval<double, std::vector<double>::iterator>::compute(
		        ci_method, scores.begin(), scores.end(), score, alpha);
		skip = false;
	}

	template <typename Writer> void serializeJSON(Writer& writer)
	{
		writer.StartObject();

		writer.String("regulator");
		writer.String(name.c_str());

		writer.String("rank");
		writer.Int(rank);

		writer.String("hits");
		writer.Int(hits);

		writer.String("score");
		writer.Double(score);

		writer.String("expected_hits");
		writer.Double(expected_hits);

		if(!skip) {
			writer.String("confidenceInterval");
			writer.StartArray();
			writer.Double(std::get<0>(ci));
			writer.Double(std::get<1>(ci));
			writer.EndArray();

			writer.String("stDev");
			writer.Double(sd);

			writer.String("mad");
			writer.Double(mad);
		}

		writer.String("pValue");
		writer.Double(p_value);

		writer.String("correctedPValue");
		writer.Double(corrected_p_value);

		if(!skip) {
			writer.String("meanCorrelation");
			writer.Double(mean_correlation);
		}

		writer.EndObject();
	}

	template <typename Value> void deserializeJSON(const Value& result)
	{
		if(!result.IsObject()) {
			throw IOError("Entry " + std::to_string(result) +
			              " is not an object.");
		}

		if(!result.HasMember("regulator")) {
			throw IOError("Found result without name.");
		}

		name = result["regulator"].GetString();

		if(!result.HasMember("rank")) {
			throw IOError("Found result without rank.");
		}

		rank = result["rank"].GetInt();

		if(!result.HasMember("hits")) {
			throw IOError("Found result without number of hits.");
		}

		hits = result["hits"].GetInt();

		if(!result.HasMember("confidenceInterval")) {
			throw IOError("Found result without confidence interval.");
		}

		auto& a = result["confidenceInterval"].GetArray();
		ci = std::make_tuple(a[0].GetDouble(), a[1].GetDouble());

		if(!result.HasMember("stDev")) {
			throw IOError("Found result without number of stDev.");
		}

		sd = result["stDev"].GetDouble();

		if(!result.HasMember("mad")) {
			throw IOError("Found result without number of mad.");
		}

		mad = result["mad"].GetDouble();

		if(!result.HasMember("pValue")) {
			throw IOError("Found result without number of pValue.");
		}

		p_value = result["pValue"].GetDouble();

		if(!result.HasMember("correctedPValue")) {
			throw IOError("Found result without number of correctedPValue.");
		}

		corrected_p_value = result["correctedPValue"].GetDouble();

		if(result.HasMember("meanCorrelation")) {
			mean_correlation = result["meanCorrelation"].GetDouble();
			;
		}
	}
};
}

#endif // GT2_REGULATION_REGULATOR_EFFECT_RESULT_H
