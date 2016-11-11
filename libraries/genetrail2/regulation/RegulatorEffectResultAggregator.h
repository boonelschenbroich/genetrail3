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

#ifndef GT2_REGULATION_REGULATOR_EFFECT_RESULT_AGGREGATOR_H
#define GT2_REGULATION_REGULATOR_EFFECT_RESULT_AGGREGATOR_H

#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <fstream>
#include <algorithm>

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/prettywriter.h>

#include "RegulatorEffectResult.h"

#include <genetrail2/core/Exception.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/macros.h>

namespace GeneTrail
{

struct GT2_EXPORT AggregatedRegulatorEffectResult
{
	std::string name = "";
	std::vector<size_t> ranks;
	size_t rank_sum = 0;
	std::vector<double> p_values;
	double aggregated_p_value = 0.0;
	double corrected_p_value = 0.0;
	std::vector<double> mean_correlations;
};

struct GT2_EXPORT sum_rank_aggregator
{
	size_t operator()(const std::vector<size_t>& values)
	{
		return std::accumulate(values.begin(), values.end(), size_t());
	}
};

struct GT2_EXPORT max_pvalue_aggregator
{
	double operator()(const std::vector<double>& values)
	{
		return *std::max_element(values.begin(), values.end());
	}
};

struct GT2_EXPORT second_order_pvalue_aggregator
{
	double operator()(const std::vector<double>& values)
	{
		std::vector<double> tmp(values);
		std::nth_element(tmp.begin(), tmp.begin() + 1, tmp.end());
		return tmp[1];
	}
};

struct GT2_EXPORT fisher_pvalue_aggregator
{
	double operator()(const std::vector<double>& pvalues)
	{
		return pvalue::fisher(pvalues, pvalue::identity());
	}
};

struct GT2_EXPORT stouffer_pvalue_aggregator
{
	double operator()(const std::vector<double>& pvalues)
	{
		std::vector<double> weights(pvalues.size(), 1.0 / pvalues.size());
		return pvalue::stouffer(pvalues, weights, pvalue::identity());
	}
};

class GT2_EXPORT RegulatorEffectResultAggregator
{
  public:
	RegulatorEffectResultAggregator() {}

	/**
	 *
	 */
	void read_results(const std::string& fname);

	template <typename Aggregator> void aggregatePValues(Aggregator aggregator)
	{
		for(auto& entry : aggregated_results_) {
			entry.second.aggregated_p_value = aggregator(entry.second.p_values);
		}
	}

	void aggregatePValues(const std::string& method);
	
	template <typename Aggregator> void aggregateRanks(Aggregator aggregator)
	{
		for(auto& entry : aggregated_results_) {
			entry.second.rank_sum = aggregator(entry.second.ranks);
		}
	}
	
	void aggregateRanks(const std::string& method);

	void adjustPValues(const std::string& method);

	void write(const std::string& fname);

  protected:
	template <typename Writer>
	void serializeJSON(Writer& writer,
	                   AggregatedRegulatorEffectResult result)
	{
		writer.StartObject();

		writer.String("regulator");
		writer.String(result.name.c_str());

		writer.String("rankSum");
		writer.Int(result.rank_sum);
		
		writer.String("ranks");
		writer.StartArray();
		for(const auto& r : result.ranks) {
			writer.Int(r);
		}
		writer.EndArray();

		writer.String("aggregatedPValue");
		writer.Double(result.aggregated_p_value);

		writer.String("correctedPValue");
		writer.Double(result.corrected_p_value);

		writer.String("pValues");
		writer.StartArray();
		for(const auto& p : result.p_values) {
			writer.Double(p);
		}
		writer.EndArray();
		writer.EndObject();
	}

	/**
	 *
	 */
	template <typename Value> void parse_result(Value& v)
	{
		if(!v.HasMember("regulator")) {
			throw IOError("Found result without name.");
		}

		std::string name = v["regulator"].GetString();

		auto it = aggregated_results_.find(name);
		if(it == aggregated_results_.end()) {
			AggregatedRegulatorEffectResult result;
			result.name = name;
			it = aggregated_results_.emplace(name,result).first;
		}

		if(!v.HasMember("rank")) {
			throw IOError("Found result without rank.");
		}

		it->second.ranks.emplace_back(v["rank"].GetInt());
		//it->second.rank_sum += v["rank"].GetInt();

		if(!v.HasMember("pValue")) {
			throw IOError("Found result without number of pValue.");
		}

		it->second.p_values.emplace_back(v["pValue"].GetDouble());
	}

	/**
	 *
	 */
	void parse_results(const std::string& line);

	private:
	size_t numberOfResults;
	std::map<std::string, AggregatedRegulatorEffectResult> aggregated_results_;
};
}

#endif // GT2_REGULATION_REGULATOR_ENRICHMENT_RESULT_AGGREGATOR_H
