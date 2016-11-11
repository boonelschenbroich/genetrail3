#include "RegulatorEffectResultAggregator.h"

namespace GeneTrail
{
void RegulatorEffectResultAggregator::read_results(const std::string& fname)
{
	std::ifstream infile(fname);
	std::string line;
	numberOfResults = 0;
	while(std::getline(infile, line)) {
		numberOfResults += 1;
		parse_results(line);
	}
}

void
RegulatorEffectResultAggregator::adjustPValues(const std::string& method)
{
	std::vector<std::pair<std::string, double>> p_values;
	p_values.reserve(aggregated_results_.size());
	for(auto& r : aggregated_results_) {
		p_values.emplace_back(
		    std::make_pair(r.first, r.second.aggregated_p_value));
	}

	// Adjust p-values
	p_values = pvalue::adjustPValues(p_values, pvalue::get_second(),
	                                 pvalue::getCorrectionMethod(method).get());

	// Update p-values
	for(const auto& pair : p_values) {
		aggregated_results_[pair.first].corrected_p_value = pair.second;
	}
}

void
RegulatorEffectResultAggregator::aggregatePValues(const std::string& method)
{
	/**
	 * This needs to be done since not every regulator might be detected by 
	 * every REAGGAE run! 
	 */
	for(auto& entry : aggregated_results_) {
		for(size_t i = entry.second.p_values.size(); i < numberOfResults; ++i){
			entry.second.p_values.emplace_back(1.0);
		}
	}
	
	if(method == "max") {
		RegulatorEffectResultAggregator::aggregatePValues(max_pvalue_aggregator());
	} else if(method == "second-order") {
		RegulatorEffectResultAggregator::aggregatePValues(
		    second_order_pvalue_aggregator());
	} else if(method == "fisher") {
		RegulatorEffectResultAggregator::aggregatePValues(
		    fisher_pvalue_aggregator());
	} else if(method == "stouffer") {
		RegulatorEffectResultAggregator::aggregatePValues(
		    stouffer_pvalue_aggregator());
	} else {
		throw NotImplemented(
		    __FILE__, __LINE__,
		    "The passed aggregation method is not yet implemented");
	}
}

void
RegulatorEffectResultAggregator::aggregateRanks(const std::string& method)
{
	if(method == "sum") {
		
		/**
		 * This needs to be done since not every regulator might be detected by 
		 * every REAGGAE run! 
		 */
		for(auto& entry : aggregated_results_) {
			for(size_t i = entry.second.ranks.size(); i < numberOfResults; ++i){
				entry.second.ranks.emplace_back(aggregated_results_.size());
			}
		}
		
		RegulatorEffectResultAggregator::aggregateRanks(sum_rank_aggregator());
	} else {
		throw NotImplemented(
		    __FILE__, __LINE__,
		    "The passed rank aggregation method is not yet implemented");
	}
}

void RegulatorEffectResultAggregator::parse_results(const std::string& line)
{
	std::ifstream f(line.c_str());
	if(!f.good()) {
		throw IOError("Cannot open result '" + line + "'.");
	}
	rapidjson::IStreamWrapper isw(f);
	rapidjson::Document d;
	d.ParseStream(isw);
	if(!d.IsArray()) {
		throw IOError("Result '" + line + "' is malformed.");
	}
	for(auto& v : d.GetArray()) {
		parse_result(v);
	}
}

void RegulatorEffectResultAggregator::write(const std::string& fname)
{
	std::vector<AggregatedRegulatorEffectResult> results;
	for(auto& r : aggregated_results_) {
		results.emplace_back(r.second);
	}
	std::sort(results.begin(), results.end(),
	          [](const AggregatedRegulatorEffectResult& lhs,
	             const AggregatedRegulatorEffectResult& rhs) {
		return lhs.rank_sum < rhs.rank_sum;
	});
	rapidjson::StringBuffer sb;
	rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(sb);

	writer.StartArray();
	for(AggregatedRegulatorEffectResult& res : results) {
		serializeJSON(writer, res);
	}
	writer.EndArray();

	std::ofstream out;
	out.open(fname);
	if(out.is_open()) {
		out << sb.GetString();
	} else {
		std::cerr << "Could not open file: " << fname << "\n";
	}
	out.close();
}
}
