#include "common.h"

bool checkIfFileExists(bpo::options_description& desc, std::string fname)
{
	std::ifstream infile(fname);
	if(!infile.good()) {
		std::cerr << "ERROR: File '" + fname + "' does not exist. \n";
		desc.print(std::cerr);
		return false;
	}
	return true;
}

std::vector<size_t>
translate_test_set(const DenseMatrix* matrix,
                   const std::vector<std::string>& test_set_names)
{
	std::vector<size_t> test_set;
	test_set.reserve(test_set_names.size());
	for(const std::string name : test_set_names) {
		test_set.emplace_back(matrix->rowIndex(name));
	}
	return test_set;
}

void writeResults(std::vector<RegulatorEffectResult>& results_,
                  const std::string& out_)
{
	std::vector<RegulatorEffectResult> results;
	for(size_t i = 0; i < results_.size(); ++i) {
		if(results_[i].name == "") {
			continue;
		}
		results.emplace_back(results_[i]);
	}
	std::sort(
	    results.begin(), results.end(),
	    [](const RegulatorEffectResult& lhs, const RegulatorEffectResult& rhs) {
		    return lhs.p_value == rhs.p_value ? abs(lhs.score) > abs(rhs.score)
		                                      : lhs.p_value < rhs.p_value;
		});

	std::ofstream out;
	out.open(out_);
	if(out.is_open()) {
		out << results_[0].header();
		size_t rank = 1;
		for(RegulatorEffectResult result : results) {
			result.rank = rank;
			result.serialize(out);
			++rank;
		}
	} else {
		std::cerr << "Could not open file: " << out_ << "\n";
	}

	out.close();
}

void writeJSONResults(std::vector<RegulatorEffectResult>& results_,
                      const std::string& out_)
{
	std::vector<RegulatorEffectResult> results;
	for(size_t i = 0; i < results_.size(); ++i) {
		if(results_[i].name == "") {
			continue;
		}
		results.emplace_back(results_[i]);
	}
	std::sort(
	    results.begin(), results.end(),
	    [](const RegulatorEffectResult& lhs, const RegulatorEffectResult& rhs) {
		    return lhs.p_value == rhs.p_value ? abs(lhs.score) > abs(rhs.score)
		                                      : lhs.p_value < rhs.p_value;
		});

	rapidjson::StringBuffer sb;
	rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(sb);

	size_t rank = 1;
	writer.StartArray();
	for(RegulatorEffectResult& res : results) {
		res.rank = rank;
		res.serializeJSON(writer);
		++rank;
	}
	writer.EndArray();

	std::ofstream out;
	out.open(out_);
	if(out.is_open()) {
		out << sb.GetString();
	} else {
		std::cerr << "Could not open file: " << out_ << "\n";
	}
	out.close();
}

void write(std::vector<GeneTrail::RegulatorEffectResult>& results,
           const std::string& out, bool json)
{
	if(json) {
		writeJSONResults(results, out);
	} else {
		writeResults(results, out);
	}
}

void adjustPValues(std::vector<RegulatorEffectResult>& results_,
                   const std::string& method)
{
	std::vector<std::pair<size_t, double>> p_values;
	p_values.reserve(results_.size());
	for(size_t i = 0; i < results_.size(); ++i) {
		if(results_[i].name == "") {
			continue;
		}
		p_values.emplace_back(std::make_pair(i, results_[i].p_value));
	}

	// Adjust p-values
	p_values = pvalue::adjustPValues(p_values, pvalue::get_second(),
	                                 pvalue::getCorrectionMethod(method).get());

	// Update p-values
	for(const auto& pair : p_values) {
		results_[pair.first].corrected_p_value = pair.second;
	}
}
