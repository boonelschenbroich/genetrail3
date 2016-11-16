#ifndef REGULATION_COMMON_H
#define REGULATION_COMMON_H

#include <boost/program_options.hpp>

#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/regulation/RegulatorEffectResult.h>

#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/prettywriter.h>

#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

namespace bpo = boost::program_options;
using namespace GeneTrail;

struct MatrixNameDatabase
{
	MatrixNameDatabase(DenseMatrix* matrix) : matrix_(matrix) {}

	std::string operator()(size_t index) { return matrix_->rowNames()[index]; }

	size_t operator()(std::string name) { return matrix_->rowIndex(name); }

	size_t size() { return matrix_->rows(); }

	DenseMatrix* matrix_;
};

bool checkIfFileExists(bpo::options_description& desc, std::string fname);

std::vector<size_t>
translate_test_set(const DenseMatrix* matrix,
                   const std::vector<std::string>& test_set_names);

/**
 * Saves all RegulatorEffectResults to the given file.
 *
 * @param out_ File to which the results should be written.
 * @param confidence_ Confidence level for the bootstrap confidence intervals.
 */
void writeResults(std::vector<RegulatorEffectResult>& results_,
                  const std::string& out_);

/**
 * Saves all RegulatorEffectResults to the given file.
 *
 * @param out_ File to which the results should be written.
 * @param confidence_ Confidence level for the bootstrap confidence intervals.
 */
void writeJSONResults(std::vector<RegulatorEffectResult>& results_,
                      const std::string& out_);

void write(std::vector<RegulatorEffectResult>& results, const std::string& out, bool json);

/**
 * Adjusts and updates all p-values to account for the multiple testing problem.
 *
 * @param method Name of the correction method that should be performed.
 */
void adjustPValues(std::vector<RegulatorEffectResult>& results_,
                   const std::string& method);

#endif // REGULATION_COMMON_H
