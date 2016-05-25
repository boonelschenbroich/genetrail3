/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *               2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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
#include "CommandLineInterface.h"

#include "Parameters.h"
#include "PermutationTest.h"

#include <boost/any.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <iostream>

namespace bpo = boost::program_options;

namespace GeneTrail {
	void addCommonCLIArgs(bpo::options_description& desc, Params& p)
	{
		using namespace bpo;

		desc.add_options()
			("help,h", "Display this message")
			("significance,t", value(&p.significance)->default_value(0.01), "The critical value for rejecting the H0 hypothesis.")
			("categories,c",   value(&p.categories_)->required(), "A .gmt file containing the categories to be tested.")
			("scores,s",       value(&p.scores_), "A whitespace separated file containing identifier and scores.")
			("minimum,n",      value(&p.minimum)->default_value(0), "Minimum number of genes allowed in categories.")
			("maximum,x",      value(&p.maximum)->default_value(1000), "Maximum number of genes allowed in categories.")
			("output,o",       value(&p.out_)->required(), "Output prefix for text files.")
			("adjustment,a",   value(&p.adjustment)->default_value(boost::none, "none"), "P-value adjustment method for multiple testing.")
			("adjust_separately,u", value(&p.adjustSeparately)->default_value(false)->zero_tokens(), "Indicates if databases are adjusted separatly or combined.")
			("pvalue_strategy,l",   value(&p.pValueMode)->default_value(PValueMode::RowWise, "row-wise"), "How should p-values be computed. Possible choices are 'row-wise', 'column-wise', and 'restandardize'")
			("permutations,p",      value(&p.numPermutations)->default_value(1000000), "If p-values are computed using a permutation test, how many permutations should be used.")
			("data_matrix_path,d",  value(&p.dataMatrixPath_), "If p-values are computed not using the row-wise strategy, a data matrix must be specified from which scores can be computed.")
			("groups,g",            value(&p.groups_), "If p-values are computed not using the 'row-wise' strategy, this file determines the samples used for sample and reference group.")
			("scoring_method,r",    value(&p.scoringMethod), "If p-values are computed not using the 'row-wise' strategy, a scoring method must be provided with which scores should be computed.")
			("seed,e",              value(&p.randomSeed), "If p-values are computed using a permutation test, this option can be used for providing a seed for the random number generator.")
		;
	}

	bool checkCLIArgs(const Params& p) {
		if(p.pValueMode == PValueMode::ColumnWise || p.pValueMode == PValueMode::Restandardize) {
			if(!p.scoringMethod) {
				std::cerr << "You must specify a valid scoring method for this pvalue strategy" << std::endl;
				return false;
			}

			if(p.dataMatrixPath() == "") {
				std::cerr << "You must specify a data matrix from which scores can be computed" << std::endl;
				return false;
			}

			if(p.groups() == "") {
				std::cerr << "You must specify a file containing the groups of the data matrix that are used during score computation." << std::endl;
				return false;
			}
		}

		if(p.significance <= 0.0 || p.significance > 1.0) {
			std::cerr << "ERROR: the significance level must be > 0.0 and <= 1.0." << std::endl;
			return false;
		}

		if(p.minimum > p.maximum) {
			std::cerr << "WARNING: minimum category size is larger than maximum category size." << std::endl;
		}

		return true;
	}

	void validate(boost::any& v, const std::vector<std::string>& values,
	              PValueMode*, int)
	{
		bpo::validators::check_first_occurrence(v);

		const auto& mode = bpo::validators::get_single_string(values);

		if(mode == "row-wise") {
			v = boost::any(PValueMode::RowWise);
		} else if(mode == "column-wise") {
			v = boost::any(PValueMode::ColumnWise);
		} else if(mode == "restandardize") {
			v = boost::any(PValueMode::Restandardize);
		} else {
			throw bpo::validation_error(
			    bpo::validation_error::invalid_option_value);
		}
	}

	void validate(boost::any& v, const std::vector<std::string>& values,
	              boost::optional<MultipleTestingCorrection>*, int)
	{
		bpo::validators::check_first_occurrence(v);

		const auto& param = bpo::validators::get_single_string(values);

		if(param == "none") {
			v = boost::any(boost::optional<MultipleTestingCorrection>());
			return;
		}

		auto method = pvalue<double>::getCorrectionMethod(param);

		if(!method) {
			throw bpo::validation_error(
			    bpo::validation_error::invalid_option_value);
		}

		v = boost::any(std::move(method));
	}

	void validate(boost::any& v, const std::vector<std::string>& values,
	              boost::optional<MatrixHTests>*, int)
	{
		bpo::validators::check_first_occurrence(v);

		const auto& param = bpo::validators::get_single_string(values);

		MatrixHTestFactory factory;
		auto method = factory.getMethod(param);

		if(!method) {
			throw bpo::validation_error(
			    bpo::validation_error::invalid_option_value);
		}

		v = boost::any(std::move(method));
	}

	void validate(boost::any& v, const std::vector<std::string>& values,
	              FilePath*, int)
	{
		bpo::validators::check_first_occurrence(v);

		const auto& param = bpo::validators::get_single_string(values);

		if(!boost::filesystem::exists(param)) {
			std::cerr << "ERROR: the file '" << param << "' does not exist." << std::endl;
			throw bpo::validation_error(
			    bpo::validation_error::invalid_option_value);
		}

		if(boost::filesystem::is_directory(param)) {
			std::cerr << "ERROR: the path '" << param << "' is a directory." << std::endl;
			throw bpo::validation_error(
			    bpo::validation_error::invalid_option_value);
		}

		v = boost::any(FilePath(param));
	}

	void validate(boost::any& v, const std::vector<std::string>& values,
	              DirectoryPath*, int)
	{
		bpo::validators::check_first_occurrence(v);

		const auto& param = bpo::validators::get_single_string(values);

		if(!boost::filesystem::exists(param)) {
			std::cerr << "ERROR: the directory '" << param << "' does not exist." << std::endl;
			throw bpo::validation_error(
			    bpo::validation_error::invalid_option_value);
		}

		if(!boost::filesystem::is_directory(param)) {
			std::cerr << "ERROR: the path '" << param << "' is not a directory." << std::endl;
			throw bpo::validation_error(
			    bpo::validation_error::invalid_option_value);
		}

		v = boost::any(DirectoryPath(param));
	}

}
