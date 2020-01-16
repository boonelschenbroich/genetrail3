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
#include "common.h"

#include "EnrichmentAlgorithm.h"
#include "EnrichmentResult.h"
#include "Parameters.h"
#include "PermutationTest.h"

#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/TextFile.h>

#include <boost/filesystem.hpp>

#include <algorithm>
#include <fstream>

static CategoryList getCategoryList(const std::string& catfile_list)
{
	CategoryList categories;

	std::ifstream input(catfile_list);
	if(!input) {
		throw GeneTrail::IOError("File (" + catfile_list +
		                         ") is not open for reading");
	}

	for(std::string line; getline(input, line);) {
		std::vector<std::string> sline(2);
		boost::split(sline, line, boost::is_any_of(" \t"));
		if(sline.size() == 2) {
			boost::trim(sline[0]);
			boost::trim(sline[1]);
			categories.emplace_back(sline[0], sline[1]);
		} else {
			throw GeneTrail::IOError("Wrong file format.");
		}
	}

	return categories;
}

static void readTestSet(GeneSet& test_set, const Params& p)
{
	GeneSetReader reader;
	if(p.scores() != "" && p.identifier() != "") {
		throw GeneTrail::IOError("Too many input files specified.");
	} else if(p.scores() != "") {
		test_set = reader.readScoringFile(p.scores());
	} else if(p.identifier() != "") {
		test_set = reader.readGeneList(p.identifier());
	} else {
		throw GeneTrail::IOError("No input file specified.");
	}
}

static PValueList resultVector(const Results& results)
{
	PValueList result;
	result.reserve(results.size());
	for(const auto& jt : results) {
		result.push_back(
		    std::make_pair(jt.first, jt.second->pvalue.convert_to<double>()));
	}
	return result;
}

static PValueList resultVector(const AllResults& results)
{
	PValueList result;
	result.reserve(results.size());
	for(const auto& it : results) {
		for(const auto& jt : it.second) {
			result.push_back(
			    std::make_pair(it.first + "\t" + jt.first,
			                   jt.second->pvalue.convert_to<double>()));
		}
	}
	return result;
}

static std::tuple<bool, size_t, std::string>
processCategory(const Category& c, const Scores& test_set, const Params& p)
{
	Scores subset = test_set.subset(c);
	subset.sortByName();

	std::string entries;
	for(const auto& entry : subset.names()) {
		entries += entry + ',';
	}

	if(!entries.empty()) {
		entries.resize(entries.size() - 1);
	}

	return std::make_tuple(p.minimum <= subset.size() && subset.size() <= p.maximum, subset.size(), std::move(entries));
}

static void writeFiles(const std::string& output_dir, const AllResults& all_results, const bool reducedOutput)
{

	size_t size = 0;
	for(auto& results_it : all_results) {
		size += results_it.second.size();
	}

	for(const auto& database : all_results) {
		std::ofstream output(output_dir + "/" + database.first + ".txt");
		if(database.second.begin() == database.second.end()) {
			output.close();
			continue;
		}
		if(!output) {
			throw GeneTrail::IOError("No input file specified.");
		}
		output << database.second.begin()->second->header(reducedOutput) << std::endl;
		for(const auto& ele : database.second) {
			ele.second->serialize(output, reducedOutput);
			output << std::endl;
		}
		output.close();
	}
}

int init(GeneSet& test_set, CategoryList& cat_list, const Params& p)
{

	try {
		readTestSet(test_set, p);
	} catch(IOError& exn) {
		std::cerr << "ERROR: Failed to read test set. Reason: " << exn.what()
		          << std::endl;
		return -1;
	}

	try {
		if (p.categories() != ""){
			// TODO: Add single category feature
			cat_list = getCategoryList(p.categories());
		} else {
			cat_list.emplace_back(p.category_name(), p.category());
		}
	} catch(IOError& exn) {
		std::cerr << "ERROR: Failed to read categories. Reason: " << exn.what()
		          << std::endl;
		return -1;
	}
	return 0;
}

static void updatePValues(Results& results, const PValueList& pvalues)
{
	for(const auto& it : pvalues) {
		results[it.first]->pvalue = it.second;
	}
}

static void updatePValues(AllResults& results, const PValueList& pvalues)
{
	for(unsigned int i = 0; i < pvalues.size(); ++i) {
		std::vector<std::string> s;
		boost::split(s, pvalues[i].first, boost::is_any_of("\t"));
		results[s[0]][s[1]]->pvalue = pvalues[i].second;
	}
}

static AllResults compute(Scores& test_set, CategoryList& cat_list,
                   EnrichmentAlgorithmPtr& algorithm, const Params& p)
{
	AllResults name_to_cat_results;
	for(const auto& cat : cat_list) {
		try {
			GMTFile input(test_set.db(), cat.second);

			if(!input) {
				std::cerr << "WARNING: Could not open database " + cat.first +
				                 " for reading! Skipping database."
				          << std::endl;
				continue;
			}

			auto category_db = input.read();
			Results name_to_result;
			for(const auto& c : category_db) {
				if(p.verbose) std::cout << "INFO: Processing - " << cat.first << " - " << c.name() << std::endl;
				auto processed = processCategory(c, test_set, p);

				std::shared_ptr<EnrichmentResult> result;

				if(!std::get<0>(processed)) {
					continue;
				}

				// TODO: get rid of this
				auto tmp_cat = std::make_shared<Category>(c);
				if(!algorithm->canUseCategory(c, std::get<1>(processed))) {
					result = std::make_shared<EnrichmentResult>(tmp_cat);
				} else {
					result = algorithm->computeEnrichment(tmp_cat);
				}

				result->hits = std::get<1>(processed);
				result->info = std::move(std::get<2>(processed));

				name_to_result.emplace(c.name(), std::move(result));
			}
			name_to_cat_results.emplace(cat.first, std::move(name_to_result));
		} catch(IOError& exn) {
			std::cerr << "WARNING: Could not process category file "
			          << cat.first << "! " << std::endl;
		}
	}

	return name_to_cat_results;
}

static void adjustCombined(AllResults& all_results, MultipleTestingCorrection correction)
{
	auto results = resultVector(all_results);
	results = pvalue::adjustPValues(results, pvalue::get_second(), correction);
	updatePValues(all_results, results);
}

static void adjustSeparately(AllResults& all_results, MultipleTestingCorrection correction)
{
	for(auto& results_it : all_results) {
		auto results = resultVector(results_it.second);
		results = pvalue::adjustPValues(results, pvalue::get_second(), correction);
		updatePValues(results_it.second, results);
	}
}

static void computeRowWisePValues(const EnrichmentAlgorithmPtr& algorithm,
                           EnrichmentResults& results, const Scores& scores, const Params& p)
{
	using Test = RowPermutationTest<double>;

	auto test =
	    algorithm->supportsIndices()
	        ? Test::IndexBased(scores, p.numPermutations, p.randomSeed)
	        : Test::CategoryBased(scores, p.numPermutations, p.randomSeed);

	test->computePValue(algorithm, results);
}

static void removeUnusedColumns(DenseMatrix& data,
                         const std::vector<std::string>& ref,
                         const std::vector<std::string>& sample)
{
	std::vector<DenseMatrix::index_type> colsToDelete;
	DenseMatrix::index_type i = 0;

	auto colIsNotUsed = [&ref, &sample](const std::string& name) {
		return std::find(ref.begin(), ref.end(), name) == ref.end() &&
		       std::find(sample.begin(), sample.end(), name) == sample.end();
	};

	for(auto&& colName : data.colNames()) {
		if(colIsNotUsed(colName)) {
			colsToDelete.push_back(i);
		}
		++i;
	}

	data.removeCols(colsToDelete);
}

static void removeUnusedRows(DenseMatrix& data, const Scores& scores) {
	std::vector<Matrix::index_type> rowsToRemove;

	Matrix::index_type i = 0;
	for(const auto& rowName : data.rowNames()) {
		if(!scores.contains(rowName)) {
			rowsToRemove.push_back(i);
		}
		++i;
	}

	data.removeRows(rowsToRemove);

	if(data.rows() != scores.size()) {
		throw std::string("Input scores and matrix are incompatible");
	}
}

static void sortMatrixRows(DenseMatrix& data, const EntityDatabase* db) {
	std::vector<size_t> matrix_indices(data.rows());
	db->transform(data.rowNames(), matrix_indices.begin());

	auto permutation = sort_permutation(matrix_indices.begin(), matrix_indices.end(), std::less<size_t>());

	std::vector<DenseMatrix::index_type> tmp(permutation.size());

	std::copy(permutation.begin(), permutation.end(), tmp.begin());

	data.shuffleRows(tmp);
}

static void computeColumnWisePValues(const EnrichmentAlgorithmPtr& algorithm,
                                     EnrichmentResults& results,
                                     const Scores& scores, const Params& p,
                                     const EntityDatabase* db)
{
	std::ifstream input(p.dataMatrixPath());
	DenseMatrixReader matrixReader;
	DenseMatrix data = matrixReader.read(input);

	TextFile t(p.groups(), ",");

	auto referenceGroup = t.read();
	auto sampleGroup = t.read();

	removeUnusedColumns(data, referenceGroup, sampleGroup);
	removeUnusedRows(data, scores);
	sortMatrixRows(data, scores.db().get());

	if(p.adjustment && p.adjustment == MultipleTestingCorrection::GSEA) {
		KSColumnPermutationTest<double> test(
		    data, p.numPermutations, referenceGroup.size(),
		    p.scoringMethod.get(), p.randomSeed, db);
		test.computePValue(algorithm, results);
	} else {
		ColumnPermutationTest<double> test(
		    data, p.numPermutations, referenceGroup.size(),
		    p.scoringMethod.get(), p.randomSeed, db);
		test.computePValue(algorithm, results);
	}
}

static void
computeRestandardizationPValues(const EnrichmentAlgorithmPtr&)
{
	throw NotImplemented(__FILE__, __LINE__,
	                     "void computeRestandardizationPValues(const "
	                     "EnrichmentAlgorithmPtr&, "
	                     "PermutationTest<double>::TestResults&)");
}

static void computePValues(EnrichmentAlgorithmPtr& algorithm,
                    const AllResults& name_to_cat_results, const Scores& scores, const Params& p)
{
	EnrichmentResults results;

	// TODO: This is really inefficient and clumsy
	//      we should have a better way of getting the necessary stuff
	//      together for pvalue computation.
	for(const auto& db : name_to_cat_results) {
		for(const auto& result : db.second) {
			results.push_back(result.second);
		}
	}

	switch(algorithm->pValueMode()) {
		case PValueMode::RowWise:
			computeRowWisePValues(algorithm, results, scores, p);
			break;
                case PValueMode::ColumnWise:
                        computeColumnWisePValues(algorithm, results, scores, p, scores.db().get());
                        break;
                case PValueMode::Restandardize:
                        computeRestandardizationPValues(algorithm);
                        break;
        }
}

void run(Scores& test_set, CategoryList& cat_list,
         EnrichmentAlgorithmPtr& algorithm, const Params& p, bool computePValue)
{
        test_set.sortByIndex();

        AllResults name_to_cat_results(compute(test_set, cat_list, algorithm, p));
        if(computePValue && !algorithm->pValuesComputed()) {
                computePValues(algorithm, name_to_cat_results, test_set, p);
        }

        if(p.adjustment && boost::get(p.adjustment) != MultipleTestingCorrection::GSEA) {
                // Checks how they should be adjusted
                if(p.adjustSeparately) {
                        adjustSeparately(name_to_cat_results, p.adjustment.get());
                } else {
                        adjustCombined(name_to_cat_results, p.adjustment.get());
                }
        }

        writeFiles(p.out(), name_to_cat_results, p.reducedOutput);
}

