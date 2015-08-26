#include "common.h"

#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/EnrichmentAlgorithm.h>
#include <genetrail2/core/PermutationTest.h>
#include <genetrail2/core/PValue.h>
#include <genetrail2/core/TextFile.h>

#include <boost/filesystem.hpp>

#include <algorithm>
#include <fstream>

void addCommonCLIArgs(bpo::options_description& desc, Params& p)
{
	using namespace bpo;

	desc.add_options()
		("help,h", "Display this message")
		("significance,t", value(&p.significance)->default_value(0.01), "The critical value for rejecting the H0 hypothesis.")
		("categories,c",   value(&p.categories)->required(), "A .gmt file containing the categories to be tested.")
		("scores,s",       value(&p.scores), "A whitespace seperated file containing identifier and scores.")
		("minimum,n",      value(&p.minimum)->default_value(0), "Minimum number of genes allowed in categories.")
		("maximum,x",      value(&p.maximum)->default_value(1000), "Maximum number of genes allowed in categories.")
		("output,o",       value(&p.out), "Output prefix for text files.")
		("adjustment,a",   value(&p.adjustment)->default_value("no"), "P-value adjustment method for multiple testing.")
		("adjust_separately,u", value(&p.adjustSeparately)->default_value(false)->zero_tokens(), "Indicates if databases are adjusted separatly or combined.")
		("pvalue_strategy,m",   value(&p.pValueMode)->default_value(PValueMode::RowWise, "row-wise"), "How should p-values be computed. Possible choices are 'row-wise', 'column-wise', and 'restandardize'")
		("permutations,p",      value(&p.numPermutations)->default_value(1000000), "If p-values are computed using a permutation test, how many permutations should be used.")
		("data_matrix_path,d",  value(&p.dataMatrixPath), "If p-values are computed not using the row-wise strategy, a data matrix must be specified from which scores can be computed.")
		("groups,g",            value(&p.groups), "If p-values are computed not using the 'row-wise' strategy, this file determines the samples used for sample and reference group.")
		("scoring_method,h",    value(&p.scoringMethod), "If p-values are computed not using the 'row-wise' strategy, a scoring method must be provided with which scores should be computed.")
		("seed,r",              value(&p.randomSeed), "If p-values are computed using a permutation test, this option can be used for providing a seed for the random number generator.")
	;
}

bool checkCLIArgs(const Params& p) {
	if(p.pValueMode == PValueMode::ColumnWise || p.pValueMode == PValueMode::Restandardize) {
		if(p.scoringMethod == "") {
			std::cerr << "You must specify a valid scoring method for this pvalue strategy" << std::endl;
			return false;
		}

		if(p.dataMatrixPath == "") {
			std::cerr << "You must specify a data matrix from which scores can be computed" << std::endl;
			return false;
		}

		if(!boost::filesystem::exists(p.dataMatrixPath)) {
			std::cerr << "The file '" << p.dataMatrixPath << "' does not exist." << std::endl;
			return false;
		}

		if(p.groups == "") {
			std::cerr << "You must specify a file containing the groups of the data matrix that are used during score computation." << std::endl;
			return false;
		}

		if(!boost::filesystem::exists(p.groups)) {
			std::cerr << "The file '" << p.groups << "' does not exist." << std::endl;
			return false;
		}
	}

	if(p.scores != "" && !boost::filesystem::exists(p.scores)) {
		std::cerr << "ERROR: the file '" << p.scores << "' does not exist." << std::endl;
	}

	if(p.identifier != "" && !boost::filesystem::exists(p.identifier)) {
		std::cerr << "ERROR: the file '" << p.scores << "' does not exist." << std::endl;
	}

	return true;
}

CategoryList getCategoryList(const std::string& catfile_list)
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

void readTestSet(GeneSet& test_set, const Params& p)
{
	GeneSetReader reader;
	if(p.scores != "" && p.identifier != "") {
		throw GeneTrail::IOError("Too many input files specified.");
	} else if(p.scores != "") {
		test_set = reader.readScoringFile(p.scores);
	} else if(p.identifier != "") {
		test_set = reader.readGeneList(p.identifier);
	} else {
		throw GeneTrail::IOError("No input file specified.");
	}
}

PValueList resultVector(const Results& results)
{
	PValueList result;
	result.reserve(results.size());
	for(const auto& jt : results) {
		result.push_back(
		    std::make_pair(jt.first, jt.second->pvalue.convert_to<double>()));
	}
	return result;
}

PValueList resultVector(const AllResults& results)
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

std::tuple<bool, size_t, std::string>
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

	return std::make_tuple(p.minimum <= c.size() && c.size() <= p.maximum, subset.size(), std::move(entries));
}

void writeFiles(const std::string& output_dir, const AllResults& all_results)
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
		output << database.second.begin()->second->header() << std::endl;
		for(const auto& ele : database.second) {
			output << *ele.second << std::endl;
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
		// TODO: Add single category feature
		cat_list = getCategoryList(p.categories);
	} catch(IOError& exn) {
		std::cerr << "ERROR: Failed to read categories. Reason: " << exn.what()
		          << std::endl;
		return -1;
	}
	return 0;
}

void updatePValues(Results& results, const PValueList& pvalues)
{
	for(const auto& it : pvalues) {
		results[it.first]->pvalue = it.second;
	}
}

void updatePValues(AllResults& results, const PValueList& pvalues)
{
	for(unsigned int i = 0; i < pvalues.size(); ++i) {
		std::vector<std::string> s;
		boost::split(s, pvalues[i].first, boost::is_any_of("\t"));
		results[s[0]][s[1]]->pvalue = pvalues[i].second;
	}
}

AllResults compute(Scores& test_set, CategoryList& cat_list,
                   EnrichmentAlgorithmPtr& algorithm, const Params& p)
{
	AllResults name_to_cat_results;
	for(const auto& cat : cat_list) {
		try {
			GMTFile input(cat.second);

			if(!input) {
				std::cerr << "WARNING: Could not open database " + cat.first +
				                 " for reading! Skipping database."
				          << std::endl;
				continue;
			}

			Results name_to_result;
			while(input) {
				auto c = std::make_shared<Category>(input.read());
				std::cout << "INFO: Processing - " << cat.first << " - "
				          << c->name() << std::endl;
				auto processed = processCategory(*c, test_set, p);

				std::shared_ptr<EnrichmentResult> result;

				if(!std::get<0>(processed)) {
					continue;
				}

				if(!algorithm->canUseCategory(*c, std::get<1>(processed))) {
					result = std::make_shared<EnrichmentResult>(c);
				} else {
					result = algorithm->computeEnrichment(c);
				}

				result->hits = std::get<1>(processed);
				result->info = std::move(std::get<2>(processed));

				name_to_result.emplace(c->name(), std::move(result));
			}
			name_to_cat_results.emplace(cat.first, std::move(name_to_result));
		} catch(IOError& exn) {
			std::cerr << "WARNING: Could not process category file "
			          << cat.first << "! " << std::endl;
		}
	}

	return name_to_cat_results;
}

void adjustCombined(AllResults& all_results, const Params& p)
{
	auto results = resultVector(all_results);
	results = pvalue<double>::adjustPValues(results, p.adjustment);
	updatePValues(all_results, results);
}

void adjustSeparately(AllResults& all_results, const Params& p)
{
	for(auto& results_it : all_results) {
		auto results = resultVector(results_it.second);
		results = pvalue<double>::adjustPValues(results, p.adjustment);
		updatePValues(results_it.second, results);
	}
}

void computeRowWisePValues(const EnrichmentAlgorithmPtr& algorithm,
                           EnrichmentResults& results, const Scores& scores, const Params& p)
{
	using Test = RowPermutationTest<double>;

	auto test =
	    algorithm->supportsIndices()
	        ? Test::IndexBased(scores, p.numPermutations, p.randomSeed)
	        : Test::CategoryBased(scores, p.numPermutations, p.randomSeed);

	test->computePValue(algorithm, results);
}

void removeUnusedColumns(DenseMatrix& data,
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

void removeUnusedRows(DenseMatrix& data, const Scores& scores) {
	std::vector<Matrix::index_type> rowsToRemove;

	Matrix::index_type i = 0;
	for(const auto& rowName : data.rowNames()) {
		if(!scores.contains(rowName)) {
			rowsToRemove.push_back(i);
		}
	}

	data.removeRows(rowsToRemove);

	if(data.rows() != scores.size()) {
		throw std::string("Input scores and matrix are incompatible");
	}
}

void sortMatrixRows(DenseMatrix& data) {
	std::vector<size_t> matrix_indices(data.rows());
	EntityDatabase::global->transform(data.rowNames(), matrix_indices.begin());

	auto permutation = sort_permutation(matrix_indices.begin(), matrix_indices.end(), std::less<size_t>());

	std::vector<DenseMatrix::index_type> tmp(permutation.size());

	std::copy(permutation.begin(), permutation.end(), tmp.begin());

	data.shuffleRows(tmp);
}

void computeColumnWisePValues(const EnrichmentAlgorithmPtr& algorithm,
                              EnrichmentResults& results, const Scores& scores, const Params& p)
{
	std::ifstream input(p.dataMatrixPath);
	DenseMatrixReader matrixReader;
	DenseMatrix data = matrixReader.read(input);

	TextFile t(p.groups, ",");

	auto referenceGroup = t.read();
	auto sampleGroup = t.read();

	removeUnusedColumns(data, referenceGroup, sampleGroup);
	removeUnusedRows(data, scores);
	sortMatrixRows(data);

	ColumnPermutationTest<double> test(data, p.numPermutations,
	                                   referenceGroup.size(), p.scoringMethod, p.randomSeed);

	test.computePValue(algorithm, results);
}

void
computeRestandardizationPValues(const EnrichmentAlgorithmPtr& algorithm)
{
	throw NotImplemented(__FILE__, __LINE__,
	                     "void computeRestandardizationPValues(const "
	                     "EnrichmentAlgorithmPtr&, "
	                     "PermutationTest<double>::TestResults&)");
}

void computePValues(EnrichmentAlgorithmPtr& algorithm,
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
			computeColumnWisePValues(algorithm, results, scores, p);
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

	// Checks how they should be adjusted
	if(p.adjustSeparately) {
		adjustSeparately(name_to_cat_results, p);
	} else {
		adjustCombined(name_to_cat_results, p);
	}

	writeFiles(p.out, name_to_cat_results);
}

namespace GeneTrail
{
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
}
