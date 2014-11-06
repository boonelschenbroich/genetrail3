#include <iostream>
#include <cstdlib>
#include <tuple>

#include <boost/program_options.hpp>

#include <genetrail2/core/MatrixHTest.h>
#include <genetrail2/core/GeneSetWriter.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixSubset.h>
#include <genetrail2/core/TextFile.h>

using namespace GeneTrail;
namespace bpo = boost::program_options;

class EmptyGroup : public std::exception
{
	public:
	EmptyGroup(const std::string& name) noexcept : groupname_(name) {}

	virtual const char* what() const noexcept {
		return (std::string("Group \"") + groupname_ + "\" does not contain any datapoint.").c_str();
	}

	private:
	std::string groupname_;
};

std::string expr1 = "", expr2 = "", output = "", method = "", groups = "";
bool binary = false;
bool additional_colname = false;
bool no_rownames = false;
bool no_colnames = false;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("expression-matrix-1,1", bpo::value<std::string>(&expr1)->required(), "Name of a text file containing expression values as a matrix.")
		("expression-matrix-2,2", bpo::value<std::string>(&expr2), "Name of a text file containing expression values as a matrix.(optional)")
		("output,o", bpo::value<std::string>(&output)->required(), "Name of the output file.")
		("groups,g", bpo::value<std::string>(&groups)->required(), "File containing two lines specifying which rownames belong to which group.")
		("no-row-names,r", bpo::value<bool>(&no_rownames)->default_value(false)->zero_tokens(), "Does the file contain row names.")
		("no-col-names,c", bpo::value<bool>(&no_colnames)->default_value(false)->zero_tokens(), "Does the file contain column names.")
		("add-col-name,a", bpo::value<bool>(&additional_colname)->default_value(false)->zero_tokens(), "File containing two lines specifying which rownames belong to which group.")
		("method,m", bpo::value<std::string>(&method)->required(), "Method used for scoring.");

	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);
		bpo::notify(vm);
	}
	catch(bpo::error& e)
	{
		std::cerr << "ERROR: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	return true;
}

DenseMatrix readDenseMatrix(std::string matrix)
{
	unsigned int opts = DenseMatrixReader::NO_OPTIONS;

	if(!no_rownames) {
		opts |= DenseMatrixReader::READ_ROW_NAMES;
	}

	if(!no_colnames) {
		opts |= DenseMatrixReader::READ_COL_NAMES;
	}

	if(additional_colname) {
		opts |= DenseMatrixReader::ADDITIONAL_COL_NAME;
	}

	DenseMatrixReader reader;

	std::ifstream strm(matrix, std::ios::binary);
	return reader.read(strm, opts);
}

std::vector<unsigned int> getIndices(const DenseMatrix& matrix, const std::vector<std::string>& colnames, const std::string& groupname)
{
	std::vector<unsigned int> indices;
	indices.reserve(colnames.size());
	for(const auto& s : colnames)
	{
		if(matrix.hasCol(s))
		{
			indices.emplace_back(matrix.colIndex(s));
		} else {
			std::cerr << "WARNING: Could not find column \"" + s + "\".\n";
		}
	}

	if(indices.empty()) {
		throw EmptyGroup(groupname);
	}

	return indices;
}

std::tuple<DenseMatrixSubset, DenseMatrixSubset> splitMatrix(DenseMatrix& matrix, const std::vector<std::string>& reference, const std::vector<std::string>& test)
{
	return std::make_tuple(
		DenseMatrixSubset::createColSubset(&matrix, getIndices(matrix, reference, "reference")),
		DenseMatrixSubset::createColSubset(&matrix, getIndices(matrix, test, "test"))
	);
}

DenseMatrix buildDenseMatrix()
{
	auto m1 = readDenseMatrix(expr1);
	if(expr2 != ""){
		auto m2 = readDenseMatrix(expr2);
		m1.cbind(m2);
	}
	return m1;
}

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv))
	{
		return -1;
	}

	if(additional_colname && no_rownames) {
		std::cerr << "Conflicting arguments. Additional colnames can only be "
		             "specified if row names are present!" << std::endl;
		return -2;
	}

	TextFile t(groups, ",", std::set<std::string>());
	DenseMatrix matrix(buildDenseMatrix());
	std::vector<std::string> reference = t.read();
	std::vector<std::string> sample = t.read();

	try {
		auto subset = splitMatrix(matrix, reference, sample);

		MatrixHTest<DenseMatrixSubset> htest;
		auto gene_set = htest.test(std::get<0>(subset), std::get<1>(subset), method);

		GeneSetWriter writer;
		writer.writeScoringFile(gene_set, output);
	} catch(EmptyGroup& e) {
		std::cerr << "ERROR: " << e.what() << "\n";
		return -3;
	}

	return 0;
}
