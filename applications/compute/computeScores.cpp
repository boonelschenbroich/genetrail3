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

std::string expr1 = "", expr2 = "", output = "", method = "", groups = "";
bool binary = false;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("expression-matrix-1,1", bpo::value<std::string>(&expr1)->required(), "Name of a text file containing expression values as a matrix.")
		("expression-matrix-2,2", bpo::value<std::string>(&expr2), "Name of a text file containing expression values as a matrix.(optional)")
		("binary,b", bpo::value(&binary)->zero_tokens(), "Flag indicating if the matrix is in binary format.")
		("output,o", bpo::value<std::string>(&output)->required(), "Name of the output file.")
		("groups,g", bpo::value<std::string>(&groups)->required(), "File containing two lines specifying which rownames belog to which group.")
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
	DenseMatrixReader reader;
	if(binary)
	{
		std::ifstream strm(matrix, std::ios::binary);
		return reader.read(strm);
	}else{
		std::ifstream strm(matrix);
		return reader.read(strm);
	}
}

std::tuple<DenseMatrixSubset, DenseMatrixSubset> splitMatrix(DenseMatrix& matrix, const std::vector<std::string>& reference, const std::vector<std::string>& test)
{
	std::vector<unsigned int> indices1;
	indices1.reserve(reference.size());
	for(const auto& s : reference)
	{
		if(matrix.hasCol(s))
		{
			indices1.emplace_back(matrix.colIndex(s));
		}
	}

	auto m1 = DenseMatrixSubset::createColSubset(&matrix, indices1);

	std::vector<unsigned int> indices2;
	indices2.reserve(test.size());
	for(const auto& s : test)
	{
		if(matrix.hasCol(s))
		{
			indices2.emplace_back(matrix.colIndex(s));
		}
	}
	auto m2 = DenseMatrixSubset::createColSubset(&matrix, indices2);
	return std::make_tuple(m1, m2);
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

	TextFile t(groups, ",", std::set<std::string>());
	DenseMatrix matrix(buildDenseMatrix());
	std::vector<std::string> reference = t.read();
	std::vector<std::string> sample = t.read();
	auto subset = splitMatrix(matrix, reference, sample);

	MatrixHTest<DenseMatrixSubset> htest;
	GeneSet<double> gene_set = htest.test(std::get<0>(subset), std::get<1>(subset), method);

	GeneSetWriter<double> writer;
	writer.writeScoringFile(gene_set, output);

	return 0;
}
