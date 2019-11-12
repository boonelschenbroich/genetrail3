#include <genetrail2/core/Exception.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixWriter.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/ORAPreprocessor.h>

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string reference, out;
size_t max_category_size_, test;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()
		("reference, r", bpo::value<std::string>(&reference)->required(), "Reference set.")
		("test-set-size, t", bpo::value<size_t>(&test)->required(), "Test set size.")
		("max-category-size, m", bpo::value<size_t>(&max_category_size_)->default_value(1000), "Maximum category size.")
		("out, o", bpo::value<std::string>(&out)->required(), "Output file.");

	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),vm);
		bpo::notify(vm);
	}
	catch(bpo::error& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	return true;
}

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv)) {
		return -1;
	}

	GeneSet reference_set;

	GeneSetReader reader;
	try
	{
		reference_set = reader.readGeneList(reference);
	}
	catch(IOError& exn)
	{
		std::cerr << "ERROR: Failed to read reference set. Reason: " << exn.what() << std::endl;
		return -1;
	}

	ORAPreprocessor proc(reference_set.size(), test, max_category_size_);
	DenseMatrixWriter writer;

    std::ofstream file;
	file.open(out);
	writer.writeBinary(file, proc.getPValues());
	file.close();

	return 0;
}
