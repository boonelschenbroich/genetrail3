#include "glasso.h"

#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrixWriter.h>

#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

using namespace GeneTrail;
namespace bpo = boost::program_options;

int main(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	std::string infile, outfile;
	float rho, thr;
	int maxit;
	bool transpose, text_out;

	desc.add_options()
		("help,h", "Display this message")
		("in,i",    bpo::value<std::string>(&infile)->required(), "Input file")
		("out,o",   bpo::value<std::string>(&outfile)->required(), "Output file")
		("rho,r",   bpo::value<float>(&rho)->required(), "The regularization parameter")
		("thr,t",   bpo::value<float>(&thr)->default_value(1.0e-6), "Convergence threshold")
		("maxit,m", bpo::value<int>(&maxit)->default_value(10000), "Maximum number of iterations")
		("transpose,t", bpo::bool_switch(&transpose)->default_value(false), "Should the input matrix be transposed.")
		("text,a", bpo::bool_switch(&text_out)->default_value(false), "Write the output as a text file.");

	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);

		if(!vm["help"].empty()) {
			desc.print(std::cout);
			return 0;
		}

		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return -1;
	}

	DenseMatrixReader reader;

	std::ifstream input(infile);

	if(!input) {
		std::cerr << "Could not open " << infile << " for reading." << std::endl;
		return -1;
	}

	unsigned int opt = DenseMatrixReader::defaultOptions();

	if(transpose) {
		opt |= DenseMatrixReader::TRANSPOSE;
	}

	std::cout << "Reading data ..." << std::endl;

	//TODO: See if we can make this casting business less wasteful
	DenseMatrix mat = reader.read(input, opt);

	auto mu = mat.matrix().rowwise().mean();
	mat.matrix().colwise() -= mu;

	std::cout << "Computing covariance matrix ..." << std::endl;
	// Compute the covariance matrix
	Eigen::MatrixXf cov = (mat.matrix() * mat.matrix().transpose() / (mat.cols() - 1)).cast<float>();

	// The glasso driver part. Mainly the creation of lots and lots of variables
	// in order to interface with Fortran...
	Eigen::MatrixXf www(mat.rows(), mat.rows());
	Eigen::MatrixXf wwwi(mat.rows(), mat.rows());

	// Wow! This is really a waste!
	// This code requires four times the memory of the input matrix...
	std::vector<float> rho_mem(mat.rows() * mat.rows(), 0.001);

	int nn = cov.rows();

	int ia = 0;
	int is = 0;
	int itr = 0;
	int ipen = 1;

	int nniter;
	float ddel;
	int jerr;

	std::cout << "Approximating precision matrix ..." << std::endl;

	glasso_(&nn, cov.data(), &rho_mem[0], &ia, &is, &itr, &ipen, &thr, &maxit,
	        www.data(), wwwi.data(), &nniter, &ddel, &jerr);

	std::ofstream out(outfile);
	if(!out) {
		std::cerr << "Could not open " << outfile << " for reading." << std::endl;
		return -1;
	}

	// We no longer need this stuff...
	cov.resize(0,0);
	www.resize(0,0);

	//TODO: Think about specializations for the matrix writer...
	//      OR make glasso double precision...
	DenseMatrix tmp(mat.rowNames(), mat.rowNames());
	tmp.matrix() = wwwi.cast<double>();

	DenseMatrixWriter writer;
	if(text_out) {
		writer.writeText(out, tmp);
	} else {
		writer.writeBinary(out, tmp);
	}

	return 0;
}
