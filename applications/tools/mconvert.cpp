/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2013 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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

#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/RMAExpressMatrixReader.h>
#include <genetrail2/core/DenseMatrixWriter.h>
#include <genetrail2/core/DenseRowSubset.h>
#include <genetrail2/core/DenseColumnSubset.h>

#include <fstream>
#include <iostream>

#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>

using namespace GeneTrail;
namespace bpo = boost::program_options;

bool readSubset(const std::string& file, std::vector<std::string>& out)
{
	std::ifstream in(file);

	if(!in) {
		std::cerr << "Could not open file " << file << " for reading!" << std::endl;
		return false;
	}

	std::string line;
	while(std::getline(in, line)) {
		boost::trim(line);

		if(line == "") {
			continue;
		}

		out.push_back(line);
	}

	return true;
}

template<class Mat>
int writeMatrix(std::ostream& ostrm, std::string out_format, const Mat& inmat)
{
	DenseMatrixWriter writer;

	if(out_format == "binary") {
		writer.writeBinary(ostrm, inmat);
	} else if(out_format == "ascii") {
		writer.writeText(ostrm, inmat);
	} else {
		std::cerr << "Invalid output format " << out_format << std::endl;
		return -1;
	}

	return 0;
}

// RAII class that handles the destruction of the output stream.
// This makes the main code much cleaner.
// TODO: Maybe move this to an auxillary header.
class Outstream {
	public:
		explicit Outstream(const std::string& strm)
			: deconstruct_(false), strm_(nullptr)
		{
			if(strm == "stdout") {
				strm_ = &std::cout;
			} else if(strm == "stderr") {
				strm_ = &std::cerr;
			} else {
				deconstruct_ = true;
			    strm_ = new std::ofstream(strm, std::ios::binary);
		    }
		}

	    ~Outstream()
	    {
		    if(deconstruct_)
			    delete strm_;
	    }

		std::ostream& operator()() {
			return *strm_;
		}

	private:
		bool deconstruct_;
		std::ostream* strm_;
};

int main(int argc, char** argv)
{
	bpo::variables_map vm;
	bpo::options_description desc;

	std::string infile, outfile, out_format, col_subset, row_subset;
	bool transpose, no_row_names, no_col_names, add_col_name, rmaexpress;

	desc.add_options()
		("help,h", "Display this message")
		("in,i",           bpo::value<std::string>(&infile)->required(), "Input file")
		("out,o",          bpo::value<std::string>(&outfile)->required(), "Output file. Use stdout or stderr to write to console.")
		("out-format,f",   bpo::value<std::string>(&out_format)->default_value("binary"), "Output format")
		("transpose,t",    bpo::value<bool>(&transpose)->default_value(false)->zero_tokens(), "Transpose the input matrix")
		("no-row-names,r", bpo::value<bool>(&no_row_names)->default_value(false)->zero_tokens(), "The input has no row names (This only affects text matrices)")
		("no-col-names,c", bpo::value<bool>(&no_col_names)->default_value(false)->zero_tokens(), "The input has no column names (This only affects text matrices)")
		("add-col-name,a", bpo::value<bool>(&add_col_name)->default_value(false)->zero_tokens(), "The input has n+1 column names (This only affects text matrices)")
		("col-subset,s",   bpo::value<std::string>(&col_subset), "An optional file containing column names that should be selected from the matrix")
		("row-subset,d",   bpo::value<std::string>(&row_subset), "An optional file containing row names that should be selected from the matrix")
		("rmaexpress,e",   bpo::value<bool>(&rmaexpress)->default_value(false)->zero_tokens(), "Read the RMAExpress format for the input matrix");

	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return -1;
	}

	std::ifstream istrm(infile, std::ios::binary);

	if(!istrm) {
		std::cerr << "Could not open " << infile << " for reading" << std::endl;
		return -1;
	}

	// Instantiate a RAII class that manages the output stream
	Outstream ostrm(outfile);

	if(!ostrm()) {
		std::cerr << "Could not open " << outfile << " for writing" << std::endl;
		return -1;
	}

	unsigned int opts = DenseMatrixReader::NO_OPTIONS;

	if(transpose) {
		opts |= DenseMatrixReader::TRANSPOSE;
	}

	if(!no_row_names) {
		opts |= DenseMatrixReader::READ_ROW_NAMES;
	}

	if(!no_col_names) {
		opts |= DenseMatrixReader::READ_COL_NAMES;
	}

	if(add_col_name) {
		opts |= DenseMatrixReader::ADDITIONAL_COL_NAME;
	}

	if(!vm["col-subset"].empty() && !vm["row-subset"].empty()) {
		std::cerr << "Currently it is not supported to subset both: columns and rows\n";
		return -2;
	}

	auto reader = rmaexpress ? std::make_shared<RMAExpressMatrixReader>()
	                         : std::make_shared<DenseMatrixReader>();

	DenseMatrix inmat = reader->read(istrm, opts);

	if(!vm["col-subset"].empty()) {
		std::vector<std::string> cs;
		if(!readSubset(col_subset, cs)) {
			return -1;
		}

		return writeMatrix(ostrm(), out_format, DenseColumnSubset(&inmat, cs));
	}

	if(!vm["row-subset"].empty()) {
		std::vector<std::string> rs;

		if(!readSubset(row_subset, rs)) {
			return -1;
		}

		return writeMatrix(ostrm(), out_format, DenseRowSubset(&inmat, rs));
	}

	return writeMatrix(ostrm(), out_format, inmat);
}
