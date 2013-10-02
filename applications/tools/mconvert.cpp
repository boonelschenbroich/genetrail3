/*
 * GeneTrail2 - An efficent library for interpreting genetic data
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

#include "../../libraries/core/src/DenseMatrix.h"
#include "../../libraries/core/src/DenseMatrixReader.h"
#include "../../libraries/core/src/RMAExpressMatrixReader.h"
#include "../../libraries/core/src/DenseMatrixWriter.h"

#include "../../libraries/core/src/CommandLineParser.h"

#include <fstream>

using namespace GeneTrail;

int main(int argc, char** argv)
{
	CommandLineParser p("Matrix Converter - Import into and conversion between GeneTrail matrix formats");

	p.addOption("help,h", "Display this message");
	p.addTypedOption<std::string>("in,i", "Input file");
	p.addTypedOption<std::string>("out,o", "Output file");
	p.addDefaultOption<std::string>("out-format,f", "binary", "Output format");
	p.addDefaultOption<bool>("transpose,t", false, "Transpose the input matrix");
	p.addDefaultOption<bool>("no-row-names,r", false, "The input has no row names (This only affects text matrices)");
	p.addDefaultOption<bool>("no-col-names,c", false, "The input has no column names (This only affects text matrices)");
	p.addDefaultOption<bool>("add-col-name,a", false, "The input has n+1 column names (This only affects text matrices)");
	p.addDefaultOption<bool>("rmaexpress,e", false, "Read the RMAExpress format for the input matrix");

	p.parse(argc, argv);

	if(p.checkParameter("help")) {
		p.printHelp();
		return 0;
	}

	if(!p.checkParameter("in")) {
		std::cerr << "Missing parameter in" << std::endl;
		p.printHelp();
		return -1;
	}

	if(!p.checkParameter("out")) {
		std::cerr << "Missing parameter in" << std::endl;
		p.printHelp();
		return -1;
	}

	std::string in;
	p.getParameter("in", in);
	std::ifstream istrm(in, std::ios::binary);

	if(!istrm) {
		std::cerr << "Could not open " << in << " for reading" << std::endl;
		return -1;
	}

	std::string out;
	p.getParameter("out", out);
	std::ofstream ostrm(out, std::ios::binary);

	if(!ostrm) {
		std::cerr << "Could not open " << out << " for writing" << std::endl;
		return -1;
	}

	DenseMatrixReader* reader;

	bool rmaexpress;
	p.getParameter("rmaexpress", rmaexpress);

	if(rmaexpress) {
		reader = new RMAExpressMatrixReader();
	} else {
		reader = new DenseMatrixReader();
	}

	DenseMatrixWriter writer;

	unsigned int opts = DenseMatrixReader::NO_OPTIONS;

	bool transpose;
	bool no_row_names;
	bool no_col_names;
	bool add_col_name;

	p.getParameter("transpose", transpose);
	p.getParameter("no-row-names", no_row_names);
	p.getParameter("no-col-names", no_col_names);
	p.getParameter("add-col-name", add_col_name);

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

	std::string out_format;
	p.getParameter("out-format", out_format);

	int errorcode = 0;
	if(out_format == "binary") {
		writer.writeBinary(ostrm, reader->read(istrm, opts));
	} else if(out_format == "ascii") {
		writer.writeText(ostrm, reader->read(istrm, opts));
	} else {
		std::cerr << "Invalid output format " << out_format << std::endl;
		errorcode = -1;
	}

	delete reader;

	return errorcode;
}