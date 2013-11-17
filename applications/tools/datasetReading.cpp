/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2013 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 * Copyright (C) 2013 Tobias Frisch <tfrisch@bioinf.uni-sb.de>
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

#include "../../libraries/core/src/Dataset.h"
#include "../../libraries/core/src/DenseDatasetImpl.h"
#include "../../libraries/core/src/Parameter.h"
#include "../../libraries/core/src/Ontology.h"
#include "../../libraries/core/src/DenseSubset.h"

#include "../../libraries/core/src/CommandLineParser.h"

#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/any.hpp>
#include <unordered_map>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <cassert>

namespace fs = boost::filesystem;

using namespace GeneTrail;


int main(int argc, char** argv)
{
	std::cerr << "Bla" << std::endl;
	CommandLineParser p("Dataset reader. Import of a dataset and export into textfile");
	
	p.addOption("help,h", "Harhar");
	
	p.addTypedOption<std::string>("in,i","Input file");
	p.addTypedOption<std::string>("out,o", "Output directory");
	p.parse(argc,argv);
	
	if(p.checkParameter("help"))
	{
		p.printHelp();
		return 0;
	}
	
	if(!p.checkParameter("in"))
	{
		std::cerr << "ERROR: datasetReading: Missing input file" << std::endl;
		p.printHelp();
		return -1;
	}
	
	if(!p.checkParameter("out"))
	{
		std::cerr << "ERROR: datasetReading: Missing output file" << std::endl;
		p.printHelp();
		return -1;
	}
	
	std::string in;
	p.getParameter("in",in);
	std::string out;
	p.getParameter("out",out);
	
	fs::is_directory(out);
	
	std::fstream fstrm (in, std::fstream::in);
	
	if(!fstrm)
	{
		std::cerr << "ERROR: datasetReading: Can not open input file: " << in;
		return -1;
	}
	
	std::cout << "INFO: datasetReading: reading file: " << in ;
	std::cout.flush();
	DenseMatrixReader reader;
	DenseMatrix matrix = reader.read(fstrm);
	std::cout <<"....DONE" << std::endl;
	
	std::string organism = "human";
	std::string url = "hallowelt";
	Ontology::date_type date;
	std::string clinial = "some information";
	
	Ontology ontology(organism, url, date, clinial);
	
	Parameter parameter;
	
	Parameter::value_type a("Hallo");
	Parameter::value_type b(10);
	Parameter::value_type c("OHOH");
	
	Dataset::label_type first_row = matrix.rowName(0);
	Dataset::label_type second_row = matrix.rowName(1);
	Dataset::label_type third_row = matrix.rowName(2);
	
	parameter.insertRow(first_row,a);
	parameter.insertRow(second_row,b);
	parameter.insertRow(third_row,c);
	
	DenseDatasetImpl d(ontology,parameter,matrix);
	
	std::cout << "INFO: datasetReading: Writing Dataset: " << out;
	std::cout.flush();
	d.writeDataset(out);
	std::cout << "....DONE" << std::endl;
	
	std::cout << "INFO: datasetReading: Reading Dataset: " << out;
	std::cout.flush();
	DenseDatasetImpl readedSet(out);
	std::cout << "....DONE"  << std::endl;
	
	assert(d.cols() == readedSet.cols());
	assert(d.rows() == readedSet.rows());
	assert(d.parameter().col_parameter().size() == readedSet.parameter().col_parameter().size());
	assert(d.parameter().row_parameter().size() == readedSet.parameter().row_parameter().size());
	
	Dataset::labelVector_type rowSubset;
	Dataset::labelVector_type test = d.colNames();
	Dataset::labelVector_type colSubset(test.begin(),test.end());
	
	rowSubset.push_back(first_row);
	rowSubset.push_back(second_row);
	rowSubset.push_back(third_row);
	
	DenseSubset subset = readedSet.createSubset(rowSubset,colSubset);
	
	Ontology o1 = subset.ontology();
	
	std::string a1 = o1.clinicalInformation();
	std::string b1 = o1.organims();
	std::string c1 = o1.url();
	boost::gregorian::date d1 = o1.creationDate();
	
	Ontology o = subset.ontology();
	
	std::string a2 = o.clinicalInformation();
	std::string b2 = o.organims();
	std::string c2 = o.url();
	boost::gregorian::date d2 = o.creationDate();
	
 	assert(subset.rows() == 3);
	assert(subset.cols() == d.cols());
	
	Dataset::file_type subsetFile = out + "/Subset";
	std::cout << "INFO: datasetReading: Writing subset to " << subsetFile;
	std::cout.flush();
	subset.writeDataset(subsetFile);
	std::cout << "....DONE" << std::endl;
	
	
	std::cout << "INFO: datasetReading: Reading subset from "  << subsetFile;
	DenseSubset readedSubset;
	
	
	
	//TODO fix segmentation fault..!
	
	readedSubset.readDataset(subsetFile);
	std::cout << "....DONE" << std::endl;
	std::cerr << "readDataset" << std::endl;	
	Dataset::file_type subsetFile2 = out + "/SubsetSecond";
	readedSubset.writeDataset(subsetFile2);
	std::cerr << "wrote dataset" << std::endl;
}
