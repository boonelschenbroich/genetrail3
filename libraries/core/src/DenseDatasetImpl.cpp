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


#include <istream>
#include <fstream>
#include <cassert>
#include <set>

#include "DenseDatasetImpl.h"
#include "DenseMatrixReader.h"
#include "DenseMatrixWriter.h"

// #include <Eigen/Core>
namespace fs = boost::filesystem;
namespace GeneTrail
{
	
	DenseDatasetImpl::DenseDatasetImpl(Dataset::file_type& directory)
	 : matrix_(readMatrix(directory))
	{
		readMetadata(directory);
	}



	DenseDatasetImpl::DenseDatasetImpl(Ontology& ontology, Parameter& parameter, DenseMatrix& denseMatrix)
		: ontology_(ontology),
			parameter_(parameter),
			matrix_(denseMatrix)
	{
	}

	DenseDatasetImpl::DenseDatasetImpl(DenseDatasetImpl&& denseDataset)
		: ontology_(std::move(denseDataset.ontology_)),
			parameter_(std::move(denseDataset.parameter_)),
			matrix_(std::move(denseDataset.matrix_)),
			subsets_(std::move(denseDataset.subsets_))
	{
	}

	DenseDatasetImpl& DenseDatasetImpl::operator=(DenseDatasetImpl&& dataset)
	{
		assert(this != &dataset);
		ontology_ = std::move(dataset.ontology_);
		parameter_ = std::move(dataset.parameter_);
		matrix_ = std::move(dataset.matrix_);
		subsets_ = std::move(dataset.subsets_);
		
		return *this;
	}

	DenseSubset& DenseDatasetImpl::createSubset(Dataset::labelVector_type& row_subset, Dataset::labelVector_type& col_subset)
	{
		DenseSubset subset(row_subset, col_subset, this);
		
		subsets_.push_back(subset);

		return subsets_.back();
	}

	bool DenseDatasetImpl::disjunct(const Dataset& d)
	{
		return disjunctColumns(d) && disjunctRows(d);
	}

	bool DenseDatasetImpl::disjunctColumns(const Dataset& d)
	{
		Dataset::labelVector_type labels = matrix_.colNames();
		std::set<label_type> ownLabels(labels.begin(), labels.end());
		labels.clear();
		for(auto& col : d.colNames())
		{
			if(ownLabels.find(col) != ownLabels.end())
			{
				return false;
			}
		}
		return true;
	}

	bool DenseDatasetImpl::disjunctRows(const Dataset& d)
	{
		Dataset::labelVector_type labels = matrix_.rowNames();
		std::set<label_type> ownLabels(labels.begin(), labels.end());
		labels.clear();
		for(auto& row : d.rowNames())
		{
			if(ownLabels.find(row) != ownLabels.end())
			{
				return false;
			}
		}
		
		return true;
	}

	DenseDataset::colXpr_type DenseDatasetImpl::col(const Dataset::label_type i)
	{
		return matrix_.col(matrix_.colIndex(i));
	}

	Dataset::value_type DenseDatasetImpl::colIndex(const Dataset::label_type key)
	{
		return matrix_.colIndex(key);
	}

	const Dataset::labelVector_type DenseDatasetImpl::colNames() const
	{
		return matrix_.colNames();
	}

	Dataset::index_type DenseDatasetImpl::cols()
	{
		return matrix_.cols();
	}

	bool DenseDatasetImpl::hasCol(const Dataset::label_type key)
	{
		return matrix_.hasCol(key);
	}

	DenseDataset::rowXpr_type DenseDatasetImpl::row(const Dataset::label_type i)
	{
		return matrix_.row(matrix_.rowIndex(i));
	}

	Dataset::value_type DenseDatasetImpl::rowIndex(const Dataset::label_type key)
	{
		return matrix_.rowIndex(key);
	}

	const Dataset::labelVector_type DenseDatasetImpl::rowNames() const
	{
		return matrix_.rowNames();
	}

	Dataset::index_type DenseDatasetImpl::rows()
	{
		return matrix_.rows();
	}

	bool DenseDatasetImpl::hasRow(const Dataset::label_type key)
	{
		return matrix_.hasRow(key);
	}

	Dataset::value_type DenseDatasetImpl::value(const Dataset::label_type i, const Dataset::label_type j)
	{
		return matrix_(rowIndex(i), colIndex(j));
	}

	Dataset::value_type DenseDatasetImpl::value(const Dataset::index_type i, const Dataset::index_type j)
	{
		return matrix_(i,j);
	}

	void DenseDatasetImpl::readDataset(Dataset::file_type& in)
	{
		if(!fs::exists(in))
		{
			std::cerr << "ERROR: DenseDatasetImpl: Given path does not exists!" << std::cerr;
		}
		
		if(!fs::is_directory(in))
		{
			std::cerr << "ERROR: DenseDatasetImpl: Given path should be directory!" << std::cerr;
			return;
		}
		readMetadata(in);
		std::string matrixPath = in + "/Matrix.txt";
		
		std::fstream fStrm(matrixPath, std::fstream::in);
		
		if(!fStrm)
		{
			std::cerr << "ERROR: DensedatasetImpl: Could not find Matrix.txt" << std::endl;
			return;
		}
		
		DenseMatrixReader reader;
		matrix_ = reader.read(fStrm);
		fStrm.close();
	}

	void DenseDatasetImpl::readMetadata(Dataset::file_type& in)
	{
		file_type parameterPath = in + "/Parameter.txt";
		file_type ontologyPath = in + "/Ontology.txt";
		
		
		std::fstream parameterStrm(parameterPath, std::fstream::in);
		std::fstream ontologyStrm(ontologyPath, std::fstream::in);
		
		if(!parameterStrm || !ontologyStrm)
		{
			std::cerr << "ERROR: DenseDatasetImpl: Could not read metadata" << std::endl;
			return;
		}
		
		parameter_.readParameter(parameterStrm);
		parameterStrm.close();
		ontology_.readOntology(ontologyStrm);
		ontologyStrm.close();
	}

	void DenseDatasetImpl::writeDataset(Dataset::file_type& out)
	{
		if(!fs::exists(out))
		{
			fs::create_directory(out);
		}
		
		if(!fs::exists(out) || !fs::is_directory(out))
		{
			std::cerr << "ERROR: DenseDatasetImpl: Could not create directory: " << out << std::endl;
		}
		
		writeMetadata(out);
		
		std::fstream fstrm(out + "/Matrix.txt", std::fstream::out);
		
		if(!fstrm)
		{
			std::cerr << "ERROR: DenseDatasetImpl: Could not write matrix!" << std::endl;
			return;
		}
		
		DenseMatrixWriter writer;
		writer.writeText(fstrm,matrix_);
		fstrm.close();
	}

	void DenseDatasetImpl::writeMetadata(Dataset::file_type& out)
	{
		file_type parameterPath = out + "/Parameter.txt";
		file_type ontologyPath = out + "/Ontology.txt";
		
		std::fstream parameterStrm(parameterPath, std::fstream::out);
		std::fstream ontologyStrm(ontologyPath, std::fstream::out);
		
		if(!parameterStrm || !ontologyStrm)
		{
			std::cerr << "ERROR: DenseDatasetImple: Could not write metadata" << std::endl;
		}
		
		parameter_.writeParameter(parameterStrm);
		parameterStrm.close();
		ontology_.writeOntology(ontologyStrm);
		ontologyStrm.close();
	}
	
DenseMatrix DenseDatasetImpl::readMatrix(Dataset::file_type& in)
{
	DenseMatrixReader reader;
	std::fstream fstrm(in + "/Matrix.txt", std::fstream::in);
	if(!fstrm)
	{
		std::cerr << "ERROR: DenseDatasetImpl: Could not read Matrix from file: " << in << std::endl;
	}
	
	return reader.read(fstrm);
}


	DenseDataset::subset_type& DenseDatasetImpl::subsets()
	{
		return subsets_;
	}

Ontology DenseDatasetImpl::ontology()
{
	return ontology_;
}

Parameter DenseDatasetImpl::parameter()
{
	return parameter_;
}



}
