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


#include <fstream>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>

#include "DenseSubset.h"
#include "DenseDatasetImpl.h"

namespace fs = boost::filesystem;

namespace GeneTrail
{

	DenseSubset::DenseSubset(Dataset::labelVector_type& row_subset, Dataset::labelVector_type& col_subset, DenseDataset* dataset)
		: row_subset_(row_subset.begin(), row_subset.end()),
		  col_subset_(col_subset.begin(), col_subset.end()),
		  dataset_(dataset)
	{
	}

	DenseSubset::DenseSubset()
		: dataset_(nullptr)
	{
	}

	DenseSubset::DenseSubset(const DenseSubset& subs)
		: row_subset_(subs.row_subset_),
		  col_subset_(subs.col_subset_),
		  subsets_(subs.subsets_)
	{
		if(subs.dataset_) {
			dataset_ = subs.dataset_->clone();
		}
	}

	DenseSubset::~DenseSubset()
	{
		delete dataset_;
	}

	DenseSubset& DenseSubset::createSubset(Dataset::labelVector_type& row_subset, Dataset::labelVector_type& col_subset)
	{
		DenseSubset subset(row_subset, col_subset, this);
		
		subsets_.push_back(subset);

		return subsets_.back();
	}

	DenseSubset& DenseSubset::operator=(DenseSubset&& dataset)
	{
		row_subset_ = std::move(dataset.row_subset_);
		col_subset_ = std::move(dataset.col_subset_);
		dataset_ = std::move(dataset.dataset_);
		
		return *this;
	}

	DenseDataset::subset_type& DenseSubset::subsets()
	{
		return subsets_;
	}

	DenseDataset::colXpr_type DenseSubset::col(const Dataset::label_type i)
	{
		return dataset_->col(i);
	}

	Dataset::value_type DenseSubset::colIndex(const Dataset::label_type key)
	{
		return dataset_->colIndex(key);
	}

	const Dataset::labelVector_type DenseSubset::colNames() const
	{
		labelVector_type res(col_subset_.begin(), col_subset_.end());
		return res;
	}

	Dataset::index_type DenseSubset::cols()
	{
		return col_subset_.size();
	}

	bool DenseSubset::hasCol(const Dataset::label_type key)
	{
		return col_subset_.find(key) != col_subset_.end();
	}

	DenseDataset::rowXpr_type DenseSubset::row(const Dataset::label_type i)
	{
		return dataset_->row(i);
	}

	Dataset::value_type DenseSubset::rowIndex(const Dataset::label_type key)
	{
		return dataset_->rowIndex(key);
	}

	const Dataset::labelVector_type DenseSubset::rowNames() const
	{
		labelVector_type res(row_subset_.begin(),row_subset_.end());
		return res;
	}

	Dataset::index_type DenseSubset::rows()
	{
		return row_subset_.size();
	}

	bool DenseSubset::hasRow(const Dataset::label_type key)
	{
		return row_subset_.find(key) != row_subset_.end();
	}

	bool DenseSubset::disjoint(const Dataset& d)
	{
		return disjointColumns(d) && disjointRows(d);
	}

	bool DenseSubset::disjointColumns(const Dataset& d)
	{
		for(auto& col : d.colNames())
		{
			if(col_subset_.find(col) != col_subset_.end())
			{
				return false;
			}
		}
		return true;
	}

	bool DenseSubset::disjointRows(const Dataset& d)
	{
		for(auto& row : d.rowNames())
		{ 
			if(row_subset_.find(row) != row_subset_.end())
			{
				return false;
			}
		}
		return true;
	}

	void DenseSubset::readDataset(Dataset::file_type& in)
	{
		delete dataset_;
		dataset_ = new DenseDatasetImpl(in);

		dataset_->readDataset(in);
		
		file_type subsetPath = in + "/Subset.txt";
		
		std::fstream fstrm(subsetPath, std::fstream::in);
		
		if(!fstrm)
		{
			std::cerr << "Could not read path: " << subsetPath << std::endl;
			return;
		}
		
		std::string line;
		
		std::getline(fstrm,line);
		
		boost::trim(line);
		
		std::vector<std::string> splittedLine;
		
		boost::split(splittedLine, line, boost::is_any_of(" \t"), boost::token_compress_on);
		for(auto& a : splittedLine)
		{
			row_subset_.insert(a);
		}
		
		std::getline(fstrm,line);
		boost::trim(line);
		splittedLine.clear();
		boost::split(splittedLine,line, boost::is_any_of(" \t"), boost::token_compress_on);
		
		for(auto& a : splittedLine)
		{
			col_subset_.insert(a);
		}
		fstrm.close();
	}

	void DenseSubset::readMetadata(Dataset::file_type& in)
	{
		//does nothing cause dataset_ reades all metadata!
	}

	void DenseSubset::writeDataset(Dataset::file_type& out)
	{
		dataset_->writeDataset(out);
		
		file_type subsetPath = out + "/Subset.txt";
		
		std::fstream fstrm(subsetPath, std::fstream::out);
		
		if(!fstrm)
		{
			std::cerr << "ERROR: DenseSubset: Could not save subset" << std::endl;
			return;
		}
		
		for(auto& row : row_subset_)
		{
			fstrm << row << "\t";
		}
		
		fstrm << "\n";
		
		for(auto& col : col_subset_)
		{
			fstrm << col << "\t";
		}
		fstrm.close();
	}
	
	Dataset::value_type DenseSubset::value(const Dataset::label_type i, const Dataset::label_type j)
	{
		return dataset_->value(i,j);
	}

	Dataset::value_type DenseSubset::value(const Dataset::index_type i, const Dataset::index_type j)
	{
		return dataset_->value(i,j);
	}

	void DenseSubset::writeMetadata(Dataset::file_type& out)
	{
		//does nothing cause dataset_ safes all metadata!
	}
	
	Ontology DenseSubset::ontology()
	{
		return dataset_->ontology();
	}

	Parameter DenseSubset::parameter()
	{
		return dataset_->parameter();
	}

	DenseSubset* DenseSubset::clone() const
	{
		return new DenseSubset(*this);
	}

}

