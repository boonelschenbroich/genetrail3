/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2013 Daniel Stöckel daniel@bioinf.uni-sb.de>
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

#include "Matrix.h"

#include <cassert>
#include <limits>

namespace GeneTrail
{

	Matrix::Matrix(index_type rows, index_type cols)
		: index_to_rowname_(rows),
		  index_to_colname_(cols)
	{
	}

	Matrix::Matrix(const std::vector<std::string>& rows, const std::vector<std::string>& cols)
		: index_to_rowname_(rows),
		  index_to_colname_(cols)
	{
		updateRowAndColNames_();
	}

	Matrix::Matrix(std::vector<std::string>&& rows, const std::vector<std::string>& cols)
		: index_to_rowname_(std::move(rows)),
		  index_to_colname_(cols)
	{
		updateRowAndColNames_();
	}

	Matrix::Matrix(const std::vector<std::string>& rows, std::vector<std::string>&& cols)
		: index_to_rowname_(rows),
		  index_to_colname_(std::move(cols))
	{
		updateRowAndColNames_();
	}

	Matrix::Matrix(std::vector<std::string>&& rows, std::vector<std::string>&& cols)
		: index_to_rowname_(std::move(rows)),
		  index_to_colname_(std::move(cols))
	{
		updateRowAndColNames_();
	}

	Matrix::Matrix(Matrix&& matrix)
		: index_to_rowname_(std::move(matrix.index_to_rowname_)),
		  rowname_to_index_(std::move(matrix.rowname_to_index_)),
		  index_to_colname_(std::move(matrix.index_to_colname_)),
		  colname_to_index_(std::move(matrix.colname_to_index_))
	{
	}

	Matrix& Matrix::operator=(Matrix&& matrix)
	{
		// Write this as an assert as moving to oneself is usually a grave client mistake...
		assert(this != &matrix);

		index_to_rowname_ = std::move(matrix.index_to_rowname_);
		rowname_to_index_ = std::move(matrix.rowname_to_index_);
		index_to_colname_ = std::move(matrix.index_to_colname_);
		colname_to_index_ = std::move(matrix.colname_to_index_);

		return *this;
	}

	void Matrix::updateRowAndColNames_()
	{
		index_type i = 0;
		for(const auto& s : index_to_rowname_) {
			rowname_to_index_.insert(std::make_pair(s, i++));
		}

		i = 0;
		for(const auto& s : index_to_colname_) {
			colname_to_index_.insert(std::make_pair(s, i++));
		}
	}


	Matrix::index_type Matrix::colIndex(const std::string& col) const
	{
		auto res = colname_to_index_.find(col);

		if(res != colname_to_index_.end()) {
			return res->second;
		}

		return std::numeric_limits<index_type>::max();
	}

	const std::string& Matrix::colName(index_type j) const
	{
		assert(j >= 0 && j < index_to_colname_.size());
		return index_to_colname_[j];
	}

	Matrix::index_type Matrix::cols() const
	{
		return index_to_colname_.size();
	}

	bool Matrix::hasCol(const std::string& name) const
	{
		return colname_to_index_.find(name) != colname_to_index_.end();
	}

	bool Matrix::hasRow(const std::string& name) const
	{
		return rowname_to_index_.find(name) != rowname_to_index_.end();
	}



	Matrix::index_type Matrix::rowIndex(const std::string& row) const
	{
		auto res = rowname_to_index_.find(row);

		if(res != rowname_to_index_.end()) {
			return res->second;
		}

		return std::numeric_limits<index_type>::max();
	}

	const std::string& Matrix::rowName(index_type i) const
	{
		assert(i >= 0 && i < index_to_rowname_.size());
		return index_to_rowname_[i];
	}

	Matrix::index_type Matrix::rows() const
	{
		return index_to_rowname_.size();
	}

	void Matrix::setColName(const std::string& old_name, const std::string& new_name)
	{
		setName_(old_name, new_name, colname_to_index_, index_to_colname_);
	}

	void Matrix::setColName(index_type j, const std::string& new_name)
	{
		setName_(j, new_name, colname_to_index_, index_to_colname_);
	}

	void Matrix::setColNames(const std::vector< std::string >& col_names)
	{
		assert(col_names.size() == cols());

		colname_to_index_.clear();

		index_to_colname_ = col_names;

		index_type i = 0;
		for(const auto& name : col_names) {
			colname_to_index_.insert(std::make_pair(name, i++));
		}
	}

	void Matrix::setName_(index_type j,
	                           const std::string& new_name,
	                           std::map<std::string, index_type>& name_to_index,
	                           std::vector<std::string>& index_to_name)
	{
		assert(j >= 0 && j < index_to_name.size());

		auto it = name_to_index.find(index_to_name[j]);
		if(it != name_to_index.end()) {
			name_to_index.erase(it);
		}

		it = name_to_index.find(new_name);

		// Steal the name from a possible previous assignment
		if(it != name_to_index.end()) {
			index_to_name[it->second] = "";
			it->second = j;
		} else {
			name_to_index.insert(std::make_pair(new_name, j));
		}

		index_to_name[j] = new_name;
	}

	void Matrix::setName_(const std::string& old_name,
	                           const std::string& new_name,
	                           std::map<std::string, index_type>& name_to_index,
	                           std::vector<std::string>& index_to_name)
	{
		auto res = name_to_index.find(old_name);

		if(res == name_to_index.end()) {
			return;
		}

		int index = res->second;

		name_to_index.erase(res);

		// Steal the index of existing rows
		res = name_to_index.find(new_name);

		if(res != name_to_index.end()) {
			index_to_name[res->second] = "";
			res->second = index;
		} else {
			name_to_index.insert(std::make_pair(new_name, index));
		}

		index_to_name[index] = new_name;
	}

	const std::vector<std::string>& Matrix::colNames() const
	{
		return index_to_colname_;
	}

	const std::vector<std::string>& Matrix::rowNames() const
	{
		return index_to_rowname_;
	}

	void Matrix::setRowName(const std::string& old_name, const std::string& new_name)
	{
		setName_(old_name, new_name, rowname_to_index_, index_to_rowname_);
	}

	void Matrix::setRowName(index_type i, const std::string& new_name)
	{
		setName_(i, new_name, rowname_to_index_, index_to_rowname_);
	}

	void Matrix::setRowNames(const std::vector< std::string >& row_names)
	{
		assert(row_names.size() == rows());

		rowname_to_index_.clear();

		index_to_rowname_ = row_names;

		index_type i = 0;
		for(const auto& name : row_names) {
			rowname_to_index_.insert(std::make_pair(name, i++));
		}
	}

	void Matrix::transpose()
	{
		std::swap(index_to_colname_, index_to_rowname_);
		std::swap(colname_to_index_, rowname_to_index_);
	}

}