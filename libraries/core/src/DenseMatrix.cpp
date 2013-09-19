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

#include "DenseMatrix.h"

#include <limits>
#include <cassert>

namespace GeneTrail
{

	DenseMatrix::DenseMatrix(index_type rows, index_type cols)
		: m_(rows, cols),
		  index_to_rowname_(rows),
		  index_to_colname_(cols)
	{
	}

	DenseMatrix::DenseMatrix(const std::vector<std::string>& rows, const std::vector<std::string>& cols)
		: m_(rows.size(), cols.size()),
		  index_to_rowname_(rows),
		  index_to_colname_(cols)
	{
		updateRowAndColNames_();
	}

	DenseMatrix::DenseMatrix(std::vector<std::string>&& rows, const std::vector<std::string>& cols)
		: m_(rows.size(), cols.size()),
		  index_to_rowname_(std::move(rows)),
		  index_to_colname_(cols)
	{
		updateRowAndColNames_();
	}

	DenseMatrix::DenseMatrix(const std::vector<std::string>& rows, std::vector<std::string>&& cols)
		: m_(rows.size(), cols.size()),
		  index_to_rowname_(rows),
		  index_to_colname_(std::move(cols))
	{
		updateRowAndColNames_();
	}

	DenseMatrix::DenseMatrix(std::vector<std::string>&& rows, std::vector<std::string>&& cols)
		: m_(rows.size(), cols.size()),
		  index_to_rowname_(std::move(rows)),
		  index_to_colname_(std::move(cols))
	{
		updateRowAndColNames_();
	}

	void DenseMatrix::updateRowAndColNames_()
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

	DenseMatrix::DenseMatrix(DenseMatrix&& matrix)
		: index_to_rowname_(std::move(matrix.index_to_rowname_)),
		  rowname_to_index_(std::move(matrix.rowname_to_index_)),
		  index_to_colname_(std::move(matrix.index_to_colname_)),
		  colname_to_index_(std::move(matrix.colname_to_index_))
	{
		//TODO Replace this by the Eigen move constructor when it gets implemented
		m_.swap(matrix.m_);
	}

	DenseMatrix& DenseMatrix::operator=(DenseMatrix&& matrix)
	{
		// Write this as an assert as moving to oneself is usually a grave client mistake...
		assert(this != &matrix);

		index_to_rowname_ = std::move(matrix.index_to_rowname_);
		rowname_to_index_ = std::move(matrix.rowname_to_index_);
		index_to_colname_ = std::move(matrix.index_to_colname_);
		colname_to_index_ = std::move(matrix.colname_to_index_);

		//TODO Replace this by the Eigen move constructor when it gets implemented
		m_.swap(matrix.m_);

		return *this;
	}

	DenseMatrix::Matrix::ColXpr DenseMatrix::col(index_type j)
	{
		return m_.col(j);
	}

	DenseMatrix::Matrix::ConstColXpr DenseMatrix::col(index_type j) const
	{
		return m_.col(j);
	}

	DenseMatrix::index_type DenseMatrix::colIndex(const std::string& col) const
	{
		auto res = colname_to_index_.find(col);

		if(res != colname_to_index_.end()) {
			return res->second;
		}

		return std::numeric_limits<index_type>::max();
	}

	const std::string& DenseMatrix::colName(index_type j) const
	{
		assert(j >= 0 && j < index_to_colname_.size());
		return index_to_colname_[j];
	}

	DenseMatrix::index_type DenseMatrix::cols() const
	{
		return m_.cols();
	}

	bool DenseMatrix::hasCol(const std::string& name) const
	{
		return colname_to_index_.find(name) != colname_to_index_.end();
	}

	bool DenseMatrix::hasRow(const std::string& name) const
	{
		return rowname_to_index_.find(name) != rowname_to_index_.end();
	}

	DenseMatrix::Matrix& DenseMatrix::matrix()
	{
		return m_;
	}

	const DenseMatrix::Matrix& DenseMatrix::matrix() const
	{
		return m_;
	}

	DenseMatrix::value_type& DenseMatrix::operator()(index_type i, index_type j)
	{
		return m_(i,j);
	}

	DenseMatrix::value_type DenseMatrix::operator()(index_type i, index_type j) const
	{
		return m_(i,j);
	}

	void DenseMatrix::remove_(const std::vector<index_type>& indices,
	                                std::map<std::string, index_type>& name_to_index,
	                                std::vector<std::string>& index_to_name,
	                                std::function<void(index_type, index_type)> copy)
	{
		size_t next_idx = 1;
		size_t write_idx = indices[0];
		name_to_index.erase(index_to_name[write_idx]);

		for(size_t read_idx = indices[0] + 1; read_idx < index_to_name.size(); ++read_idx)
		{
			// If we reached another column that should be deleted, just advance the read pointer
			// the next index
			if((size_t)next_idx < indices.size() && read_idx == indices[next_idx])
			{
				assert(indices[next_idx - 1] < indices[next_idx]);
				// Update the name <-> index mapping
				name_to_index.erase(index_to_name[read_idx]);
				++next_idx;
			}
			else
			{
				// Update the name <-> index mapping
				auto it = name_to_index.find(index_to_name[write_idx]);

				if(it != name_to_index.end() && it->second == write_idx) {
					name_to_index.erase(it);
				}

				const std::string name = index_to_name[read_idx];
				index_to_name[write_idx] = name;
				name_to_index[name] = write_idx;

				copy(write_idx, read_idx);
				++write_idx;
			}
		}

		// Check consistency of the name maps
		assert(name_to_index.size() == index_to_name.size() - indices.size());
	}

	void DenseMatrix::removeCols(const std::vector< index_type >& indices)
	{
		//TODO consolidate with removeRows
		if(indices.empty()) {
			return;
		}

		// Call the generalized remove function
		remove_(indices, colname_to_index_, index_to_colname_, [this](index_type i, index_type j) {
			m_.col(i) = m_.col(j);
		});

		// Free the memory of the unneeded columns
		m_.conservativeResize(Eigen::NoChange, m_.cols() - indices.size());
		index_to_colname_.resize(m_.cols());
	}


	void DenseMatrix::removeRows(const std::vector< index_type >& indices)
	{
		//TODO consolidate with removeRows
		if(indices.empty()) {
			return;
		}

		// Call the generalized remove function
		remove_(indices, rowname_to_index_, index_to_rowname_, [this](index_type i, index_type j) {
			m_.row(i) = m_.row(j);
		});

		// Free the memory of the unneeded columns
		m_.conservativeResize(m_.rows() - indices.size(), Eigen::NoChange);
		index_to_rowname_.resize(m_.rows());
	}

	DenseMatrix::Matrix::RowXpr DenseMatrix::row(index_type i)
	{
		return m_.row(i);
	}

	DenseMatrix::Matrix::ConstRowXpr DenseMatrix::row(index_type i) const
	{
		return m_.row(i);
	}

	DenseMatrix::index_type DenseMatrix::rowIndex(const std::string& row) const
	{
		auto res = rowname_to_index_.find(row);

		if(res != rowname_to_index_.end()) {
			return res->second;
		}

		return std::numeric_limits<index_type>::max();
	}

	const std::string& DenseMatrix::rowName(index_type i) const
	{
		assert(i >= 0 && i < index_to_rowname_.size());
		return index_to_rowname_[i];
	}

	DenseMatrix::index_type DenseMatrix::rows() const
	{
		return m_.rows();
	}

	void DenseMatrix::setCol(const std::string& name, const DenseMatrix::Vector& v)
	{
		auto res = colname_to_index_.find(name);

		if(res != colname_to_index_.end()) {
			m_.col(res->second) = v;
		}
	}

	void DenseMatrix::setCol(index_type j, const DenseMatrix::Vector& v)
	{
		m_.col(j) = v;
	}

	void DenseMatrix::setColName(const std::string& old_name, const std::string& new_name)
	{
		setName_(old_name, new_name, colname_to_index_, index_to_colname_);
	}

	void DenseMatrix::setColName(index_type j, const std::string& new_name)
	{
		setName_(j, new_name, colname_to_index_, index_to_colname_);
	}

	void DenseMatrix::setColNames(const std::vector< std::string >& col_names)
	{
		assert(col_names.size() == (size_t)m_.cols());

		colname_to_index_.clear();

		index_to_colname_ = col_names;

		index_type i = 0;
		for(const auto& name : col_names) {
			colname_to_index_.insert(std::make_pair(name, i++));
		}
	}

	void DenseMatrix::setName_(index_type j,
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

	void DenseMatrix::setName_(const std::string& old_name,
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

	const std::vector<std::string>& DenseMatrix::colNames() const
	{
		return index_to_colname_;
	}

	const std::vector<std::string>& DenseMatrix::rowNames() const
	{
		return index_to_rowname_;
	}

	void DenseMatrix::setRow(const std::string& name, const DenseMatrix::Vector& v)
	{
		auto res = rowname_to_index_.find(name);

		if(res != rowname_to_index_.end()) {
			m_.row(res->second) = v.transpose();
		}
	}

	void DenseMatrix::setRow(index_type i, const DenseMatrix::Vector& v)
	{
		m_.row(i) = v.transpose();
	}

	void DenseMatrix::setRowName(const std::string& old_name, const std::string& new_name)
	{
		setName_(old_name, new_name, rowname_to_index_, index_to_rowname_);
	}

	void DenseMatrix::setRowName(index_type i, const std::string& new_name)
	{
		setName_(i, new_name, rowname_to_index_, index_to_rowname_);
	}

	void DenseMatrix::setRowNames(const std::vector< std::string >& row_names)
	{
		assert(row_names.size() == (size_t)m_.rows());

		rowname_to_index_.clear();

		index_to_rowname_ = row_names;

		index_type i = 0;
		for(const auto& name : row_names) {
			rowname_to_index_.insert(std::make_pair(name, i++));
		}
	}

	void DenseMatrix::shuffleCols(const std::vector< index_type >& perm)
	{
		//TODO implement
		assert(false);
	}

	void DenseMatrix::shuffleRows(const std::vector< index_type >& perm)
	{
		//TODO implement
		assert(false);
	}

	void DenseMatrix::transpose()
	{
		m_.transposeInPlace();
		std::swap(index_to_colname_, index_to_rowname_);
		std::swap(colname_to_index_, rowname_to_index_);
	}

}
