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
		: AbstractMatrix(rows, cols),
		  m_(rows, cols)
	{
	}

	DenseMatrix::DenseMatrix(std::vector<std::string> rows, std::vector<std::string> cols)
		: AbstractMatrix(std::move(rows), std::move(cols)),
		  m_(this->rows(), this->cols())
	{
	}

	DenseMatrix::DenseMatrix(DenseMatrix&& matrix)
		: AbstractMatrix(std::move(matrix))
	{
		//TODO Replace this by the Eigen move constructor when it gets implemented
		m_.swap(matrix.m_);
	}

	DenseMatrix& DenseMatrix::operator=(DenseMatrix&& matrix)
	{
		//TODO Replace this by the Eigen move constructor when it gets implemented
		m_.swap(matrix.m_);
		AbstractMatrix::operator=(std::move(matrix));

		return *this;
	}

	DenseMatrix::DMatrix::ColXpr DenseMatrix::col(index_type j)
	{
		return m_.col(j);
	}

	DenseMatrix::DMatrix::ConstColXpr DenseMatrix::col(index_type j) const
	{
		return m_.col(j);
	}

	DenseMatrix::DMatrix& DenseMatrix::matrix()
	{
		return m_;
	}

	const DenseMatrix::DMatrix& DenseMatrix::matrix() const
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

	DenseMatrix::DMatrix::RowXpr DenseMatrix::row(index_type i)
	{
		return m_.row(i);
	}

	DenseMatrix::DMatrix::ConstRowXpr DenseMatrix::row(index_type i) const
	{
		return m_.row(i);
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
		AbstractMatrix::transpose();
		m_.transposeInPlace();
	}

}
