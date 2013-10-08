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

#include "SparseMatrix.h"

namespace GeneTrail
{
	SparseMatrix::SparseMatrix(Matrix::index_type rows, Matrix::index_type cols)
		: Matrix(rows, cols),
		  m_(this->rows(), this->cols())
	{
	}

	SparseMatrix::SparseMatrix(const std::vector< std::string >& rows, const std::vector< std::string >& cols)
		: Matrix(rows, cols),
		  m_(this->rows(), this->cols())
	{
	}

	SparseMatrix::SparseMatrix(std::vector< std::string >&& rows, const std::vector< std::string >& cols)
		: Matrix(cols, cols),
		  m_(this->rows(), this->cols())
	{
	}

	SparseMatrix::SparseMatrix(const std::vector< std::string >& rows, std::vector< std::string >&& cols)
		: Matrix(rows, rows),
		  m_(this->rows(), this->cols())
	{
	}

	SparseMatrix::SparseMatrix(std::vector< std::string >&& rows, std::vector< std::string >&& cols)
		: Matrix(rows, cols),
		  m_(this->rows(), this->cols())
	{
	}

	SparseMatrix::SparseMatrix(SparseMatrix&& matrix)
		: Matrix(std::move(matrix))
	{
		m_.swap(matrix.m_);
	}

	SparseMatrix& SparseMatrix::operator=(SparseMatrix&& matrix)
	{
		m_.swap(matrix.m_);

		Matrix::operator=(std::move(matrix));

		return *this;
	}

	SparseMatrix::SMatrix& SparseMatrix::matrix()
	{
		return m_;
	}

	const SparseMatrix::SMatrix& SparseMatrix::matrix() const
	{
		return m_;
	}

	Matrix::value_type& SparseMatrix::operator()(Matrix::index_type i, Matrix::index_type j)
	{
		return m_.coeffRef(i,j);
	}

	Matrix::value_type SparseMatrix::operator()(Matrix::index_type i, Matrix::index_type j) const
	{
		return m_.coeff(i,j);
	}

	void SparseMatrix::removeCols(const std::vector< Matrix::index_type >& indices)
	{
		assert(false);
	}

	void SparseMatrix::removeRows(const std::vector< Matrix::index_type >& indices)
	{
		assert(false);
	}

	void SparseMatrix::shuffleCols(const std::vector< Matrix::index_type >& perm)
	{
		assert(false);
	}

	void SparseMatrix::shuffleRows(const std::vector< Matrix::index_type >& perm)
	{
		assert(false);
	}

	void SparseMatrix::transpose()
	{
		Matrix::transpose();
	}
}
