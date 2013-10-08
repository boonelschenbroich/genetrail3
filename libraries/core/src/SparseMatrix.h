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

#ifndef GT2_SPARSE_MATRIX_H
#define GT2_SPARSE_MATRIX_H

#include "Matrix.h"

#include <Eigen/Sparse>

namespace GeneTrail
{
	class SparseMatrix : public Matrix
	{
		public:
			/// The Eigen class used for representing the internal matrix
			typedef Eigen::SparseMatrix<value_type> SMatrix;

			/**
			 * Default constructor
			 *
			 * Reserves the storage for a rows * cols matrix.
			 * The matrix will not be initialised.
			 */
			SparseMatrix(index_type rows, index_type cols);

			/**
			 * Row-name constructor
			 *
			 * This constructs a rows.size() x cols.size() matrix and sets the row and column names
			 * to rows and cols respectively
			 */
			SparseMatrix(const std::vector<std::string>&  rows, const std::vector<std::string>&  cols);

			/**
			 * @see SparseMatrix(const std::vector<std::string>&, const std::vector<std::string>&)
			 */
			SparseMatrix(std::vector<std::string>&& rows, const std::vector<std::string>&  cols);

			/**
			 * @see SparseMatrix(const std::vector<std::string>&, const std::vector<std::string>&)
			 */
			SparseMatrix(const std::vector<std::string>&  rows, std::vector<std::string>&& cols);

			/**
			 * @see SparseMatrix(const std::vector<std::string>&, const std::vector<std::string>&)
			 */
			SparseMatrix(std::vector<std::string>&& rows, std::vector<std::string>&& cols);

			/**
			 * Default copy constructor
			 */
			SparseMatrix(const SparseMatrix&) = default;
			/**
			 * Move Constructor
			 */
			SparseMatrix(SparseMatrix&& matrix);

			/**
			 * Returns a reference to the internal Eigen matrix
			 */
			SMatrix& matrix();

			/**
			 * Returns a reference to the internal Eigen matrix
			 */
			const SMatrix& matrix() const;

			///@}
			///\ingroup Operators
			///@{

			/**
			 * Returns a reference to the matrix coefficient at position (i,j)
			 *
			 * \warning This method is in O(1) for sparse matrices
			 */
			value_type& operator()(index_type i, index_type j);

			/**
			 * Returns the matrix coefficient at position (i,j)
			 *
			 * \warning This method is in O(1) for sparse matrices
			 */
			value_type  operator()(index_type i, index_type j) const;

			/**
			 * Move assignment operator
			 */
			SparseMatrix& operator=(SparseMatrix&& matrix);
			///@}

			///@}
			///\ingroup Matrix Operations
			///@{

			/**
			 * Reorder the rows of the matrix according to a given
			 * permutation.
			 *
			 * \param perm A vector of matrix indices. The new matrix
			 *             will consist of the rows at the indices
			 *             contained in perm in the order they were
			 *             specified.
			 */
			void shuffleRows(const std::vector<index_type>& perm);

			/**
			 * Reorder the columns of the matrix according to a given
			 * permutation.
			 *
			 * \see SparseMatrix::shuffleRows
			 */
			void shuffleCols(const std::vector<index_type>& perm);

			/**
			 * Remove the rows identified by the passed indices
			 *
			 * \param indices a vector of row indices
			 * \warning the passed indices must be sorted in ascending order
			 */
			void removeRows(const std::vector<index_type>& indices);

			/**
			 * Remove the columns identified by the passed indices
			 *
			 * \see SparseMatrix::removeRows
			 */
			void removeCols(const std::vector<index_type>& indices);

			/**
			 * Transpose the matrix.
			 *
			 * Row names and col names will also be switched.
			 *
			 * \note If you just need the transpose of the matrix in a computation use
			 *       .matrix().transpose()
			 */
			void transpose();

			///@}

		private:
			SMatrix m_;
	};
}

#endif // GT2_SPARSE_MATRIX_H