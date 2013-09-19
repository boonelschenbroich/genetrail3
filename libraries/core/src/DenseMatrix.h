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

#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H

#include <Eigen/Core>

#include <vector>
#include <map>

namespace GeneTrail
{
	/**
	 * A wrapper around MatrixXd that attaches row and column names.
	 *
	 * Additionally some convenience matrix manipulation and reshaping
	 * functions have been added. Serialize and deserialize this class
	 * using DenseMatrixWriter/Reader.
	 *
	 * \note DenseMatrix implements move constructor and assignment operators
	 * this means it is efficient to return temporary objects from a
	 * function.
	 */
	class DenseMatrix
	{
		public:
			/// The precision used in the matrix
			typedef double       value_type;

			/// The index type used in the matrix
			typedef unsigned int index_type;

			/// The Eigen class used for representing the internal matrix
			typedef Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> Matrix;

			/// The Eigen class used for representing rows and columns
			typedef Eigen::Matrix<value_type, Eigen::Dynamic, 1> Vector;

			/**
			 * Default constructor
			 *
			 * Reserves the storage for a rows * cols matrix.
			 * The matrix will not be initialised.
			 */
			DenseMatrix(index_type rows, index_type cols);

			/**
			 * Row-name constructor
			 *
			 * This constructs a rows.size() x cols.size() matrix and sets the row and column names
			 * to rows and cols respectively
			 */
			DenseMatrix(const std::vector<std::string>&  rows, const std::vector<std::string>&  cols);

			/**
			 * @see DenseMatrix(const std::vector<std::string>&, const std::vector<std::string>&)
			 */
			DenseMatrix(std::vector<std::string>&& rows, const std::vector<std::string>&  cols);

			/**
			 * @see DenseMatrix(const std::vector<std::string>&, const std::vector<std::string>&)
			 */
			DenseMatrix(const std::vector<std::string>&  rows, std::vector<std::string>&& cols);

			/**
			 * @see DenseMatrix(const std::vector<std::string>&, const std::vector<std::string>&)
			 */
			DenseMatrix(std::vector<std::string>&& rows, std::vector<std::string>&& cols);

			/**
			 * Default copy constructor
			 */
			DenseMatrix(const DenseMatrix&) = default;
			/**
			 * Move Constructor
			 */
			DenseMatrix(DenseMatrix&& matrix);

			///\ingroup Accessors
			///@{

			/**
			 * Sets the row names of the matrix to the values specified in row_names
			 */
			void setRowNames(const std::vector<std::string>& row_names);

			/**
			 * Sets the column names of the matrix to the values specified in col_names
			 */
			void setColNames(const std::vector<std::string>& col_names);

			/**
			 * Return the column names
			 */
			const std::vector<std::string>& colNames() const;

			/**
			 * Return the row names
			 */
			const std::vector<std::string>& rowNames() const;

			/**
			 * Rename row old_name to new_name
			 *
			 * This operation is a noop if old_name is not present in the matrix
			 *
			 * \warning If a row with name "new_name" is alread present
			 *          its row name will be set to the empty string
			 */
			void setRowName(const std::string& old_name, const std::string& new_name);

			/**
			 * Rename row i to new_name
			 *
			 * \warning If a row with name "new_name" is alread present
			 *          its row name will be set to the empty string
			 */
			void setRowName(index_type i, const std::string& new_name);

			/**
			 * Get the name of row i
			 */
			const std::string& rowName(index_type i) const;

			/**
			 * Rename column old_name to new_name
			 *
			 * This operation is a noop if old_name is not present in the matrix
			 *
			 * \warning If a column with name "new_name" is alread present
			 *          its colum name will be set to the empty string
			 */
			void setColName(const std::string& old_name, const std::string& new_name);

			/**
			 * Rename column j to new_name
			 *
			 * \warning If a column with name "new_name" is alread present
			 *          its colum name will be set to the empty string
			 */
			void setColName(index_type j, const std::string& new_name);

			/**
			 * Get the name of column j
			 */
			const std::string& colName(index_type j) const;

			/**
			 * Return the row index of the column identified by
			 * row.
			 *
			 * \warning the result of this method is undefined if
			 *          row is not a valid row name
			 */
			index_type rowIndex(const std::string& row) const;

			/**
			 * Return the column index of the column identified by
			 * column.
			 *
			 * \warning the result of this method is undefined if
			 *          column is not a valid column name
			 */
			index_type colIndex(const std::string& col) const;

			/**
			 * Return true if the matrix contains a row with row name
			 * "name". False otherwise.
			 */
			bool hasRow(const std::string& name) const;

			/**
			 * Return true if the matrix contains a column with column name
			 * "name". False otherwise.
			 */
			bool hasCol(const std::string& name) const;

			/**
			 * Return the i-th row.
			 */
			Matrix::RowXpr row(index_type i);

			/**
			 * Return the i-th row. Const version.
			 */
			Matrix::ConstRowXpr row(index_type i) const;

			/**
			 * Set the row identified by "name" to v
			 *
			 * This method is a noop if "name" is not present in the matrix
			 */
			void setRow(const std::string& name, const Vector& v);

			/**
			 * Set the i-th row to v
			 */
			void setRow(index_type i, const Vector& v);

			/**
			 * Return the j-th column
			 */
			Matrix::ColXpr col(index_type j);

			/**
			 * Return the j-th column
			 */
			Matrix::ConstColXpr col(index_type j) const;

			/**
			 * Set the column identified by "name" to v
			 *
			 * This method is a noop if "name" is not present in the matrix
			 */
			void setCol(const std::string& name, const Vector& v);

			/**
			 * Set the i-th column to v
			 */
			void setCol(index_type j, const Vector& v);

			/**
			 * The number of columns in the matrix
			 */
			index_type cols() const;

			/**
			 * The number of rows in the matrix
			 */
			index_type rows() const;

			/**
			 * Returns a reference to the internal Eigen matrix
			 */
			Matrix& matrix();

			/**
			 * Returns a reference to the internal Eigen matrix
			 */
			const Matrix& matrix() const;

			///@}
			///\ingroup Operators
			///@{

			/**
			 * Returns a reference to the matrix coefficient at position (i,j)
			 */
			value_type& operator()(index_type i, index_type j);

			/**
			 * Returns the matrix coefficient at position (i,j)
			 */
			value_type  operator()(index_type i, index_type j) const;

			/**
			 * Move assignment operator
			 */
			DenseMatrix& operator=(DenseMatrix&& matrix);
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
			 * \see DenseMatrix::shuffleRows
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
			 * \see DenseMatrix::removeRows
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
			// Actual matrix payload
			Matrix m_;

			// Containers for mapping row names
			std::vector<std::string> index_to_rowname_;
			std::map<std::string, index_type> rowname_to_index_;

			// Containers for mapping column names
			std::vector<std::string> index_to_colname_;
			std::map<std::string, index_type> colname_to_index_;

			// Helper for removing rows or columns
			void remove_(const std::vector<index_type>& indices,
			                   std::map<std::string, index_type>& name_to_index,
			                   std::vector<std::string>& index_to_name,
			                   std::function<void(index_type, index_type)> copy);

			// Helper for renaming rows or columns
			void setName_(index_type j,
			              const std::string& new_name,
			              std::map<std::string, index_type>& name_to_index,
			              std::vector<std::string>& index_to_name);

			void setName_(const std::string& old_name,
			              const std::string& new_name,
			              std::map<std::string, index_type>& name_to_index,
			              std::vector<std::string>& index_to_name);

			void updateRowAndColNames_();
	};
}

#endif // DENSEMATRIX_H
