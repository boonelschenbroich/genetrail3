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

#ifndef GT2_DENSE_MATRIX_SUBSET_H
#define GT2_DENSE_MATRIX_SUBSET_H

#include "DenseMatrix.h"

#include "config.h"

namespace GeneTrail
{
	/**
	 * This class represents a subset of a DenseMatrix.
	 *
	 * The reason we derive from DenseMatrix is that we
	 * want to be able to use DenseMatrix and DenseMatrixSubset
	 * interchangably.
	 *
	 * Lets see if this works...
	 */
	class GT2_EXPORT DenseMatrixSubset : public Matrix
	{
		public:
			typedef std::vector<DenseMatrix::index_type> ISubset;
			typedef std::vector<std::string> SSubset;

			static DenseMatrixSubset createRowSubset(DenseMatrix* mat, ISubset rows);
			static DenseMatrixSubset createRowSubset(DenseMatrix* mat, const SSubset& rows);
			static DenseMatrixSubset createColSubset(DenseMatrix* mat, ISubset cols);
			static DenseMatrixSubset createColSubset(DenseMatrix* mat, const SSubset& cols);
			static DenseMatrixSubset createSubset   (DenseMatrix* mat, const SSubset& rows, const SSubset& cols);

			DenseMatrixSubset(DenseMatrix* mat, ISubset  rows, ISubset  cols);

			DenseMatrixSubset(const DenseMatrixSubset& subs) = default;
			DenseMatrixSubset(DenseMatrixSubset&& subs);

			DenseMatrixSubset& operator=(const DenseMatrixSubset& subs) = default;
			DenseMatrixSubset& operator=(DenseMatrixSubset&& subs);

			value_type& operator()(index_type i, index_type j);
			value_type  operator()(index_type i, index_type j) const;

			virtual const std::string& colName(index_type j) const override;
			virtual const std::string& rowName(index_type i) const override;

			virtual index_type colIndex(const std::string& col) const override;
			virtual index_type rowIndex(const std::string& row) const override;

			virtual index_type cols() const override;
			virtual index_type rows() const override;

			virtual bool hasCol(const std::string& name) const override;
			virtual bool hasRow(const std::string& name) const override;

			virtual void setColName(index_type j, const std::string& new_name) override;
			virtual void setColName(const std::string& old_name, const std::string& new_name) override;
			virtual void setColNames(const std::vector< std::string >& col_names) override;
			virtual void setRowName(index_type i, const std::string& new_name) override;
			virtual void setRowName(const std::string& old_name, const std::string& new_name) override;
			virtual void setRowNames(const std::vector< std::string >& row_names) override;

			virtual void removeCols(const std::vector< index_type >& indices) override;
			virtual void removeRows(const std::vector< index_type >& indices) override;
			virtual void shuffleCols(const std::vector< index_type >& perm) override;
			virtual void shuffleRows(const std::vector< index_type >& perm) override;
			virtual void transpose() override;

			virtual const std::vector< std::string >& colNames() const override;
			virtual const std::vector< std::string >& rowNames() const override;

		private:
			DenseMatrix* mat_;
			ISubset row_subset_;
			ISubset col_subset_;

			// TODO: Think of something smart to fix the hack below
			mutable SSubset row_names_cache_;
			mutable SSubset col_names_cache_;

			void remove_(const std::vector<Matrix::index_type>& indices, ISubset& subset);
	};
}

#endif // GT2_DENSE_MATRIX_SUBSET_H