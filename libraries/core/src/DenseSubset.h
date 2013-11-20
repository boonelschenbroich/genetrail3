/*
 * GeneTrail2 - An efficent library for interpreting genetic data
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


#ifndef DENSE_SUBSET_H
#define DENSE_SUBSET_H

#include "DenseDataset.h"
#include <set>

namespace GeneTrail
{
	class DenseSubset : public DenseDataset
	{
		public:
			
			typedef std::set<std::string> label_subset_type;
			
			DenseSubset();
			DenseSubset(const DenseSubset& );

			virtual ~DenseSubset();
			
			DenseSubset(labelVector_type& row_subset, labelVector_type& col_subset, DenseDataset* dataset);
			
			DenseSubset& createSubset(labelVector_type& row_subset, labelVector_type& col_subset) override;
			
			DenseSubset& operator=(DenseSubset&& dataset);
			
			
			colXpr_type col(const label_type i) override;
			value_type colIndex(const label_type key) override;
			index_type cols() override;
			bool hasCol(const label_type key) override;
			
			
			rowXpr_type row(const label_type i) override;
			value_type rowIndex(const label_type key) override;
			index_type rows() override;
			bool hasRow(const label_type key) override;
			
			subset_type& subsets() override;
			
			void readDataset(file_type& in) override;
			void writeDataset(file_type& out) override;
			
			bool disjoint(const Dataset& d) override;
			bool disjointColumns(const Dataset& d) override;
			bool disjointRows(const Dataset& d) override;
			
			const labelVector_type colNames() const override;
			const labelVector_type rowNames() const override;
			
			value_type value(const label_type i, const label_type j) override;
			value_type value(const index_type i, const index_type j) override;
		
			Ontology ontology() override;
			Parameter parameter() override;

			virtual DenseSubset* clone() const override;
			
			
		protected:
			void readMetadata(file_type& in) override;
			void writeMetadata(file_type& out) override;

		private:
			label_subset_type row_subset_;
			label_subset_type col_subset_;
			DenseDataset* dataset_;
			subset_type subsets_;
	};
}

#endif // DENSE_SUBSET_H
