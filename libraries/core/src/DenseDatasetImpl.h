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

#ifndef DENSE_DATASET_IMPL_H
#define DENSE_DATASET_IMPL_H


#include "DenseDataset.h"
#include "Ontology.h"
#include "Parameter.h"
#include "DenseSubset.h"

namespace GeneTrail
{
	class DenseDatasetImpl : public DenseDataset
	{
		public:
			
			DenseDatasetImpl(file_type& directory);
			
			DenseDatasetImpl(Ontology& ontology, Parameter& parameter, DenseMatrix& denseMatrix);
			
			DenseDatasetImpl(DenseDatasetImpl&& denseDataset);
			
			DenseDatasetImpl(const DenseDatasetImpl& denseDataset) = default;
			
			DenseDatasetImpl& operator=(DenseDatasetImpl&& dataset);
			
			DenseSubset& createSubset(labelVector_type& row_subset, labelVector_type& col_subset) override;
			
			
			colXpr_type col(const label_type i) override;
			value_type colIndex(const label_type key) override;
			const labelVector_type colNames() const override;
			index_type cols() override;
			bool hasCol(const label_type key) override;
			
			bool hasRow(const label_type key) override;
			rowXpr_type row(const label_type i) override;
			value_type rowIndex(const label_type key) override;
			const labelVector_type rowNames() const override;
			index_type rows() override;
			
			value_type value(const label_type i, const label_type j) override;
			value_type value(const index_type i, const index_type j) override;
			
			void readDataset(file_type& in) override;
			void writeDataset(file_type& out) override;
			
			bool disjunct(const Dataset& d) override;
			bool disjunctColumns(const Dataset& d) override;
			bool disjunctRows(const Dataset& d) override;
			
			subset_type& subsets() override;
			
			Ontology ontology() override;
			Parameter parameter() override;
			
			
			virtual DenseDatasetImpl* clone() const override;
			
			
		protected:
			void writeMetadata(file_type& out) override;
			void readMetadata(file_type& in) override;
			DenseMatrix readMatrix(file_type& in);
			
		private:
			
			Ontology ontology_;
			
			Parameter parameter_;
			
			DenseMatrix matrix_;
			
			subset_type subsets_;
	};
	
}


#endif // DENSE_DATASET_IMPL_H
