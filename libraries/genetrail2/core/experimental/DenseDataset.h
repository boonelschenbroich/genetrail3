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


#ifndef DENSE_DATASET_H
#define DENSE_DATASET_H

#include "Dataset.h"
#include "../DenseMatrix.h"

namespace GeneTrail
{
  class DenseSubset;
	class DenseDataset : public Dataset
	{
		
		public:
			virtual ~DenseDataset() {};
			
			typedef DenseMatrix::DMatrix matrix_type;
			typedef matrix_type::RowXpr rowXpr_type;
			typedef matrix_type::ColXpr colXpr_type;
			typedef std::vector<DenseSubset> subset_type;
			
			/**
			* 
			**/
			virtual value_type value(const index_type i,const index_type j) = 0;
			
			/**
			* 
			**/
			virtual value_type value(const label_type i, const label_type j) = 0;
			
			/**
			* 
			**/
			virtual colXpr_type col(const label_type i) = 0;
			
			/**
			* 
			**/
			virtual rowXpr_type row(const label_type i) = 0;
			
			/**
			* 
			**/
			virtual subset_type& subsets() = 0;
			
			
			virtual DenseDataset* clone() const override = 0;
			
		
		protected:

		private:
	};
}


#endif // DENSE_DATASET_H

