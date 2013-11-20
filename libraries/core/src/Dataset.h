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

#ifndef DATASET_H
#define DATASET_H



#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>

#include "Parameter.h"
#include "Ontology.h"

namespace GeneTrail{
	
	class Dataset{
		
		
		public:
			
			typedef std::string label_type;
			typedef double value_type;
			typedef unsigned int index_type;
			typedef std::string file_type;
			typedef std::vector<label_type> labelVector_type;
			
			
			/**
			 * This operator writes all information into the stream 
			 **/
			virtual void writeDataset (file_type& out) = 0; 
			
			/**
			 * This operator reads all information about the dataset out of the stream
			 **/
			virtual void readDataset (file_type& in) = 0;
			
			/**
			 * Generates a subset
			 **/
			virtual Dataset& createSubset(labelVector_type& row_subset, labelVector_type& col_subset) = 0;
			
			/**
			 * Return the column names
			 **/
			virtual const labelVector_type colNames() const = 0;
			
			/**
			 * Return the row names
			 **/
			virtual const labelVector_type rowNames() const = 0;
			
			/**
			 * 
			 **/
			virtual index_type cols() = 0;
			
			/**
			 * 
			 **/
			virtual index_type rows() = 0;
			
			/**
			* 
			**/
			virtual value_type rowIndex(const label_type key) = 0;
			
			/**
			* 
			**/
			virtual value_type colIndex(const label_type key) = 0;
			
			/**
			*
			**/
			virtual bool hasRow(const label_type key) = 0;
			
			/**
			* 
			**/
			virtual bool hasCol(const label_type key) = 0;
			
			/**
			 * 
			 **/
			virtual bool disjoint(const Dataset& d) = 0;
			
			/**
			 * 
			 **/
			virtual bool disjointRows(const Dataset& d) = 0;
			
			/**
			 * 
			 **/
			virtual bool disjointColumns(const Dataset& d) = 0;
			
			/**
			 * 
			 **/
			virtual Ontology ontology() = 0;
			
			/**
			 * 
			 **/
			virtual Parameter parameter() = 0;
			
			/**
			 * Creates a copy of this dataset
			 */
			virtual Dataset* clone() const = 0;
		protected:
			virtual void readMetadata(file_type& in) = 0;
			virtual void writeMetadata(file_type& out) = 0;
			
		private:
	};
	
}





#endif // DATASET_H
