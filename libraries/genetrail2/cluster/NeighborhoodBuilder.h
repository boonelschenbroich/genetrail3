/*
 * GeneTrail2 - An efficient library for interpreting genetic data
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
 */

#ifndef GT2_NEIGHBORHOOD_BUILDER_H
#define GT2_NEIGHBORHOOD_BUILDER_H

#include <genetrail2/core/macros.h>

#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/SparseMatrix.h>

namespace GeneTrail
{
	/**
	 * This class constructs a sparse neighborhood graph
	 * from a matrix of datapoints. This is useful for
	 * subsequently clustering the most similar datapoints
	 * using e.g. the METISClusterer
	 */
	class GT2_EXPORT NeighborhoodBuilder
	{
		public:
			/**
			 * Constructs the actual neighborhood graph of the matrix rows.
			 *
			 * @param mat The data matrix. Currently rows are expected to hold
			 *            the variables. The columns represent the samples.
			 * @returns a neighborhood graph constructed in the following fashion:
			 *          for each variable the k most similar datapoints are chosen.
			 *          This leads to a maximum number of edges of 2 * |V| * k
			 */
			SparseMatrix build(DenseMatrix mat) const;

			/**
			 * Set the number of neighbors to k
			 */
			void setNumNeighbors(int k);

		private:
			typedef Eigen::Triplet<SparseMatrix::value_type> T;

			int k_;
			void insert_(T, std::vector<T>::iterator it) const;
			SparseMatrix buildMatrix_(std::vector<T>& entries, const DenseMatrix& mat) const;
	};
}

#endif // GT2_NEIGHBORHOOD_BUILDER_H

