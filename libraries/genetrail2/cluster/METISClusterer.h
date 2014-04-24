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
 */

#ifndef GT2_METIS_CLUSTERER_H
#define GT2_METIS_CLUSTERER_H

#include <genetrail2/core/macros.h>

#include <genetrail2/core/SparseMatrix.h>
#include "SparseClusterer.h"

#include <vector>

namespace GeneTrail
{
	/**
	 * Computes a (for now balanced) partitioning of the input graph.
	 *
	 * This class uses the METIS library for computing the graph partition.
	 */
	class GT2_EXPORT METISClusterer : public SparseClusterer
	{
		public:
			/**
			 * The available algorithms for computing
			 * the partition.
			 */
			enum Algorithm {
				Recursive,
				KWay
			};

			/**
			 * Create a clusterer
			 */
			explicit METISClusterer(Algorithm alg = Recursive);

			/**
			 * Convenience function @see computeClusters(SparseMatrix::SMatrix& matrix)
			 */
			void computeClusters(SparseMatrix& matrix);

			/**
			 * Compute a partitioning of the graph represented by a sparse eigen matrix.
			 *
			 * @warning The matrix _must_ be symmetric.
			 */
			void computeClusters(SparseMatrix::SMatrix& matrix);

			/**
			 * Set the algorithm used for computing the partition.
			 */
			void setAlgorithm(Algorithm alg);

			/**
			 * Set the number of clusters that should be computed.
			 */
			void setNumClusters(int k);

			/**
			 * Obtain the assignement of vertices to clusters.
			 *
			 * This object is only valid if @see computeClusters was called previously
			 */
			const std::vector<int>& grouping() const;

		private:
			Algorithm algorithm_;

			int k_;
			std::vector<int> grouping_;
	};
}

#endif // GT2_METIS_CLUSTERER_H

