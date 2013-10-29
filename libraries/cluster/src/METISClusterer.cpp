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

#include "METISClusterer.h"

#include <metis.h>

#include <iostream>

namespace GeneTrail
{
	// METIS uses integer weights for the cut calculation. Remap our
	// floating point weights to the full range of edge weights.
	std::vector<idx_t> fixedPointWeights_(SparseMatrix::SMatrix& mat)
	{
		SparseMatrix::value_type max = -std::numeric_limits<SparseMatrix::value_type>::infinity();

		for(int k = 0; k < mat.outerSize(); ++k) {
			for (SparseMatrix::SMatrix::InnerIterator it(mat, k); it; ++it)
			{
				max = std::max(max, it.value());
			}
		}

		// As range we use the 16bit integers. This should
		// prevent possible overflows in METIS (and be precise
		// enough for our purposes)
		const idx_t imax = std::numeric_limits<int16_t>::max();

		std::vector<idx_t> result(mat.nonZeros());
		int i = 0;
		for(int k = 0; k < mat.outerSize(); ++k) {
			for (SparseMatrix::SMatrix::InnerIterator it(mat, k); it; ++it)
			{
				result[i] = std::round(imax * (it.value() / max));
				++i;
			}
		}

		return result;
	}

	METISClusterer::METISClusterer(METISClusterer::Algorithm alg)
		: algorithm_(alg),
		  k_(0)
	{
	}

	void METISClusterer::computeClusters(SparseMatrix& matrix)
	{
		computeClusters(matrix.matrix());
	}

	void METISClusterer::computeClusters(SparseMatrix::SMatrix& mat)
	{
		if(mat.rows() != mat.cols()) {
			//TODO: Proper error handling
			return;
		}

		grouping_.resize(mat.rows());

		//TODO: Evalutate whether this should be done by the user.
		//      This would allow passing const matrices
		mat.makeCompressed();

		idx_t nvtxs = mat.rows();
		idx_t ncon = 1;

		idx_t options[METIS_NOPTIONS];
		METIS_SetDefaultOptions(options);

		std::vector<idx_t> eweights = fixedPointWeights_(mat);

		idx_t objval = -1;

		switch(algorithm_) {
			case Recursive:
				METIS_PartGraphRecursive(&nvtxs, &ncon, mat.outerIndexPtr(), mat.innerIndexPtr(),
				                         NULL, NULL, &eweights[0], &k_, NULL, NULL, options,
				                         &objval, &grouping_[0]);
				break;
			case KWay:
				METIS_PartGraphKway(&nvtxs, &ncon, mat.outerIndexPtr(), mat.innerIndexPtr(),
				                    NULL, NULL, &eweights[0], &k_, NULL, NULL, options,
				                    &objval, &grouping_[0]);
				break;
		}
	}

	const std::vector< int >& METISClusterer::grouping() const
	{
		return grouping_;
	}

	void METISClusterer::setAlgorithm(METISClusterer::Algorithm alg)
	{
		algorithm_ = alg;
	}

	void METISClusterer::setNumClusters(int k)
	{
		k_ = k;
	}
}
