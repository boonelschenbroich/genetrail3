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

#include "NeighborhoodBuilder.h"

#include <deque>
#include <vector>

#include <iostream>

namespace GeneTrail
{
	typedef Eigen::Triplet<SparseMatrix::value_type> T;

	bool triple_equal(const T& a, const T& b) {
		return (a.row() == b.row()) && (a.col() == b.col());
	}

	bool triple_less(const T& a, const T& b) {
		return a.row() < b.row() || ((a.row() == (b.row())) && a.col() < b.col());
	}

	SparseMatrix NeighborhoodBuilder::build(DenseMatrix mat) const
	{
		/*
		 * This matrix contains the k_ best neighbors for every vertex.
		 * The candidates are sorted in ascending order.
		 *
		 * The matrix is updated for entry as correlations are computed.
		 * This allows us to cut the computation time in half!
		 *
		 * We allocate twice the (worst-case) storage needed, as we have
		 * to symmetrize the matrix afterwards. Alternatively we could
		 * resize later, which might lead to a full copy.
		 */
		std::vector<T> entries(2 * k_ * mat.rows());

		std::vector<DenseMatrix::value_type> sd(mat.rows());

		DenseMatrix::Vector mu = mat.matrix().rowwise().mean();

		for(unsigned int i = 0; i < mat.rows(); ++i) {
			mat.row(i) = mat.row(i).array() - mu[i];
			sd[i] = mat.row(i).norm();
		}

		for(unsigned int i = 0; i < mat.rows(); ++i) {
			for(unsigned int j = i + 1; j < mat.rows(); ++j) {
				double cov = mat.row(i).dot(mat.row(j));
				cov = fabs(cov) / (sd[i] * sd[j]);

				// Insert into the entries vector
				insert_(T(i, j, cov), entries.begin() + i * k_);
				insert_(T(i, j, cov), entries.begin() + j * k_);
			}
		}

		return buildMatrix_(entries, mat);
	}

	SparseMatrix NeighborhoodBuilder::buildMatrix_(std::vector<T>& entries, const DenseMatrix& mat) const
	{
		// Make sure, that there are no duplicate entries
		std::sort(entries.begin(), entries.begin() + k_ * mat.rows(), triple_less);
		auto new_end = std::unique(entries.begin(), entries.begin() + k_ * mat.rows(), triple_equal);

		// An iterator for inserting the transposed entries
		// This will also serve as the new end pointer of the
		// entries array.
		auto new_it = new_end;

		// Symmetrize the matrix
		for(auto old_it = entries.begin(); old_it != new_end; ++old_it, ++new_it) {
			*new_it = T(old_it->col(), old_it->row(), old_it->value());
		}

		// Build the temporary matrix for the results
		SparseMatrix result(mat.rowNames(), mat.rowNames());

		// Fill the matrix and convert to CCS
		// new_it points to the end of valid matrix entries.
		result.matrix().setFromTriplets(entries.begin(), new_it);
		result.matrix().makeCompressed();

		return result;
	}

	void NeighborhoodBuilder::setNumNeighbors(int k)
	{
		k_ = k;
	}

	// Bubble-sort the entry to the right position
	void NeighborhoodBuilder::insert_(T v, std::vector<T>::iterator it) const
	{
		// Is the new value better than one current value
		if(v.value() <= it->value()) {
			return;
		}

		*it = v;
		auto prev = it;
		++it;

		for(int i = 1; i < k_; ++i) {
			// Is this the right place?
			if(prev->value() < it->value()) {
				break;
			}

			//Nope, the new value has to be advanced to the next field
			std::swap(*prev, *it);
			++it;
			++prev;
		}
	}
}
