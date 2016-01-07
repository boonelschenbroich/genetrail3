/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#include "misc_algorithms.h"

namespace GeneTrail
{
	namespace internal
	{
		void invert_cycle(size_t start, const std::vector<size_t>& perm,
		                  std::vector<size_t>& inv_perm)
		{
			auto prev_pos = start;
			auto current_pos = perm[start];

			for(; current_pos != start;
			    prev_pos = current_pos, current_pos = perm[current_pos]) {
				inv_perm[current_pos] = prev_pos;
			}

			inv_perm[current_pos] = prev_pos;
		}
	}

	void invert_permutation(const std::vector<size_t>& perm,
	                        std::vector<size_t>& inv_perm)
	{
		// Fill the inverse permutation with a place holder and allocate
		// enough space.
		inv_perm.assign(perm.size(), perm.size());

		for(size_t i = 0; i < perm.size(); ++i) {
			if(inv_perm[i] == perm.size()) {
				internal::invert_cycle(i, perm, inv_perm);
			}
		}
	}

	std::vector<size_t> invert_permutation(const std::vector<size_t>& perm)
	{
		std::vector<size_t> result;

		invert_permutation(perm, result);

		return result;
	}
}
