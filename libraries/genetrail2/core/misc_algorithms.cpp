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
