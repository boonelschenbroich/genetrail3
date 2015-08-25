#ifndef GT2_MISC_ALGORITHMS_H
#define GT2_MISC_ALGORITHMS_H

#include <algorithm>
#include <vector>

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

	template <typename InputIterator, typename Compare>
	void sort_permutation(std::vector<size_t>& p, InputIterator begin,
	                      InputIterator end, Compare compare)
	{
		p.resize(std::distance(begin, end));
		std::iota(p.begin(), p.end(), 0);
		std::sort(p.begin(), p.end(),
		          [&](int i, int j) { return compare(begin[i], begin[j]); });
	}

	template <typename InputIterator, typename Compare>
	std::vector<size_t> sort_permutation(InputIterator begin, InputIterator end,
	                                     Compare compare)
	{
		std::vector<size_t> p;
		sort_permutation(p, begin, end, compare);

		return p;
	}
}

#endif // GT2_MISC_ALGORITHMS_H
