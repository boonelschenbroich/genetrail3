#ifndef GT2_MISC_ALGORITHMS_H
#define GT2_MISC_ALGORITHMS_H

#include "macros.h"

#include <algorithm>
#include <vector>

namespace GeneTrail
{
	GT2_EXPORT void invert_permutation(const std::vector<size_t>& perm,
	                        std::vector<size_t>& inv_perm);

	GT2_EXPORT std::vector<size_t> invert_permutation(const std::vector<size_t>& perm);

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
