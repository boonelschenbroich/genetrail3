/*
 * GeneTrail2 - An efficient library for interpreting genetic data
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
#ifndef GT2_MISC_ALGORITHMS_H
#define GT2_MISC_ALGORITHMS_H

#include "macros.h"

#include <algorithm>
#include <numeric>
#include <vector>

namespace GeneTrail
{

/**
 * This method is the same as std::shuffle, however it only shuffles a part of
 * the input range. This is useful for e.g. computing permutation tests or
 * sampling a fixed where only a small portion of the input range needs to be
 * returned.
 *
 * Note that although the range [until, last) gets modified, it is not
 * guaranteed to be fully shuffled. This only holds for [first, until).
 *
 * @param first Start of the range that should be shuffled
 * @param until End of the range that should be shuffled and start of the range
 *              from which elements can be taken.
 * @param last  End of the range from which elements can be taken.
 */
template <typename RandomAccessIterator, typename RandomNumberGenerator>
void partial_shuffle(RandomAccessIterator first, RandomAccessIterator until,
                     RandomAccessIterator last, RandomNumberGenerator&& gen)
{
	if(first == until) {
		return;
	}

	if(until == last) {
		std::shuffle(first, last, std::forward<RandomNumberGenerator>(gen));
		return;
	}

	using uid = std::uniform_int_distribution<size_t>;
	using pt = uid::param_type;

	uid dist;
	auto d = (last - first) - 1;
	for(RandomAccessIterator i = first; i != until; ++i, --d) {
		// Generate a random index in the range [0, remaining_size - 1]
		std::iter_swap(i, i + dist(gen, pt(0, d)));
	}
}

	GT2_EXPORT void invert_permutation(const std::vector<size_t>& perm,
	                        std::vector<size_t>& inv_perm);

	GT2_EXPORT std::vector<size_t> invert_permutation(const std::vector<size_t>& perm);

	template <typename InputIterator, typename Compare>
	void sort_permutation(std::vector<size_t>& p, InputIterator begin,
	                      InputIterator end, Compare compare)
	{
		p.resize(std::distance(begin, end));
		std::iota(p.begin(), p.end(), static_cast<std::size_t>(0));
		std::sort(p.begin(), p.end(),
		          [&](std::size_t i, std::size_t j) { return compare(begin[i], begin[j]); });
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
