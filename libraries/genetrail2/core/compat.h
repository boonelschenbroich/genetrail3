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
#ifndef GENETRAIL2_COMPAT_H
#define GENETRAIL2_COMPAT_H

#include <genetrail2/core/config.h>

// Add a simplified implementation for C++ 14's
// std::make_unique. We put it in the standard namespace
// so it can be replaced by the proper implementation.
#ifndef GT2_HAS_MAKE_UNIQUE
#include <memory>

namespace std
{
	template<typename T, typename... Ts>
	unique_ptr<T> make_unique(Ts&&... params) {
		return unique_ptr<T>(new T(std::forward<Ts>(params)...));
	}
}

#endif //GT2_HAS_MAKE_UNIQUE

#endif //GENETRAIL2_COMPAT_H
