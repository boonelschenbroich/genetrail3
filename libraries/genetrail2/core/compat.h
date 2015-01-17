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
