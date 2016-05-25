/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2014 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_MULTIPRECISION_H
#define GT2_CORE_MULTIPRECISION_H

#include <genetrail2/core/config.h>

#if defined GENETRAIL2_HAS_MPFR
#include <boost/multiprecision/mpfr.hpp>
namespace GeneTrail
{
	template <int N>
	using big_float_tpl =
	    boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<N>>;
}
#elif defined GENETRAIL2_HAS_GMP
#include <boost/multiprecision/gmp.hpp>
namespace GeneTrail
{
	template <int N>
	using big_float_tpl =
	    boost::multiprecision::number<boost::multiprecision::gmp_float<N>>;
}
#else
#include <boost/multiprecision/cpp_dec_float.hpp>
namespace GeneTrail
{
	template <int N>
	using big_float_tpl =
	    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<N>>;
}
#endif

namespace GeneTrail
{
	using big_float_50   = big_float_tpl<50>;
	using big_float_100  = big_float_tpl<100>;
	using big_float_500  = big_float_tpl<500>;
	using big_float_1000 = big_float_tpl<1000>;
	using big_float      = big_float_50;
}

#endif // GT2_CORE_MULTIPRECISION_H

