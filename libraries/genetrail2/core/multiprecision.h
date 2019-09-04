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
#include <boost/math/bindings/detail/big_lanczos.hpp>

#if defined GENETRAIL2_HAS_MPFR
#include <boost/multiprecision/mpfr.hpp>
namespace GeneTrail
{
namespace bmp = boost::multiprecision;
template <int N> using big_float_tpl = bmp::number<bmp::mpfr_float_backend<N>>;

using rational = bmp::mpq_rational;
}
#elif defined GENETRAIL2_HAS_GMP
#include <boost/multiprecision/gmp.hpp>
namespace GeneTrail
{
namespace bmp = boost::multiprecision;
template <int N> using big_float_tpl = bmp::number<bmp::gmp_float<N>>;

using rational = bmp::mpq_rational;
}
#else
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
namespace GeneTrail
{
namespace bmp = boost::multiprecision;
template <int N> using big_float_tpl = bmp::number<bmp::cpp_dec_float<N>>;

using rational = bmp::cpp_rational;
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

// Defining the lanczos approximation that we want for our precision. The
// default lookup from boost (version 1.65 and above) wants undefined_lanczos
// which is incredibly slow. (More details in github.com/boostorg/math/issues/247)
namespace boost{ namespace math{ namespace lanczos{
	// Define template for our own lanczos class that uses the closest existing
	// lanczos for the given type
        template <typename return_type, typename lanczos_type>
	struct gt2_lanczos_tpl{
                static return_type lanczos_sum(const return_type& z){ return lanczos_type::lanczos_sum(z); }
                static return_type lanczos_sum_expG_scaled(const return_type& z){ return lanczos_type::lanczos_sum_expG_scaled(z); }
                static return_type lanczos_sum_near_1(const return_type& z){ return lanczos_type::lanczos_sum_near_1(z); }
                static return_type lanczos_sum_near_2(const return_type& z){ return lanczos_type::lanczos_sum_near_2(z); }
                static return_type g(){ return lanczos_type::g(); }
        };

	// Specializations of the general lanczos struct. Now we can safely use
	// boost::math::binomial_coefficient() and it will choose one of the
	// lanczos specializations defined here.
        template<class Policy>
        struct lanczos<GeneTrail::big_float_50, Policy>{
                typedef gt2_lanczos_tpl<GeneTrail::big_float_50, lanczos31UDT> type;
        };
        template<class Policy>
        struct lanczos<GeneTrail::big_float_100, Policy>{
                typedef gt2_lanczos_tpl<GeneTrail::big_float_100, lanczos61UDT> type;
        };
        template<class Policy>
        struct lanczos<GeneTrail::big_float_500, Policy>{
                typedef gt2_lanczos_tpl<GeneTrail::big_float_500, lanczos61UDT> type;
        };
        template<class Policy>
        struct lanczos<GeneTrail::big_float_1000, Policy>{
                typedef gt2_lanczos_tpl<GeneTrail::big_float_1000, lanczos61UDT> type;
        };

}}}


#endif // GT2_CORE_MULTIPRECISION_H

