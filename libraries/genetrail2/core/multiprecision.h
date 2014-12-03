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

