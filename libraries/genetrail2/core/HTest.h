#ifndef GT2_CORE_HTEST_H
#define GT2_CORE_HTEST_H

#include <boost/math/distributions.hpp>

namespace GeneTrail
{
	namespace HTest
	{
		template<template<typename T> class HTest, typename value_type, typename InputIterator1, typename InputIterator2>
		value_type test(HTest<value_type>& t, const InputIterator1& first_begin, const InputIterator1& first_end, const InputIterator2& second_begin, const InputIterator2& second_end){
			return t.test(first_begin, first_end, second_begin, second_end);
		}

		template<template<typename T> class HTest, typename value_type>
		value_type twoSidedPValue(HTest<value_type>& t, const value_type& score) {
			return 2.0 * boost::math::cdf(boost::math::complement(t.distribution(), fabs(score)));
		}

		template<template<typename T> class HTest, typename value_type>
		std::pair<value_type, value_type> confidenceInterval(HTest<value_type>& t, const value_type& alpha){
			return t.confidenceInterval(alpha);
		}

		template<template<typename T> class HTest, typename value_type, typename InputIterator>
		value_type test(HTest<value_type>& t, const InputIterator& begin, const InputIterator& end){
			return t.test(begin, end);
		}

		template<template<typename T> class HTest, typename value_type>
		value_type lowerTailedPValue(HTest<value_type>& t, const value_type& score) {
			return boost::math::cdf(t.distribution(), score);
		}

		template<template<typename T> class HTest, typename value_type>
		value_type upperTailedPValue(HTest<value_type>& t, const value_type& score) {
			return boost::math::cdf(boost::math::complement(t.distribution(), score));
		}
	}
}

#endif //GT2_CORE_HTEST_H

