#ifndef GT2_CORE_GENE_LABEL_PERMUTATION_TEST_H
#define GT2_CORE_GENE_LABEL_PERMUTATION_TEST_H

#include <algorithm>

#include "macros.h"

namespace GeneTrail
{
	template<typename value_type, typename InputIterator>
	class GT2_EXPORT GeneLabelPermutationTest
	{
		public:

		GeneLabelPermutationTest() = default;

		GeneLabelPermutationTest(InputIterator begin, InputIterator end, unsigned int permutations = 1000)
		: begin_(begin), end_(end), permutations_(permutations)
		{}

		template<typename Comparator, typename Function>
		value_type computePValue(unsigned int sampleSize, value_type score, Comparator comp, Function f){
			unsigned int n = 0;
			for(int i=0; i< permutations_; ++i){
				std::random_shuffle(begin_, end_);
				if(comp(score, f(begin_, begin_ + sampleSize))){
					++n;
				}
			}
			return ((value_type)n)/((value_type)permutations_);
		}

		template <typename Function>
		value_type computeLowerTailedPValue(unsigned int sampleSize, value_type score, Function f){
			return computePValue(sampleSize, score, less, f);
		}

		template <typename Function>
		value_type computeUpperTailedPValue(unsigned int sampleSize, value_type score, Function f){
			return computePValue(sampleSize, score, greater, f);
		}

		template <typename Function>
		value_type twoSidedPValue(unsigned int sampleSize, value_type score, Function f){
			return computePValue(sampleSize, score, extremer, f);
		}

		private:
		static bool less(value_type a, value_type b){
			return a >= b;
		}

		static bool greater(value_type a, value_type b){
			return a <= b;
		}

		static bool extremer(value_type a, value_type b){
			return  a <= b || a >= b ;
		}

		InputIterator begin_;
		InputIterator end_;
		unsigned int permutations_;
	};
}

#endif //GT2_CORE_GENE_LABEL_PERMUTATION_TEST_H

