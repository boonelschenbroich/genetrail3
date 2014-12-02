#ifndef GT2_CORE_WGSEA_PERMUTATION_TEST_H
#define GT2_CORE_WGSEA_PERMUTATION_TEST_H

#include <algorithm>
#include <vector>
#include <utility>
#include <random>
#include <iostream>

#include "macros.h"
#include "PermutationTest.h"
#include "WeightedGeneSetEnrichmentAnalysis.h"

namespace GeneTrail
{
	template <typename value_type> class GT2_EXPORT WGSEAPermutationTest
	{
		using PValue = std::pair<std::string, double>;
		using PValues = std::vector<PValue>;
		using TestResults = std::vector<TestResult<value_type>>;

		private:
		TestResults tests_;
		WeightedGeneSetEnrichmentAnalysis<value_type>& wgsea_;
		const size_t sizeOfList_;
		const size_t permutations_;
		PValues pvalues_;
		std::mt19937 twister_;
		std::vector<size_t> indices_;

		public:
		WGSEAPermutationTest(
		    TestResults& tests,
		    WeightedGeneSetEnrichmentAnalysis<value_type>& wgsea,
		    const size_t sizeOfList, const size_t permutations)
		    : tests_(tests),
		      wgsea_(wgsea),
		      sizeOfList_(sizeOfList),
		      permutations_(permutations),
		      twister_{std::random_device{}()},
		      indices_(sizeOfList)
		{
			pvalues_.reserve(tests.size());
			// Sort vector of TestResults
			// This makes computation faster
			sort(tests_.begin(), tests_.end(),
			     [](const TestResult<value_type>& a,
			        const TestResult<value_type>& b)
			         ->bool { return a.sampleSize < b.sampleSize; });

			std::iota(indices_.begin(), indices_.end(), 0);
		}

		void getRandomIndices(const size_t n)
		{
			for(size_t i = 0; i < n; ++i) {
				std::swap(indices_[i],
				          indices_[i + twister_() % (indices_.size() - i)]);
			}

			sort(indices_.begin(), indices_.begin() + n);
		}

		PValues computePValue()
		{
			for(size_t i = 0; i < permutations_; ++i) {
				std::cout << "INFO: Running - Permutation test " << i << "/"
				          << permutations_ << std::endl;
				// As we don't want to recompute any values,
				// we save them here.
				size_t currentSampleSize = 0;
				value_type currentScore = 0.0;
				for(size_t i = 0; i < tests_.size(); ++i) {
					// Check if we need to compute new values
					if(tests_[i].sampleSize != currentSampleSize) {
						currentSampleSize = tests_[i].sampleSize;
						getRandomIndices(currentSampleSize);
						currentScore = wgsea_.computeRunningSum(
						    indices_.begin(),
						    indices_.begin() + currentSampleSize);
					}
					if(tests_[i].enriched) {
						if(tests_[i].score <= currentScore) {
							++tests_[i].counter;
						}
					} else {
						if(tests_[i].score >= currentScore) {
							++tests_[i].counter;
						}
					}
				}
			}
			for(size_t i = 0; i < tests_.size(); ++i) {
				pvalues_.push_back(std::make_pair(
				    tests_[i].name, tests_[i].computePValue(permutations_)));
			}
			return pvalues_;
		}
	};
}

#endif // GT2_CORE_WGSEA_PERMUTATION_TEST_H

