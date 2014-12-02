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
		const WeightedGeneSetEnrichmentAnalysis<value_type>& wgsea_;
		const size_t sizeOfList_;
		const size_t permutations_;
		PValues pvalues_;
		std::mt19937 twister_;
		std::vector<size_t> indices_;

		// Workspace for sorting the indices, its content
		// is irrelevant and it is only a member as we want to avoid
		// frequent reallocations.
		std::vector<size_t> tmp_indices_;

		/**
		 * Randomly select n indices from our index list
		 */
		void shuffleIndices_(const size_t n)
		{
			for(size_t i = 0; i < n; ++i) {
				std::swap(indices_[i],
				          indices_[i + twister_() % (indices_.size() - i)]);
			}
		}

		/**
		 * This sorts the indices vector until position b.
		 * It is assumed, that the vector is sorted up to position a.
		 */
		void sortIndices_(size_t a, size_t b)
		{
			// Sort [a, b)
			std::sort(indices_.begin() + a, indices_.begin() + b);

			// Merge [0, a) and [a, b) into a temporary vector.
			// Note that inplace_merge would reallocate a new vector
			// every time it is called.
			std::merge(indices_.begin(), indices_.begin() + a,
			           indices_.begin() + a, indices_.begin() + b,
			           tmp_indices_.begin());

			// Copy the sorted range [0, b) from the temporary vector
			// to the indices_ vector
			std::copy(tmp_indices_.begin(), tmp_indices_.begin() + b,
			          indices_.begin());
		}

		/**
		 * Perform a single permutation for all catagories
		 */
		void performSinglePermutation_()
		{
			size_t currentSampleSize = 0;
			value_type currentScore = 0.0;

			// Shuffle the first few indices
			shuffleIndices_(tests_.back().sampleSize);

			for(auto& t : tests_) {
				// Check if the sampleSize has changed. As the tests_ vector is
				// sorted we can use one running sum value for all categories of
				// the same size.
				if(t.sampleSize != currentSampleSize) {
					// Our sampleSize changed, we reorder the indices
					sortIndices_(currentSampleSize, t.sampleSize);

					// Update currentSampleSize and currentScore
					currentSampleSize = t.sampleSize;
					currentScore = wgsea_.computeRunningSum(
					    indices_.begin(), indices_.begin() + currentSampleSize);
				}

				// Check whether the category is enriched and then
				// increase the counter if the score is non-significant
				if(t.enriched) {
					if(t.score <= currentScore) {
						++t.counter;
					}
				} else {
					if(t.score >= currentScore) {
						++t.counter;
					}
				}
			}
		}

		public:
		WGSEAPermutationTest(
		    const TestResults& tests,
		    const WeightedGeneSetEnrichmentAnalysis<value_type>& wgsea,
		    const size_t sizeOfList, const size_t permutations)
		    : tests_(tests),
		      wgsea_(wgsea),
		      sizeOfList_(sizeOfList),
		      permutations_(permutations),
		      twister_{std::random_device{}()},
		      indices_(sizeOfList),
		      tmp_indices_(sizeOfList)
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

		const PValues& computePValue()
		{
			for(size_t k = 1; k <= permutations_; ++k) {
				std::cout << "INFO: Running - Permutation test " << k << "/"
				          << permutations_ << "\n";

				performSinglePermutation_();
			}

			for(auto& t : tests_) {
				pvalues_.emplace_back(t.name, t.computePValue(permutations_));
			}

			return pvalues_;
		}
	};
}

#endif // GT2_CORE_WGSEA_PERMUTATION_TEST_H

