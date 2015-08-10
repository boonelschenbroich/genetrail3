#ifndef GT2_CORE_PERMUTATION_TEST_H
#define GT2_CORE_PERMUTATION_TEST_H

#include "macros.h"
#include "DenseMatrix.h"
#include "DenseMatrixSubset.h"
#include "EnrichmentAlgorithm.h"
#include "MatrixHTest.h"

#include <algorithm>
#include <vector>
#include <utility>
#include <random>
#include <iostream>

namespace GeneTrail
{
	template <typename value_type> struct TestResult
	{
		TestResult(const Category* cat, value_type s, size_t ss)
		    : category(cat), score(s), counter(0), sampleSize(ss), enriched(false)
		{
		}

		const Category* category;
		value_type score;
		size_t counter;
		size_t sampleSize;
		bool enriched;

		double computePValue(size_t permutations) const
		{
			// Here we add a pseudo count to avoid p-values of 0.
			// Reference: Fewer permutations, more accurate P-values.
			// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2687965/

			if(permutations == 0){
				return 1;
			}

			return ((double)counter + 1) / ((double)permutations);
		}
	};

	namespace internal {
		template <typename value_type> class PermutationBase
		{
			public:
			using TestResults = std::vector<TestResult<value_type>>;
			using PValue = std::pair<std::string, double>;
			using PValues = std::vector<PValue>;

			protected:
			template <typename Iterator>
			void
			performSinglePermutation_(TestResults& tests,
			                          const EnrichmentAlgorithmPtr& algorithm,
			                          const Category& c)
			{
			}
		};
	}

	template <typename value_type> class GT2_EXPORT PermutationTest : public internal::PermutationBase<value_type>
	{
		public:
		using Base = internal::PermutationBase<value_type>;
		using TestResults = typename Base::TestResults;
		using PValues = typename Base::PValues;

		template <typename InputIterator>
		PermutationTest(const TestResults& tests, InputIterator begin,
		                InputIterator end, size_t permutations)
		    : tests_(tests), permutations_(permutations), names_(begin, end)
		{
			// Sort vector of TestResults
			// This makes computation faster
			std::sort(tests_.begin(), tests_.end(),
			     [](const TestResult<value_type>& a,
			        const TestResult<value_type>& b)
			         ->bool { return a.sampleSize < b.sampleSize; });
		}

		PValues computePValue(const EnrichmentAlgorithmPtr& algorithm)
		{
			PValues pvalues;
			pvalues.reserve(tests_.size());

			for(size_t i = 0; i < permutations_; ++i) {
				std::cout << "INFO: Running - Permutation test " << i << "/"
				          << permutations_ << std::endl;
				// As we don't want to recompute any values,
				// we save them here.
				shuffle_();

				size_t currentSampleSize = -1;
				value_type currentScore = 0.0;
				Category c;
				for(size_t i = 0; i < tests_.size(); ++i) {
					// Check if we need to compute new values
					if(tests_[i].sampleSize != currentSampleSize) {
						c = Category("", names_.begin(),
						             names_.end() + currentSampleSize);
						currentSampleSize = tests_[i].sampleSize;
						currentScore = algorithm->computeEnrichmentScore(c);
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

			for(const auto& test : tests_) {
				pvalues.emplace_back(test.category->name(),
				                     test.computePValue(permutations_));
			}

			return pvalues;
		}

		private:
		void shuffle_()
		{
			std::shuffle(names_.begin(), names_.end(), twister_);
		}

		TestResults tests_;
		Category category_;
		size_t permutations_;
		std::mt19937 twister_;
		std::vector<std::string> names_;
	};

	template <typename value_type> class SamplePermutationTest : public internal::PermutationBase<value_type>
	{
		public:
		using Base = internal::PermutationBase<value_type>;
		using TestResults = typename Base::TestResults;
		using PValues = typename Base::PValues;

		SamplePermutationTest(TestResults tests, const DenseMatrix& data,
		                size_t permutations, size_t sample_size, size_t reference_size, const std::string& method)
			: tests_(tests),
			  permutations_(permutations),
			  data_(data),
			  sample_size_(sample_size),
			  reference_size_(reference_size),
			  method_(method)
		{
			//TODO: Initalize the random number generator.
			//      Allow saving seeds, etc.

			// Sort vector of TestResults
			// This makes computation faster
			std::sort(tests_.begin(), tests_.end(),
			     [](const TestResult<value_type>& a,
			        const TestResult<value_type>& b)
			         ->bool { return a.sampleSize < b.sampleSize; });
		}

		PValues computePValue(const EnrichmentAlgorithmPtr& algorithm)
		{
			MatrixHTest<DenseMatrixSubset> scoring_;
			PValues pvalues;
			pvalues.reserve(tests_.size());
			std::vector<size_t> indices(data_.cols());
			std::iota(indices.begin(), indices.end(), static_cast<size_t>(0));

			for(size_t i = 0; i < permutations_; ++i) {
				std::cout << "INFO: Running - Permutation test " << i << "/"
				          << permutations_ << std::endl;
				std::shuffle(indices.begin(), indices.end(), twister_);

				auto reference = DenseMatrixSubset::createColSubset(&data_, indices.begin(), indices.begin() + reference_size_);
				auto sample    = DenseMatrixSubset::createColSubset(&data_, indices.begin() + reference_size_, indices.end());
				auto gene_set = scoring_.test(reference, sample, method_);

				value_type currentScore = 0.0;
				for(auto& test : tests_) {
					currentScore =
					    algorithm->computeEnrichmentScore(*test.category);
					if(test.enriched) {
						if(test.score <= currentScore) {
							++test.counter;
						}
					} else {
						if(test.score >= currentScore) {
							++test.counter;
						}
					}
				}
			}

			for(auto&& test : tests_) {
				pvalues.emplace_back(test.category->name(),
				                     test.computePValue(permutations_));
			}

			return pvalues;
		}

		private:
		TestResults tests_;
		size_t permutations_;
		DenseMatrix data_;
		size_t sample_size_;
		size_t reference_size_;
		std::string method_;
		std::mt19937 twister_;
	};
}

#endif // GT2_CORE_PERMUTATION_TEST_H


