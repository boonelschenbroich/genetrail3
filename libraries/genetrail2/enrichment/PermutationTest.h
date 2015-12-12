#ifndef GT2_CORE_PERMUTATION_TEST_H
#define GT2_CORE_PERMUTATION_TEST_H

#include "EnrichmentAlgorithm.h"
#include "EnrichmentResult.h"

#include <genetrail2/core/macros.h>
#include <genetrail2/core/misc_algorithms.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixSubset.h>
#include <genetrail2/core/MatrixHTest.h>

#include <boost/iterator/counting_iterator.hpp>

#include <algorithm>
#include <vector>
#include <utility>
#include <random>
#include <iostream>

namespace GeneTrail
{
namespace internal
{
template <typename value_type> class PermutationBase
{
  protected:
	double computePValue_(size_t permutations, size_t counter) const
	{
		// Here we add a pseudo count to avoid p-values of 0.
		// Reference: Fewer permutations, more accurate P-values.
		// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2687965/

		if(permutations == 0) {
			return 1;
		}

		return ((double)counter + 1) / ((double)permutations);
	}

	void sortResults_(EnrichmentResults& results)
	{
		// Sort vector of EnrichmentResults
		// This makes computation faster
		std::sort(
		    results.begin(), results.end(),
		    [](const EnrichmentResultPtr& a, const EnrichmentResultPtr& b) {
			    return a->hits < b->hits;
			});
	}

	void updateCounter_(const EnrichmentResultPtr& test, size_t& counter,
	                    double score)
	{
		if(test->enriched) {
			if(test->score <= score) {
				++counter;
			}
		} else {
			if(test->score >= score) {
				++counter;
			}
		}
	}

	void updatePValues_(EnrichmentResults& tests,
	                    const std::vector<size_t>& counter, size_t permutations)
	{
		for(size_t i = 0; i < tests.size(); ++i) {
			tests[i]->pvalue = computePValue_(permutations, counter[i]);
		}
	}
};
}

template <typename value_type>
class GT2_EXPORT RowPermutationTest
    : public internal::PermutationBase<value_type>
{
  public:
	static std::unique_ptr<RowPermutationTest>
	IndexBased(const Scores& s, size_t permutations, uint64_t randomSeed)
	{
		return std::unique_ptr<RowPermutationTest>(
		    new RowPermutationTest(boost::counting_iterator<size_t>(0),
		                           boost::counting_iterator<size_t>(s.size()),
		                           permutations, randomSeed));
	}

	static std::unique_ptr<RowPermutationTest>
	CategoryBased(const Scores& s, size_t permutations, uint64_t randomSeed)
	{
		return std::unique_ptr<RowPermutationTest>(new RowPermutationTest(
		    s.indices().begin(), s.indices().end(), permutations, randomSeed));
	}

	void computePValue(const EnrichmentAlgorithmPtr& algorithm,
	                   EnrichmentResults& tests)
	{
		if(tests.empty()) {
			return;
		}

		this->sortResults_(tests);

		std::vector<size_t> counter(tests.size());
		for(size_t i = 0; i < permutations_; ++i) {
			std::cout << "INFO: Running - Permutation test " << (i + 1) << "/"
			          << permutations_ << std::endl;

			performSinglePermutation_(algorithm, tests, counter);
		}

		this->updatePValues_(tests, counter, permutations_);
	}

  private:
	template <typename InputIterator>
	RowPermutationTest(InputIterator begin, InputIterator end,
	                   size_t permutations, uint64_t randomSeed)
	    : category_(nullptr),
	      permutations_(permutations),
	      twister_(randomSeed),
	      indices_(begin, end),
	      tmp_indices_(std::distance(begin, end))
	{
	}

	std::tuple<double, double>
	computeEnrichmentScore_(const EnrichmentAlgorithmPtr& algorithm,
	                        size_t currentSampleSize)
	{
		auto start = indices_.begin();
		auto stop = indices_.begin() + currentSampleSize;

		if(algorithm->supportsIndices()) {
			return algorithm->computeEnrichmentScore(start, stop);
		} else {
			category_.replaceAll(start, stop);
			return algorithm->computeEnrichmentScore(category_);
		}
	}

	void performSinglePermutation_(const EnrichmentAlgorithmPtr& algorithm,
	                               const EnrichmentResults& tests,
	                               std::vector<size_t>& counter)
	{
		size_t currentSampleSize = 0;
		value_type currentScore = 0.0;

		// Shuffle the indices. The tests are sorted
		// by the number of hits for every category, so we only
		// need to shuffle tests.back()->hits many.
		shuffle_(tests.back()->hits);

		for(size_t i = 0; i < tests.size(); ++i) {
			// Check if the sampleSize has changed. As the tests_ vector is
			// sorted we can use one running sum value for all categories of
			// the same size.
			if(tests[i]->hits != currentSampleSize) {
				presort_(currentSampleSize, tests[i]->hits);
				currentSampleSize = tests[i]->hits;

				currentScore = std::get<0>(
				    computeEnrichmentScore_(algorithm, currentSampleSize));
			}

			this->updateCounter_(tests[i], counter[i], currentScore);
		}
	}

	void shuffle_(size_t n)
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
	void presort_(size_t a, size_t b)
	{
		// Sort [a, b)
		std::sort(indices_.begin() + a, indices_.begin() + b);

		// Merge [0, a) and [a, b) into a temporary vector.
		// Note that inplace_merge would reallocate a new vector
		// every time it is called.
		std::merge(indices_.begin(), indices_.begin() + a, indices_.begin() + a,
		           indices_.begin() + b, tmp_indices_.begin());

		// Copy the sorted range [0, b) from the temporary vector
		// to the names_ vector
		std::copy(tmp_indices_.begin(), tmp_indices_.begin() + b,
		          indices_.begin());
	}

	Category category_;
	size_t permutations_;
	std::mt19937_64 twister_;
	std::vector<size_t> indices_;
	std::vector<size_t> tmp_indices_;
};

template <typename value_type>
class ColumnPermutationTest : public internal::PermutationBase<value_type>
{
  public:
	ColumnPermutationTest(const DenseMatrix& data, size_t permutations,
	                      size_t reference_size, MatrixHTests method,
	                      uint64_t randomSeed)
	    : permutations_(permutations),
	      data_(data),
	      reference_size_(reference_size),
	      method_(method),
	      twister_(randomSeed)
	{
		assert(rowNamesStrictlySorted_(data));
	}

	void initScoring_()
	{
		std::vector<size_t> row_db_indices(data_.rows());

		EntityDatabase::global->transform(data_.rowNames().begin(),
		                                  data_.rowNames().end(),
		                                  row_db_indices.begin());

		scoring.setRowDBIndices(row_db_indices);
	}

	void computePValue(const EnrichmentAlgorithmPtr& algorithm,
	                   EnrichmentResults& tests)
	{
		positions_.clear();

		initScoring_();

		std::vector<size_t> counter(tests.size());
		std::vector<size_t> column_indices(data_.cols());
		std::iota(column_indices.begin(), column_indices.end(),
		          static_cast<size_t>(0));

		for(size_t i = 1; i <= permutations_; ++i) {
			std::cout << "INFO: Running - Permutation test " << i << "/"
			          << permutations_ << std::endl;

			if(algorithm->supportsIndices()) {
				performSinglePermutationIndices_(algorithm, tests, counter,
				                                 column_indices);
			} else {
				performSinglePermutation_(algorithm, tests, counter,
				                          column_indices);
			}
		}

		this->updatePValues_(tests, counter, permutations_);
	}

  private:
	bool rowNamesStrictlySorted_(const DenseMatrix& data)
	{
		Matrix::index_type i = 0;

		for(auto j = i + 1; j < data.rows(); ++j, ++i) {
			if(data.rowName(i) >= data.rowName(j)) {
				return false;
			}
		}

		return true;
	};

	Scores computeScores_(std::vector<size_t>::iterator begin,
	                      std::vector<size_t>::iterator end)
	{
		std::shuffle(begin, end, twister_);

		auto mid = begin + reference_size_;

		auto ref = DenseMatrixSubset::createColSubset(&data_, begin, mid);
		auto sam = DenseMatrixSubset::createColSubset(&data_, mid, end);

		return Scores(scoring.test(method_, ref, sam));
	}

	void setupPositons_(const Scores& scores, const EnrichmentResults& tests)
	{
		positions_.resize(tests.size());

		for(size_t i = 0; i < tests.size(); ++i) {
			positions_[i] = scores.subsetIndices(*tests[i]->category);
		}
	}

	Scores updateScores(const EnrichmentResults& tests,
	                    std::vector<size_t>& column_indices)
	{
		Scores scores =
		    computeScores_(column_indices.begin(), column_indices.end());

		// We first need to sort the scores by index, as we know the
		// position of the category genes in the sorted list.
		scores.sortByIndex();

		// If we do not yet know the positions compute them.
		// We can only do this here as we do not have the scores
		// available earlier.
		if(positions_.empty()) {
			setupPositons_(scores, tests);
		}

		// We now obtain the permutation of the genes that is used
		// for sorting the scores by value.
		sort_permutation(permutation_, scores.scores().begin(),
		                 scores.scores().end(), std::less<double>());

		// After that we invert the permutation so that we
		// can use it as a lookup table to get the new position.
		invert_permutation(permutation_, inv_permutation_);

		return scores;
	}

	double computeEnrichmentScore_(const EnrichmentAlgorithmPtr& algorithm,
	                               const EnrichmentResults& tests, size_t i)
	{
		// Setup the vector that will hold the positions of the category
		// genes in the new score vector.
		intersection_.resize(tests[i]->hits);

		// Fill the intersection vector using the inverse permutation
		std::transform(positions_[i].begin(), positions_[i].end(),
		               intersection_.begin(),
		               [&](size_t entry) { return inv_permutation_[entry]; });

		// Sort the intersection_ vector
		std::sort(intersection_.begin(), intersection_.end());

		// Compute the enrichment score for this permutation
		return std::get<0>(algorithm->computeEnrichmentScore(
		    intersection_.begin(), intersection_.end()));
	}

	void performSinglePermutation_(const EnrichmentAlgorithmPtr& algorithm,
	                               const EnrichmentResults& tests,
	                               std::vector<size_t>& counter,
	                               std::vector<size_t>& column_indices)
	{
		Scores scores =
		    computeScores_(column_indices.begin(), column_indices.end());

		algorithm->setScores(scores);

		for(size_t i = 0; i < tests.size(); ++i) {
			auto score = std::get<0>(
			    algorithm->computeEnrichmentScore(*tests[i]->category));

			this->updateCounter_(tests[i], counter[i], score);
		}
	}

	void performSinglePermutationIndices_(
	    const EnrichmentAlgorithmPtr& algorithm, const EnrichmentResults& tests,
	    std::vector<size_t>& counter, std::vector<size_t>& column_indices)
	{
		Scores scores = updateScores(tests, column_indices);

		// Now pass the scores to the algorithm
		algorithm->setScores(scores);

		for(size_t i = 0; i < tests.size(); ++i) {
			auto score = computeEnrichmentScore_(algorithm, tests, i);

			// Update the counter with the newly computed value
			this->updateCounter_(tests[i], counter[i], score);
		}
	}

	size_t permutations_;
	DenseMatrix data_;
	size_t reference_size_;
	MatrixHTests method_;
	std::mt19937 twister_;

	MatrixHTest scoring;
	std::vector<size_t> permutation_;
	std::vector<size_t> inv_permutation_;
	std::vector<size_t> intersection_;
	std::vector<std::vector<size_t>> positions_;
};
}

#endif // GT2_CORE_PERMUTATION_TEST_H
