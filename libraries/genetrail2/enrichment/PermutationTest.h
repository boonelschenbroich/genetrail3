/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#ifndef GT2_ENRICHMENT_PERMUTATION_TEST_H
#define GT2_ENRICHMENT_PERMUTATION_TEST_H

#include "EnrichmentAlgorithm.h"
#include "EnrichmentResult.h"

#include <genetrail2/core/macros.h>
#include <genetrail2/core/misc_algorithms.h>
#include <genetrail2/core/DenseColumnSubset.h>
#include <genetrail2/core/DenseMatrix.h>
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
	void printStatus_(size_t i, size_t count)
	{
		std::cout << "INFO: Running - Permutation test " << (i + 1) << "/"
		          << count << std::endl;
	}

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
		    [](const EnrichmentResultPtr& a,
		       const EnrichmentResultPtr& b) { return a->hits < b->hits; });
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
			this->printStatus_(i, permutations_);

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
class ColumnPermutationBase : public internal::PermutationBase<value_type>
{
  public:
	ColumnPermutationBase(const DenseMatrix& data, size_t permutations,
	                      size_t reference_size, MatrixHTests method,
	                      uint64_t randomSeed, const EntityDatabase* db)
	    : permutations_(permutations),
	      data_(data),
	      reference_size_(reference_size),
	      method_(method),
	      twister_(randomSeed),
	      db_(db)
	{
	}

	void initScoring_()
	{
		std::vector<size_t> row_db_indices(data_.rows());

		db_->transform(data_.rowNames().begin(), data_.rowNames().end(),
		               row_db_indices.begin());

		// We need to make sure, that the EntityDatabase indices of
		// the data are sorted in strictly ascending order, as we
		// exploit this property later to find the genes belonging to a
		// category.
		assert(rowIndicesStrictlySorted_(row_db_indices));

		scoring.setRowDBIndices(row_db_indices);
	}

  protected:
	bool rowIndicesStrictlySorted_(const std::vector<size_t>& data)
	{
		Matrix::index_type i = 0;

		for(auto j = i + 1; j < data.size(); ++j, ++i) {
			if(data[i] >= data[j]) {
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

		auto ref = DenseColumnSubset(&data_, begin, mid);
		auto sam = DenseColumnSubset(&data_, mid, end);

		return Scores(scoring.test(method_, ref, sam));
	}

	void setupPositons_(const Scores& scores, const EnrichmentResults& tests)
	{
		positions_.resize(tests.size());

		for(size_t i = 0; i < tests.size(); ++i) {
			positions_[i] = scores.subsetIndices(*tests[i]->category);
		}
	}

	Scores updateLookupTables_(const EnrichmentResults& tests, Scores& scores,
	                           Order order)
	{
		// We first need to sort the scores by index, as we know the
		// position of the category genes in the sorted list.
		scores.sortByIndex();

		// If we do not yet know the positions compute them.
		// We can only do this here as we do not have the scores
		// available earlier.
		if(positions_.empty()) {
			setupPositons_(scores, tests);
		}

		// We now obtain the permution of the genes that is used
		// for sorting the scores by value.
		switch(order) {
			case Order::Decreasing:
				sort_permutation(permutation_, scores.scores().begin(),
				                 scores.scores().end(), std::less<double>());
			case Order::Increasing:
				sort_permutation(permutation_, scores.scores().begin(),
				                 scores.scores().end(), std::greater<double>());
		}

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
	const EntityDatabase* db_;
};

template <typename value_type>
class ColumnPermutationTest : public ColumnPermutationBase<value_type>
{
  public:
	ColumnPermutationTest(const DenseMatrix& data, size_t permutations,
	                      size_t reference_size, MatrixHTests method,
	                      uint64_t randomSeed, const EntityDatabase* db)
	    : ColumnPermutationBase<value_type>(data, permutations, reference_size,
	                                        method, randomSeed, db)
	{
	}

	void computePValue(const EnrichmentAlgorithmPtr& algorithm,
	                   EnrichmentResults& tests)
	{
		this->positions_.clear();

		this->initScoring_();

		std::vector<size_t> counter(tests.size());
		std::vector<size_t> column_indices(this->data_.cols());
		std::iota(column_indices.begin(), column_indices.end(),
		          static_cast<size_t>(0));

		for(size_t i = 0; i < this->permutations_; ++i) {
			this->printStatus_(i, this->permutations_);

			if(algorithm->supportsIndices()) {
				performSinglePermutationIndices_(algorithm, tests, counter,
				                                 column_indices);
			} else {
				performSinglePermutation_(algorithm, tests, counter,
				                          column_indices);
			}
		}

		this->updatePValues_(tests, counter, this->permutations_);
	}

	void performSinglePermutation_(const EnrichmentAlgorithmPtr& algorithm,
	                               const EnrichmentResults& tests,
	                               std::vector<size_t>& counter,
	                               std::vector<size_t>& column_indices)
	{
		Scores scores =
		    this->computeScores_(column_indices.begin(), column_indices.end());

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
		// Create a new permutation and compute new scores.
		Scores scores =
		    this->computeScores_(column_indices.begin(), column_indices.end());

		// Update all the required lookup tables needed for finding
		// category members quickly.
		this->updateLookupTables_(tests, scores, algorithm->getOrder());

		// Now pass the scores to the algorithm
		algorithm->setScores(scores);

		for(size_t i = 0; i < tests.size(); ++i) {
			auto score = this->computeEnrichmentScore_(algorithm, tests, i);

			// Update the counter with the newly computed value
			this->updateCounter_(tests[i], counter[i], score);
		}
	}
};

template <typename value_type>
class KSColumnPermutationTest : public ColumnPermutationBase<value_type>
{
  public:
	KSColumnPermutationTest(const DenseMatrix& data, size_t permutations,
	                        size_t reference_size, MatrixHTests method,
	                        uint64_t randomSeed, const EntityDatabase* db)
	    : ColumnPermutationBase<value_type>(data, permutations, reference_size,
	                                        method, randomSeed, db)
	{
	}

	void computePValue(const EnrichmentAlgorithmPtr& algorithm,
	                   EnrichmentResults& tests)
	{
		this->positions_.clear();

		this->initScoring_();

		std::vector<double> permuted_values(tests.size() * this->permutations_);
		std::vector<size_t> column_indices(this->data_.cols());
		std::iota(column_indices.begin(), column_indices.end(),
		          static_cast<size_t>(0));

		for(size_t i = 0; i < this->permutations_; ++i) {
			this->printStatus_(i, this->permutations_);

			if(algorithm->supportsIndices()) {
				performSinglePermutationIndices_(
				    i, algorithm, tests, permuted_values, column_indices);
			} else {
				performSinglePermutation_(i, algorithm, tests, permuted_values,
				                          column_indices);
			}
		}

		updatePValues_(tests, permuted_values);
	}

  protected:
	ptrdiff_t computeGreaterPos(const std::vector<double>& v,
	                            const std::vector<double>::const_iterator zero,
	                            double value)
	{
		return v.end() - std::lower_bound(zero, v.end(), value);
	}

	ptrdiff_t computeGreaterNeg(const std::vector<double>& v,
	                            const std::vector<double>::const_iterator zero,
	                            double value)
	{
		return std::lower_bound(v.begin(), zero, value) - v.begin();
	}

	void updatePValues_(EnrichmentResults& tests, std::vector<double>& scores)
	{
		std::vector<double> mean_pos(tests.size());
		std::vector<double> mean_neg(tests.size());
		std::vector<double> normalized_scores(tests.size());

		size_t total_perm_pos = 0;
		size_t total_perm_neg = 0;
		size_t total_score_pos = 0;
		size_t total_score_neg = 0;
		for(size_t i = 0; i < tests.size(); ++i) {
			size_t num_pos = 0;
			size_t num_neg = 0;
			for(size_t j = 0; j < this->permutations_; ++j) {
				const auto& s = scores[i * this->permutations_ + j];

				if(s > 0) {
					++num_pos;
					mean_pos[i] += s;
				} else {
					++num_neg;
					mean_neg[i] -= s;
				}
			}

			total_perm_pos += num_pos;
			total_perm_neg += num_neg;

			mean_pos[i] /= num_pos;
			mean_neg[i] /= num_neg;

			std::cout << mean_pos[i] << ' ' << mean_neg[i] << '\n';

			for(size_t j = 0; j < this->permutations_; ++j) {
				auto& s = scores[i * this->permutations_ + j];

				s /= s > 0 ? mean_pos[i] : mean_neg[i];
			}

			// normalized_scores[i] = tests[i]->score / glob_mean;

			if(tests[i]->score > 0) {
				normalized_scores[i] = tests[i]->score / mean_pos[i];
				++total_score_pos;
			} else {
				normalized_scores[i] = tests[i]->score / mean_neg[i];
				++total_score_neg;
			}

			std::cout << tests[i]->category->name() << ' ' << tests[i]->score
			          << ' ' << normalized_scores[i] << '\n';
		}

		std::sort(scores.begin(), scores.end());
		std::sort(normalized_scores.begin(), normalized_scores.end());

		auto zero_p = std::lower_bound(scores.begin(), scores.end(), 0.0);
		auto zero_s = std::lower_bound(normalized_scores.begin(),
		                               normalized_scores.end(), 0.0);

		for(size_t i = 0; i < tests.size(); ++i) {
			size_t greater_p = 0;
			size_t greater_s = 0;
			if(tests[i]->score > 0) {
				greater_p = computeGreaterPos(scores, zero_p,
				                              tests[i]->score / mean_pos[i]);
				greater_s = computeGreaterPos(normalized_scores, zero_s,
				                              tests[i]->score / mean_pos[i]);
			} else {
				greater_p = computeGreaterNeg(scores, zero_p,
				                              tests[i]->score / mean_neg[i]);
				greater_s = computeGreaterNeg(normalized_scores, zero_s,
				                              tests[i]->score / mean_neg[i]);
			}

			if(greater_s == 0) {
				tests[i]->pvalue = 1.0;
			} else {
				if(tests[i]->score > 0) {
					tests[i]->pvalue =
					    static_cast<double>(greater_p * total_score_pos) /
					    (greater_s * total_perm_pos);
				} else {
					tests[i]->pvalue =
					    static_cast<double>(greater_p * total_score_neg) /
					    (greater_s * total_perm_neg);
				}
			}

			if(tests[i]->pvalue > 1.0) {
				tests[i]->pvalue = 1.0;
			}
		}
	}

	void performSinglePermutation_(size_t j,
	                               const EnrichmentAlgorithmPtr& algorithm,
	                               const EnrichmentResults& tests,
	                               std::vector<double>& permuted_values,
	                               std::vector<size_t>& column_indices)
	{
		Scores scores =
		    this->computeScores_(column_indices.begin(), column_indices.end());

		algorithm->setScores(scores);

		for(size_t i = 0; i < tests.size(); ++i) {
			auto score = std::get<0>(
			    algorithm->computeEnrichmentScore(*tests[i]->category));

			permuted_values[i * this->permutations_ + j] = score;
		}
	}

	void performSinglePermutationIndices_(
	    size_t j, const EnrichmentAlgorithmPtr& algorithm,
	    const EnrichmentResults& tests, std::vector<double>& permuted_values,
	    std::vector<size_t>& column_indices)
	{
		// Compute scores for a new permutation
		Scores scores =
		    this->computeScores_(column_indices.begin(), column_indices.end());

		// Update all the required lookup tables needed for finding
		// category members quickly.
		this->updateLookupTables_(tests, scores, algorithm->getOrder());

		// Now pass the scores to the algorithm
		algorithm->setScores(scores);

		for(size_t i = 0; i < tests.size(); ++i) {
			auto score = this->computeEnrichmentScore_(algorithm, tests, i);

			permuted_values[i * this->permutations_ + j] = score;
		}
	}
};
}

#endif // GT2_ENRICHMENT_PERMUTATION_TEST_H
