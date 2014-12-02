#ifndef GT2_CORE_MATRIX_HTEST_H
#define GT2_CORE_MATRIX_HTEST_H

#include "HTest.h"
#include "FTest.h"
#include "IndependentTTest.h"
#include "DependentTTest.h"
#include "WilcoxonRankSumTest.h"
#include "WilcoxonMatchedPairsSignedRankTest.h"
#include "IndependentShrinkageTTest.h"
#include "SignalToNoiseRatio.h"
#include "IndexIterator.h"
#include "GeneSet.h"

#include "macros.h"

#include <vector>
#include <utility>
#include <iterator>
#include <map>
#include <set>
#include <functional>
#include <tuple>
#include <initializer_list>
#include <cassert>

namespace GeneTrail
{
	template <class Matrix> class GT2_EXPORT MatrixHTest
	{
		public:
		MatrixHTest() {};

		private:
		typedef std::vector<double> _v;
		typedef std::vector<std::vector<double>> _vv;
		typedef std::vector<double>::iterator _viter;
		typedef std::function<double(_v, _v)> return_type;
		typedef std::function<GeneSet(Matrix&, Matrix&)> matrix_func;

		return_type independent_t_test = [](_v a, _v b) {
			IndependentTTest<double> independent_t;
			return HTest::test(independent_t, a.begin(), a.end(), b.begin(),
			                   b.end());
		};

		return_type dependent_t_test = [](_v a, _v b) {
			DependentTTest<double> dependent_t;
			return HTest::test(dependent_t, a.begin(), a.end(), b.begin(),
			                   b.end());
		};

		return_type f_test = [](_v a, _v b) {
			FTest<double> f;
			return HTest::test(f, a.begin(), a.end(), b.begin(), b.end());
		};

		return_type wilcoxon = [](_v a, _v b) {
			WilcoxonRankSumTest<double> t;
			return HTest::test(t, a.begin(), a.end(), b.begin(), b.end());
		};

		return_type dependent_wilcoxon = [](_v a, _v b) {
			WilcoxonMatchedPairsSignedRankTest<double> t;
			return HTest::test(t, a.begin(), a.end(), b.begin(), b.end());
		};

		return_type signal_to_noise_ratio = [](_v a, _v b) {
			SignalToNoiseRatio<double> t;
			return HTest::test(t, a.begin(), a.end(), b.begin(), b.end());
		};

		return_type log_mean_fold_quotient = [](_v a, _v b) {
			return statistic::log_mean_fold_quotient<double>(
			    a.begin(), a.end(), b.begin(), b.end());
		};

		return_type mean_fold_quotient = [](_v a, _v b) {
			return statistic::mean_fold_quotient<double>(a.begin(), a.end(),
			                                             b.begin(), b.end());
		};

		return_type mean_fold_difference = [](_v a, _v b) {
			return statistic::mean_fold_difference<double>(a.begin(), a.end(),
			                                               b.begin(), b.end());
		};

		return_type pearson_correlation = [](_v x, _v y) {
			std::vector<double> a(x.begin(), x.end());
			a.insert(a.end(), y.begin(), y.end());
			std::vector<double> b(a.size(), 0.0);
			std::fill(b.begin() + x.size(), a.end(), 1.0);

			return statistic::pearson_correlation<double>(a.begin(), a.end(),
			                                              b.begin(), b.end());
		};

		return_type spearman_correlation = [](_v x, _v y) {
			std::vector<double> a(x.begin(), x.end());
			a.insert(a.end(), y.begin(), y.end());
			std::vector<double> b(a.size(), 0.0);
			std::fill(b.begin() + x.size(), a.end(), 1.0);

			return statistic::spearman_correlation<double>(a.begin(), a.end(),
			                                               b.begin(), b.end());
		};

		return_type z_score = [](_v a, _v b) {
			assert(a.size() == 1);
			double x = a[0];
			return statistic::z_score<double>(x, b.begin(), b.end());
		};

		std::map<std::string, return_type> rowWiseMethods{
		    {"independent-t-test", independent_t_test},
		    {"dependent-t-test", dependent_t_test},
		    {"f-test", f_test},
		    {"wilcoxon", wilcoxon},
		    {"dependent-wilcoxon", dependent_wilcoxon},
		    {"signal-to-noise-ratio", signal_to_noise_ratio},
		    {"log-mean-fold-quotient", log_mean_fold_quotient},
		    {"mean-fold-quotient", mean_fold_quotient},
		    {"mean-fold-difference", mean_fold_difference},
		    {"pearson_correlation", pearson_correlation},
		    {"spearman_correlation", spearman_correlation},
		    {"z-score", z_score}};

		std::set<std::string> dependentTests{"dependent-t-test",
		                                     "dependent-wilcoxon"};

		/**
		 * This functions removes all NANs from the given row and returns a
		 *vector containing all values.
		 *
		 * @param m Matrix
		 * @param r Row index
		 * @return Vector containing all values that are not NAN.
		 */
		std::vector<double> removeNANs(Matrix& m, unsigned int r)
		{
			std::vector<double> result;
			result.reserve(m.cols());
			for(unsigned int c = 0; c < m.cols(); ++c) {
				// This checks if m(r,c) is NAN, NA or Inf...
				if(m(r, c) != m(r, c)) {
					continue;
				}
				result.push_back(m(r, c));
			}
			return result;
		}

		/**
		 * This functions removes all NANs from the given matrix and returns a
		 * vector containing all values.
		 * @param m Matrix
		 * @return Vector containing all values that are not NAN.
		 */
		std::vector<std::vector<double>> removeAllNANs(Matrix& m)
		{
			std::vector<std::vector<double>> result;
			result.reserve(m.rows());
			for(unsigned int r = 0; r < m.rows(); ++r) {
				result.emplace_back(removeNANs(m, r));
			}
			return result;
		}

		/**
		 * This functions removes all NANs from the given ranges and returns two
		 *vectors.
		 *
		 * @param r Matrix
		 * @param s Matrix
		 * @param r Row index
		 * @return Pair of Vectors containing all values that are not NAN.
		 */
		std::tuple<std::vector<double>, std::vector<double>>
		removePairwiseNANs(Matrix& ref, Matrix& s, unsigned int r)
		{
			std::vector<double> result1;
			result1.reserve(ref.cols());
			std::vector<double> result2;
			result2.reserve(s.cols());
			for(unsigned int c = 0; c < ref.cols(); ++c) {
				// This checks if m(r,c) is NAN, NA or Inf...
				if(ref(r, c) != ref(r, c) || s(r, c) != s(r, c)) {
					continue;
				}
				result1.push_back(ref(r, c));
				result2.push_back(s(r, c));
			}
			return make_tuple(result1, result2);
		}

		/**
		 * This functions removes all NANs from the given matrix and returns a
		 * vector containing all values.
		 * @param m Matrix
		 * @return Vector containing all values that are not NAN.
		 */
		std::tuple<std::vector<std::vector<double>>,
		           std::vector<std::vector<double>>>
		removeAllPairwiseNANs(Matrix& reference, Matrix& sample)
		{
			std::vector<std::vector<double>> result1;
			result1.reserve(reference.rows());
			std::vector<std::vector<double>> result2;
			result2.reserve(sample.rows());
			for(unsigned int r = 0; r < reference.rows(); ++r) {
				auto tmp = removePairwiseNANs(reference, sample, r);
				result1.emplace_back(std::get<0>(tmp));
				result2.emplace_back(std::get<1>(tmp));
			}
			return std::make_tuple(result1, result2);
		}

		/**
		 * This functions applies thegiven test to the two given vectors.
		 * The vectors represent the different groups.
		 *
		 * @param matrix A matrix object
		 * @param test The test to apply
		 * @returns the score for each row
		 */
		template <typename Test>
		double apply(std::vector<double> reference, std::vector<double> sample,
		             Test test)
		{
			return test(reference, sample);
		}

		/**
		 * This functions applies the given test row-wise to the given matrices.
		 * The matrices represent the different groups.
		 *
		 * @param reference A matrix object
		 * @param sample A matrix object
		 * @param test The test to apply
		 * @returns the score for each row
		 */
		GeneSet independentRApply(Matrix& reference, Matrix& sample,
		                          std::string test)
		{
			GeneSet result;
			auto genes = reference.rowNames();
			for(unsigned int r = 0; r < reference.rows(); ++r) {
				std::vector<double> v1 = removeNANs(reference, r);
				std::vector<double> v2 = removeNANs(sample, r);
				result.insert(genes[r], apply(v1, v2, rowWiseMethods[test]));
			}
			return result;
		}

		/**
		 * This functions applies the given test row-wise to the given matrices.
		 * The matrices represent the different groups.
		 *
		 * @param reference A matrix object
		 * @param sample A matrix object
		 * @param test The test to apply
		 * @returns the score for each row
		 */
		GeneSet dependentRApply(Matrix& reference, Matrix& sample,
		                        std::string test)
		{
			GeneSet result;
			auto genes = reference.rowNames();
			for(unsigned int r = 0; r < reference.rows(); ++r) {
				std::vector<double> v1;
				std::vector<double> v2;
				std::tie(v1, v2) = removePairwiseNANs(reference, sample, r);
				result.insert(genes[r], apply(v1, v2, rowWiseMethods[test]));
			}
			return result;
		}

		/**
		 *
		 */
		GeneSet independentShrinkageTTest(Matrix& reference, Matrix& sample)
		{
			GeneSet result;
			IndependentShrinkageTTest<double> t;
			auto ref = removeAllNANs(reference);
			auto sam = removeAllNANs(sample);
			auto v = t.test(ref.begin(), ref.end(), sam.begin(), sam.end());
			for(unsigned int r = 0; r < reference.rows(); ++r) {
				result.insert(reference.rowName(r), v[r]);
			}
			return result;
		}

		public:
		/**
		 * This functions applies the given test row-wise to the given matrices.
		 * The matrices represent the different groups.
		 *
		 * @param reference A matrix object
		 * @param sample A matrix object
		 * @param test The test to apply
		 * @returns A GeneSet object containing a score for each row
		 */
		GeneSet test(Matrix& reference, Matrix& sample, const std::string& test)
		{
			if(test == "independent-shrinkage-t-test") {
				return independentShrinkageTTest(reference, sample);
			} else {
				if(dependentTests.find(test) != dependentTests.end()) {
					assert(reference.cols() == sample.cols());
					return dependentRApply(reference, sample, test);
				}
				return independentRApply(reference, sample, test);
			}
		}
	};
}

#endif // GT2_CORE_MATRIX_HTEST_H

