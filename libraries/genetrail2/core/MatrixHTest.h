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
#include "MatrixIterator.h"
#include "GeneSet.h"
#include "Scores.h"

#include "macros.h"

#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include <vector>
#include <utility>
#include <iterator>
#include <map>
#include <set>
#include <functional>
#include <tuple>
#include <initializer_list>
#include <cassert>
#include <cmath>

#include "compat.h"
#include "Exception.h"

namespace GeneTrail
{
	enum class NanMode { Ignore, Remove };

	struct Independent {
	};
	struct Dependent {
	};

	struct Scalar {
		using value_type = double;
	};
	struct Vectorized {
		using value_type = std::vector<double>;
	};

	struct IgnoreNaN {
	};
	struct RemoveNaN {
	};

	enum class MatrixHTests {
		IndependentTTest = 0,
		DependentTTest,
		IndependentShrinkageTTest,
		DependentShrinkageTTest,
		IndependentWilcoxonTest,
		DependentWilcoxonTest,
		FTest,
		SignalToNoiseRatio,
		LogMeanFoldQuotient,
		MeanFoldQuotient,
		MeanFoldDifference,
		PearsonCorrelation,
		SpearmanCorrelation,
		ZScore
	};

	class GT2_EXPORT MatrixHTestFactory
	{
		public:
		MatrixHTestFactory();

		template <class Iterator1, class Iterator2, class Vectorized = Scalar>
		class Test
		{
			public:
			virtual typename Vectorized::value_type
			    test(Iterator1, Iterator1, Iterator2, Iterator2) const = 0;
		};

		template <class Iterator1, class Iterator2>
		class Test<Iterator1, Iterator2, Vectorized>
		{
			public:
			virtual Vectorized::value_type test(Iterator1, Iterator1, Iterator2,
			                                    Iterator2) const = 0;
			virtual Vectorized::value_type testRemoveNaN(Iterator1, Iterator1,
			                                             Iterator2,
			                                             Iterator2) const = 0;
		};

		enum class IsDependent { Dependent, Independent };
		enum class SupportsVectors { Scalar, Vectorized };

		struct MethodDescriptor {
			MatrixHTests id;
			IsDependent is_dependent;
			SupportsVectors supports_vectors;
		};

		/**
		 * Factory method to obtain a
		 *
		 */
		template <class Iterator1, class Iterator2, class Dependent,
		          class Vectorized>
		std::unique_ptr<Test<Iterator1, Iterator2, Vectorized>>
		create(MatrixHTests method) const;

		template <class Iterator1, class Iterator2>
		std::unique_ptr<Test<Iterator1, Iterator2, Vectorized>>
		create(MatrixHTests method, Dependent, Vectorized) const
		{
			switch(method) {
				case MatrixHTests::DependentShrinkageTTest:
					throw NotImplemented(__FILE__, __LINE__,
					                     "Dependent Shrinkage t-Test");
				default:
					throw std::invalid_argument("Reqested method does not "
					                            "support the proper "
					                            "interface.");
			}
		}

		template <class Iterator1, class Iterator2>
		std::unique_ptr<Test<Iterator1, Iterator2, Vectorized>>
		create(MatrixHTests method, Independent, Vectorized) const
		{
			switch(method) {
				case MatrixHTests::IndependentShrinkageTTest:
					return std::make_unique<
					    independent_shrinkage_t_test<Iterator1, Iterator2>>();
				default:
					throw std::invalid_argument("Requested method does not "
					                            "support the proper "
					                            "interface.");
			}
		}

		template <class Iterator1, class Iterator2>
		std::unique_ptr<Test<Iterator1, Iterator2>>
		create(MatrixHTests method, Dependent, Scalar) const
		{
			switch(method) {
				case MatrixHTests::DependentTTest:
					return std::make_unique<
					    dependent_t_test<Iterator1, Iterator2>>();
				case MatrixHTests::DependentWilcoxonTest:
					return std::make_unique<
					    dependent_wilcoxon<Iterator1, Iterator2>>();
				default:
					throw std::invalid_argument("Reqested method does not "
					                            "support the proper "
					                            "interface.");
			}
		}

		template <class Iterator1, class Iterator2>
		std::unique_ptr<Test<Iterator1, Iterator2>>
		create(MatrixHTests method, Independent, Scalar) const
		{
			switch(method) {
				case MatrixHTests::IndependentTTest:
					return std::make_unique<
					    independent_t_test<Iterator1, Iterator2>>();
				case MatrixHTests::IndependentWilcoxonTest:
					return std::make_unique<
					    independent_wilcoxon<Iterator1, Iterator2>>();
				case MatrixHTests::FTest:
					return std::make_unique<f_test<Iterator1, Iterator2>>();
				case MatrixHTests::SignalToNoiseRatio:
					return std::make_unique<
					    signal_to_noise_ratio<Iterator1, Iterator2>>();
				case MatrixHTests::LogMeanFoldQuotient:
					return std::make_unique<
					    log_mean_fold_quotient<Iterator1, Iterator2>>();
				case MatrixHTests::MeanFoldQuotient:
					return std::make_unique<
					    mean_fold_quotient<Iterator1, Iterator2>>();
				case MatrixHTests::MeanFoldDifference:
					return std::make_unique<
					    mean_fold_difference<Iterator1, Iterator2>>();
				case MatrixHTests::PearsonCorrelation:
					return std::make_unique<
					    pearson_correlation<Iterator1, Iterator2>>();
				case MatrixHTests::SpearmanCorrelation:
					return std::make_unique<
					    spearman_correlation<Iterator1, Iterator2>>();
				case MatrixHTests::ZScore:
					return std::make_unique<z_score<Iterator1, Iterator2>>();
				default:
					throw std::invalid_argument("Reqested method does not "
					                            "support the proper "
					                            "interface.");
			}
		}

		boost::optional<MatrixHTests> getMethod(const std::string& method) const
		{
			auto res = ids_.find(method);

			if(res == ids_.end()) {
				return boost::none;
			}

			return boost::make_optional(res->second);
		};

		MethodDescriptor getDescriptor(const std::string& method) const
		{
			auto res = ids_.find(method);

			if(res == ids_.end()) {
				throw std::invalid_argument(method);
			}

			return getDescriptor(res->second);
		};

		MethodDescriptor getDescriptor(MatrixHTests method) const
		{
			return descriptors_[static_cast<size_t>(method)];
		};

		bool hasTest(const std::string& method) const
		{
			return ids_.find(method) != ids_.end();
		}

		private:
		template <class Iterator1, class Iterator2>
		class independent_t_test : public Test<Iterator1, Iterator2>
		{
			double test(Iterator1 fst_begin, Iterator1 fst_end,
			            Iterator2 snd_begin, Iterator2 snd_end) const override
			{
				IndependentTTest<double> test;
				return HTest::test(test, fst_begin, fst_end, snd_begin,
				                   snd_end);
			}
		};

		template <class Iterator1, class Iterator2>
		class dependent_t_test : public Test<Iterator1, Iterator2>
		{
			double test(Iterator1 fst_begin, Iterator1 fst_end,
			            Iterator2 snd_begin, Iterator2 snd_end) const override
			{
				DependentTTest<double> test;
				return HTest::test(test, fst_begin, fst_end, snd_begin,
				                   snd_end);
			}
		};

		template <class Iterator1, class Iterator2>
		class independent_shrinkage_t_test
		    : public Test<Iterator1, Iterator2, Vectorized>
		{
			Vectorized::value_type test(Iterator1 fst_begin, Iterator1 fst_end,
			                            Iterator2 snd_begin,
			                            Iterator2 snd_end) const override
			{
				IndependentShrinkageTTest<double> test;
				return test.test(fst_begin, fst_end, snd_begin, snd_end);
			}

			Vectorized::value_type
			testRemoveNaN(Iterator1 fst_begin, Iterator1 fst_end,
			              Iterator2 snd_begin, Iterator2 snd_end) const override
			{
				IndependentShrinkageTTest<double> test;
				return test.testRemoveNaN(fst_begin, fst_end, snd_begin,
				                          snd_end);
			}
		};

		template <class Iterator1, class Iterator2>
		class f_test : public Test<Iterator1, Iterator2>
		{
			double test(Iterator1 fst_begin, Iterator1 fst_end,
			            Iterator2 snd_begin, Iterator2 snd_end) const override
			{
				FTest<double> test;
				return HTest::test(test, fst_begin, fst_end, snd_begin,
				                   snd_end);
			}
		};

		template <class Iterator1, class Iterator2>
		class independent_wilcoxon : public Test<Iterator1, Iterator2>
		{
			double test(Iterator1 fst_begin, Iterator1 fst_end,
			            Iterator2 snd_begin, Iterator2 snd_end) const override
			{
				WilcoxonRankSumTest<double> test;
				return HTest::test(test, fst_begin, fst_end, snd_begin,
				                   snd_end);
			}
		};

		template <class Iterator1, class Iterator2>
		class dependent_wilcoxon : public Test<Iterator1, Iterator2>
		{
			double test(Iterator1 fst_begin, Iterator1 fst_end,
			            Iterator2 snd_begin, Iterator2 snd_end) const override
			{
				WilcoxonMatchedPairsSignedRankTest<double> test;
				return HTest::test(test, fst_begin, fst_end, snd_begin,
				                   snd_end);
			}
		};

		template <class Iterator1, class Iterator2>
		class signal_to_noise_ratio : public Test<Iterator1, Iterator2>
		{
			double test(Iterator1 fst_begin, Iterator1 fst_end,
			            Iterator2 snd_begin, Iterator2 snd_end) const override
			{
				SignalToNoiseRatio<double> test;
				return HTest::test(test, fst_begin, fst_end, snd_begin,
				                   snd_end);
			}
		};

		template <class Iterator1, class Iterator2>
		class log_mean_fold_quotient : public Test<Iterator1, Iterator2>
		{
			double test(Iterator1 fst_begin, Iterator1 fst_end,
			            Iterator2 snd_begin, Iterator2 snd_end) const override
			{
				return statistic::log_mean_fold_quotient<double>(
				    fst_begin, fst_end, snd_begin, snd_end);
			}
		};

		template <class Iterator1, class Iterator2>
		class mean_fold_quotient : public Test<Iterator1, Iterator2>
		{
			double test(Iterator1 fst_begin, Iterator1 fst_end,
			            Iterator2 snd_begin, Iterator2 snd_end) const override
			{
				return statistic::mean_fold_quotient<double>(
				    fst_begin, fst_end, snd_begin, snd_end);
			}
		};

		template <class Iterator1, class Iterator2>
		class mean_fold_difference : public Test<Iterator1, Iterator2>
		{
			double test(Iterator1 fst_begin, Iterator1 fst_end,
			            Iterator2 snd_begin, Iterator2 snd_end) const override
			{
				return statistic::mean_fold_difference<double>(
				    fst_begin, fst_end, snd_begin, snd_end);
			}
		};

		template <class Iterator1, class Iterator2>
		class pearson_correlation : public Test<Iterator1, Iterator2>
		{
			double test(Iterator1 fst_begin, Iterator1 fst_end,
			            Iterator2 snd_begin, Iterator2 snd_end) const override
			{
				const auto n = std::distance(fst_begin, fst_end);

				std::vector<double> a(fst_begin, fst_end);
				a.insert(a.end(), snd_begin, snd_end);
				std::vector<double> b(a.size(), 0.0);
				std::fill(b.begin() + n, a.end(), 1.0);

				return statistic::pearson_correlation<double>(
				    a.begin(), a.end(), b.begin(), b.end());
			}
		};

		template <class Iterator1, class Iterator2>
		class spearman_correlation : public Test<Iterator1, Iterator2>
		{
			double test(Iterator1 fst_begin, Iterator1 fst_end,
			            Iterator2 snd_begin, Iterator2 snd_end) const override
			{
				const auto n = std::distance(fst_begin, fst_end);

				std::vector<double> a(fst_begin, fst_end);
				a.insert(a.end(), snd_begin, snd_end);
				std::vector<double> b(a.size(), 0.0);
				std::fill(b.begin() + n, a.end(), 1.0);

				return statistic::spearman_correlation<double>(
				    a.begin(), a.end(), b.begin(), b.end());
			}
		};

		template <class Iterator1, class Iterator2>
		class z_score : public Test<Iterator1, Iterator2>
		{
			double test(Iterator1 fst_begin, Iterator1 fst_end,
			            Iterator2 snd_begin, Iterator2 snd_end) const override
			{
				(void)fst_end;
				assert(fst_begin != fst_end);
				double x = *fst_begin;
				return statistic::z_score<double>(x, snd_begin, snd_end);
			}
		};

		std::unordered_map<std::string, MatrixHTests> ids_;
		std::vector<MethodDescriptor> descriptors_;

		bool descriptorsProperlyInitialized_() const;
	};

	class MatrixHTest
	{
		public:
		void setRowDBIndices(const std::vector<size_t>& row_db_indices) {
			row_db_indices_ = row_db_indices;
		}

		template <typename Matrix>
		Scores test(const std::string& method, const Matrix& ref,
		            const Matrix& sam, NanMode mode = NanMode::Ignore)
		{
			return test_(factory.getDescriptor(method), ref, sam, mode);
		}

		template <typename Matrix>
		Scores test(MatrixHTests method, const Matrix& ref,
		            const Matrix& sam, NanMode mode = NanMode::Ignore)
		{
			return test_(factory.getDescriptor(method), ref, sam, mode);
		}

		private:
		MatrixHTestFactory factory;
		std::vector<size_t> row_db_indices_;

		template <typename Matrix>
		Scores test_(const MatrixHTestFactory::MethodDescriptor& descriptor,
		             const Matrix& ref, const Matrix& sam, NanMode mode) const
		{
			assert(ref.rows() == sam.rows());
			assert(ref.rowNames() == sam.rowNames());

			switch(mode) {
				case NanMode::Ignore:
					return test_(descriptor, ref, sam, IgnoreNaN());
				case NanMode::Remove:
					return test_(descriptor, ref, sam, RemoveNaN());
				default:
					throw NotImplemented(__FILE__, __LINE__, "The passed NaNMode is not yet implemented");
			};
		}

		template <typename Matrix, typename NaNStrategy>
		Scores test_(const MatrixHTestFactory::MethodDescriptor& descriptor,
		             const Matrix& ref, const Matrix& sam, NaNStrategy) const
		{
			switch(descriptor.is_dependent) {
				case MatrixHTestFactory::IsDependent::Dependent:
					return test_(descriptor, ref, sam, NaNStrategy(),
					             Dependent());
				case MatrixHTestFactory::IsDependent::Independent:
					return test_(descriptor, ref, sam, NaNStrategy(),
					             Independent());
			}

			throw std::logic_error("Unreachable code");
		}

		template <typename Matrix, typename NaNStrategy, typename Dependency>
		Scores test_(const MatrixHTestFactory::MethodDescriptor& descriptor,
		             const Matrix& ref, const Matrix& sam, NaNStrategy,
		             Dependency) const
		{
			switch(descriptor.supports_vectors) {
				case MatrixHTestFactory::SupportsVectors::Scalar:
					return test_(descriptor, ref, sam, NaNStrategy(),
					             Dependency(), Scalar());
				case MatrixHTestFactory::SupportsVectors::Vectorized:
					return test_(descriptor, ref, sam, NaNStrategy(),
					             Dependency(), Vectorized());
			}

			throw std::logic_error("Unreachable code");
		}

		template <typename Matrix, typename NaNStrategy, typename Dependency,
		          typename Scalar>
		Scores test_(const MatrixHTestFactory::MethodDescriptor& descriptor,
		             const Matrix& ref, const Matrix& sam, NaNStrategy,
		             Dependency, Scalar) const;

		template <typename Matrix>
		Scores test_(const MatrixHTestFactory::MethodDescriptor& descriptor,
		             const Matrix& ref, const Matrix& sam, RemoveNaN, Dependent,
		             Scalar) const
		{
			Scores scores(ref.rows());

			RowMajorMatrixIterator<Matrix> ref_it(&ref, 0), sam_it(&sam, 0);

			using value_type = typename Matrix::value_type;
			using Paired = boost::tuple<value_type, value_type>;

			// Predicate that is true if a floating point value
			// is not NaN
			auto both_not_nan = [](const Paired& x) {
				return !std::isnan(boost::get<0>(x)) &&
				       !std::isnan(boost::get<1>(x));
			};

			for(const auto& name : ref.rowNames()) {
				auto zip_begin = boost::make_zip_iterator(
				    boost::make_tuple(ref_it->begin(), sam_it->begin()));
				auto zip_end = boost::make_zip_iterator(
				    boost::make_tuple(ref_it->end(), sam_it->end()));

				auto filter_zip_begin = boost::make_filter_iterator(
				    both_not_nan, zip_begin, zip_end);
				auto filter_zip_end =
				    boost::make_filter_iterator(both_not_nan, zip_end, zip_end);

				auto fst = [](const Paired& x) { return boost::get<0>(x); };
				auto snd = [](const Paired& x) { return boost::get<1>(x); };

				auto ref_nan_begin =
				    boost::make_transform_iterator(filter_zip_begin, fst);
				auto ref_nan_end =
				    boost::make_transform_iterator(filter_zip_end, fst);
				auto sam_nan_begin =
				    boost::make_transform_iterator(filter_zip_begin, snd);
				auto sam_nan_end =
				    boost::make_transform_iterator(filter_zip_end, snd);

				using Iterator1 = decltype(ref_nan_begin);
				using Iterator2 = decltype(sam_nan_begin);

				auto method = factory.create<Iterator1, Iterator2>(
				    descriptor.id, Independent(), Scalar());

				auto score = method->test(ref_nan_begin, ref_nan_end,
				                          sam_nan_begin, sam_nan_end);

				scores.emplace_back(name, score);

				++ref_it;
				++sam_it;
			}

			return scores;
		}

		template <typename Matrix>
		Scores test_(const MatrixHTestFactory::MethodDescriptor& descriptor,
		             const Matrix& ref, const Matrix& sam, RemoveNaN,
		             Independent, Scalar) const
		{
			Scores scores(ref.rows());

			RowMajorMatrixIterator<Matrix> ref_it(&ref, 0), sam_it(&sam, 0);

			// Predicate that is true if a floating point value
			// is not NaN
			auto is_not_nan =
			    [](typename Matrix::value_type x) { return !std::isnan(x); };

			for(const auto& name : ref.rowNames()) {
				auto ref_nan_begin = boost::make_filter_iterator(
				    is_not_nan, ref_it->begin(), ref_it->end());
				auto ref_nan_end = boost::make_filter_iterator(
				    is_not_nan, ref_it->end(), ref_it->end());
				auto sam_nan_begin = boost::make_filter_iterator(
				    is_not_nan, sam_it->begin(), sam_it->end());
				auto sam_nan_end = boost::make_filter_iterator(
				    is_not_nan, sam_it->end(), sam_it->end());

				using Iterator1 = decltype(ref_nan_begin);
				using Iterator2 = decltype(sam_nan_begin);

				auto method = factory.create<Iterator1, Iterator2>(
				    descriptor.id, Independent(), Scalar());

				auto score = method->test(ref_nan_begin, ref_nan_end,
				                          sam_nan_begin, sam_nan_end);
				scores.emplace_back(name, score);

				++ref_it;
				++sam_it;
			}

			return scores;
		}

		template <typename Matrix, typename Dep>
		Scores test_(const MatrixHTestFactory::MethodDescriptor& descriptor,
		             const Matrix& ref, const Matrix& sam, RemoveNaN,
		             Dep, Vectorized) const
		{
			Scores scores(ref.rows());

			RowMajorMatrixIterator<Matrix> ref_begin(&ref, 0),
			    ref_end(&ref, ref.rows()), sam_begin(&sam, 0),
			    sam_end(&sam, sam.rows());

			using Iterator1 = decltype(ref_begin);
			using Iterator2 = decltype(sam_begin);

			auto method = factory.create<Iterator1, Iterator2>(
			    descriptor.id, Dep(), Vectorized());

			auto v =
			    method->testRemoveNaN(ref_begin, ref_end, sam_begin, sam_end);

			assignScores_(scores, v, ref);

			return scores;
		}

		template <typename Matrix, typename Dep>
		Scores test_(const MatrixHTestFactory::MethodDescriptor& descriptor,
		             const Matrix& ref, const Matrix& sam, IgnoreNaN, Dep,
		             Scalar) const
		{
			Scores scores(ref.rows());

			RowMajorMatrixIterator<Matrix> ref_it(&ref, 0), sam_it(&sam, 0);

			using Iterator1 = decltype(ref_it->begin());
			using Iterator2 = decltype(sam_it->begin());

			auto method = factory.create<Iterator1, Iterator2>(descriptor.id,
			                                                   Dep(), Scalar());

			for(size_t r = 0; r < ref.rows(); ++r, ++ref_it, ++sam_it) {
				auto score = method->test(ref_it->begin(), ref_it->end(),
				                          sam_it->begin(), sam_it->end());
				if(row_db_indices_.empty()) {
					scores.emplace_back(ref.rowName(r), score);
				} else {
					scores.emplace_back(row_db_indices_[r], score);
				}
			}

			return scores;
		}

		template <typename Matrix, typename Dep>
		Scores test_(const MatrixHTestFactory::MethodDescriptor& descriptor,
		             const Matrix& ref, const Matrix& sam, IgnoreNaN, Dep,
		             Vectorized) const
		{
			Scores scores(ref.rows());

			RowMajorMatrixIterator<Matrix> ref_begin(&ref, 0),
			    ref_end(&ref, ref.rows()), sam_begin(&sam, 0),
			    sam_end(&sam, sam.rows());

			using Iterator1 = decltype(ref_begin);
			using Iterator2 = decltype(sam_begin);

			auto method = factory.create<Iterator1, Iterator2>(
			    descriptor.id, Dep(), Vectorized());

			auto v = method->test(ref_begin, ref_end, sam_begin, sam_end);

			assignScores_(scores, v, ref);

			return scores;
		}

		template<typename Matrix>
		void assignScores_(Scores& scores, const std::vector<typename Matrix::value_type>& v, const Matrix& ref) const {
			if(row_db_indices_.empty()) {
				for(unsigned int r = 0; r < ref.rows(); ++r) {
					scores.emplace_back(ref.rowName(r), v[r]);
				}
			} else {
				for(unsigned int r = 0; r < ref.rows(); ++r) {
					scores.emplace_back(row_db_indices_[r], v[r]);
				}
			}
		}

	};
}

#endif // GT2_CORE_MATRIX_HTEST_H
