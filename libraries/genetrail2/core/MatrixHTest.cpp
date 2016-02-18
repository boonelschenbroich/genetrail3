#include "MatrixHTest.h"

namespace GeneTrail
{
	MatrixHTestFactory::MatrixHTestFactory()
	{
		MethodDescriptor independent_t_test_descriptor{
			MatrixHTests::IndependentTTest,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor dependent_t_test_descriptor{
			MatrixHTests::DependentTTest,
			IsDependent::Dependent,
			SupportsVectors::Scalar};

		MethodDescriptor independent_shrinkage_t_test_descriptor{
			MatrixHTests::IndependentShrinkageTTest,
			IsDependent::Independent,
			SupportsVectors::Vectorized};

		// This is currently a dummy and might need to be adjusted
		// in the future.
		MethodDescriptor dependent_shrinkage_t_test_descriptor{
			MatrixHTests::DependentShrinkageTTest,
			IsDependent::Dependent,
			SupportsVectors::Vectorized};

		MethodDescriptor f_test_descriptor{
			MatrixHTests::FTest,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor independent_wilcoxon_descriptor{
			MatrixHTests::IndependentWilcoxonTest,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor dependent_wilcoxon_descriptor{
			MatrixHTests::DependentWilcoxonTest,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor signal_to_noise_ratio_descriptor{
			MatrixHTests::SignalToNoiseRatio,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor log_mean_fold_quotient_descriptor{
			MatrixHTests::LogMeanFoldQuotient,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor log_median_fold_quotient_descriptor{
			MatrixHTests::LogMedianFoldQuotient,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor mean_fold_quotient_descriptor{
			MatrixHTests::MeanFoldQuotient,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor mean_fold_difference_descriptor{
			MatrixHTests::MeanFoldDifference,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor median_fold_quotient_descriptor{
			MatrixHTests::MedianFoldQuotient,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor median_fold_difference_descriptor{
			MatrixHTests::MedianFoldDifference,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor pearson_correlation_descriptor{
			MatrixHTests::PearsonCorrelation,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor spearman_correlation_descriptor{
			MatrixHTests::SpearmanCorrelation,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor z_score_descriptor{
			MatrixHTests::ZScore,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		ids_.emplace("independent-t-test", MatrixHTests::IndependentTTest);
		ids_.emplace("dependent-t-test", MatrixHTests::DependentTTest);
		ids_.emplace("independent-shrinkage-t-test", MatrixHTests::IndependentShrinkageTTest);
		ids_.emplace("dependent-shrinkage-t-test", MatrixHTests::DependentShrinkageTTest);
		ids_.emplace("f-test", MatrixHTests::FTest);
		ids_.emplace("wilcoxon", MatrixHTests::IndependentWilcoxonTest);
		ids_.emplace("dependent-wilcoxon", MatrixHTests::DependentWilcoxonTest);
		ids_.emplace("signal-to-noise-ratio", MatrixHTests::SignalToNoiseRatio);
		ids_.emplace("log-mean-fold-quotient", MatrixHTests::LogMeanFoldQuotient);
		ids_.emplace("log-median-fold-quotient", MatrixHTests::LogMedianFoldQuotient);
		ids_.emplace("mean-fold-quotient", MatrixHTests::MeanFoldQuotient);
		ids_.emplace("mean-fold-difference", MatrixHTests::MeanFoldDifference);
		ids_.emplace("median-fold-quotient", MatrixHTests::MedianFoldQuotient);
		ids_.emplace("median-fold-difference", MatrixHTests::MedianFoldDifference);
		ids_.emplace("pearson_correlation", MatrixHTests::PearsonCorrelation);
		ids_.emplace("spearman_correlation", MatrixHTests::SpearmanCorrelation);
		ids_.emplace("z-score", MatrixHTests::ZScore);

		descriptors_ = std::vector<MethodDescriptor>{
			independent_t_test_descriptor,
			dependent_t_test_descriptor,
			independent_shrinkage_t_test_descriptor,
			dependent_shrinkage_t_test_descriptor,
			independent_wilcoxon_descriptor,
			dependent_wilcoxon_descriptor,
			f_test_descriptor,
			signal_to_noise_ratio_descriptor,
			log_mean_fold_quotient_descriptor,
			log_median_fold_quotient_descriptor,
			mean_fold_quotient_descriptor,
			mean_fold_difference_descriptor,
			median_fold_quotient_descriptor,
			median_fold_difference_descriptor,
			pearson_correlation_descriptor,
			spearman_correlation_descriptor,
			z_score_descriptor};

		assert(descriptorsProperlyInitialized_());
	}

	bool MatrixHTestFactory::descriptorsProperlyInitialized_() const
	{
		for(size_t i = 0; i < descriptors_.size(); ++i) {
			if(static_cast<size_t>(descriptors_[i].id) != i) {
				return false;
			}
		}

		return true;
	}
}
