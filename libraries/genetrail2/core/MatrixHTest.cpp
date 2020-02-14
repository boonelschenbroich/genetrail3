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

		MethodDescriptor mean_first_group_descriptor{
			MatrixHTests::MeanFirstGroup,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor mean_larger_zero_descriptor{
			MatrixHTests::MeanLargerZero,
			IsDependent::Independent,
			SupportsVectors::Scalar};

		MethodDescriptor larger_zero_descriptor{
			MatrixHTests::LargerZero,
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
		ids_.emplace("mean-first-group", MatrixHTests::MeanFirstGroup);
		ids_.emplace("mean-larger-zero", MatrixHTests::MeanLargerZero);
		ids_.emplace("larger-zero", MatrixHTests::LargerZero);

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
			z_score_descriptor,
			mean_first_group_descriptor,
			mean_larger_zero_descriptor,
			larger_zero_descriptor
		};

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
