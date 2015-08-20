#include "OverRepresentationAnalysis.h"

using namespace GeneTrail;

OverRepresentationAnalysis::OverRepresentationAnalysis(
    const Category& reference_set, const Category& test_set)
    : reference_set_(reference_set), test_set_(test_set)
{
	m_ = reference_set_.size();
	n_ = test_set_.size();
	useHypergeometricTest_ =
	    categoryContainsAllGenes(reference_set_, test_set_);
}

bool
OverRepresentationAnalysis::categoryContainsAllGenes(const Category& reference,
                                                     const Category& testSet)
{
	for(const auto& gene : testSet) {
		if(!reference.contains(gene)) {
			return false;
		}
	}

	return true;
}

double OverRepresentationAnalysis::expectedNumberOfHits(const Category& category) const {
	auto l = Category::intersect("null", category, reference_set_).size();
	return (l * n_) / static_cast<double>(m_);
}

double OverRepresentationAnalysis::computePValue(const Category& category) const
{
	// GeneTrail 1
	// size_t l = category.size();

	size_t k = Category::intersect("null", category, test_set_).size();
	size_t l = Category::intersect("null", category, reference_set_).size();

	auto expected_k = ((double)l * n_) / ((double)m_);
	bool enriched = expected_k < k;

	big_float p;

	if(useHypergeometricTest_) {
		if(enriched) {
			p = hyperTest_.upperTailedPValue(m_, l, n_, k);
		} else {
			p = hyperTest_.lowerTailedPValue(m_, l, n_, k);
		}
	} else {
		if(enriched) {
			p = fisherTest_.upperTailedPValue(m_, l, n_, k);
		} else {
			p = fisherTest_.lowerTailedPValue(m_, l, n_, k);
		}
	}

	return p.convert_to<double>();
}

double OverRepresentationAnalysis::computeScore(const Category& category) const
{
	// GeneTrail 1
	// size_t l = category.size();

	size_t k = Category::intersect("null", category, test_set_).size();
	size_t l = Category::intersect("null", category, reference_set_).size();

	if(useHypergeometricTest_) {
		return hyperTest_.compute(m_, l, n_, k).convert_to<double>();
	} else {
		return fisherTest_.compute(m_, l, n_, k).convert_to<double>();
	}
}
