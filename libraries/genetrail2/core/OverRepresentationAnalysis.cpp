#include "OverRepresentationAnalysis.h"

using namespace GeneTrail;

OverRepresentationAnalysis::OverRepresentationAnalysis(const Category& reference_set,const Category& test_set)
:reference_set_(reference_set),test_set_(test_set)
{
	m_ = reference_set_.size();
	n_ = test_set_.size();
	useHypergeometricTest_ = categoryContainsAllGenes(reference_set_, test_set_);
}

bool OverRepresentationAnalysis::categoryContainsAllGenes(const Category& reference, const Category& testSet){
    for (auto gene = testSet.begin(); gene != testSet.end(); ++gene) {
        if (!reference.contains(*gene)) {
            return false;
        }
    }
    return true;
}

ORAResult OverRepresentationAnalysis::computePValue(const Category& category) {

	ORAResult result;
	result.reference = category.reference();
	result.name = category.name();

	// GeneTrail 1
	//uint64_t l = category.size();

	uint64_t l = Category::intersect("null", reference_set_, category).size();
	Category intersection = Category::intersect("null", test_set_, category);
	uint64_t k = intersection.size();
	result.hits = k;
	auto expected_k = ((double)l * n_) / ((double)m_);
	result.expected_hits = expected_k;
	bool enriched = expected_k < k;
	result.enriched = enriched;

	big_float p;

	if (useHypergeometricTest_) {
		if (enriched) {
			p = hyperTest_.upperTailedPValue(m_, l, n_, k);
		} else {
			p = hyperTest_.lowerTailedPValue(m_, l, n_, k);
		}
	} else {
		if (enriched) {
			p = fisherTest_.upperTailedPValue(m_, l, n_, k);
		} else {
			p = fisherTest_.lowerTailedPValue(m_, l, n_, k);
		}
	}

	result.pvalue = p;

	std::string genes = "";
	for (auto gene = intersection.begin(); gene != intersection.end(); ++gene) {
		genes += (gene == intersection.begin()) ? *gene : "," + *gene;
	}
	result.info = genes;

    	return result;
}
