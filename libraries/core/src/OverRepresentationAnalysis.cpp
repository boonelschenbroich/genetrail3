#include "OverRepresentationAnalysis.h"

using namespace GeneTrail;
using namespace boost::multiprecision;

bool OverRepresentationAnalysis::categoryContainsAllGenes(const Category& reference, const Category& testSet){
    for (auto gene = testSet.begin(); gene != testSet.end(); ++gene) {
        if (!reference.contains(*gene)) {
            return false;
        }
    }
    return true;
}

std::tuple<double, double, std::string> OverRepresentationAnalysis::computePValue(const Category& category, const Category& referenceSet, const Category& testSet) {

	FishersExactTest<uint64_t, cpp_dec_float_50> fisherTest_;
	HypergeometricTest<uint64_t, cpp_dec_float_50> hyperTest_;

	uint64_t m = referenceSet.size();
	// GeneTrail 1
	//uint64_t l = category.size();
	uint64_t l = intersect("null", referenceSet, category).size();
    uint64_t n = testSet.size();
	Category intersection = intersect("null", testSet, category);
    uint64_t k = intersection.size();
	auto expected_k = ((double)l * n) / ((double)m);

	cpp_dec_float_50 p;

    if (categoryContainsAllGenes(referenceSet, testSet)) {
        if (expected_k < k) {
            p = hyperTest_.upperTailedPValue(m, l, n, k);
        } else {
            p = hyperTest_.lowerTailedPValue(m, l, n, k);
        }
    } else {
        if (expected_k < k) {
            p = fisherTest_.upperTailedPValue(m, l, n, k);
        } else {
            p = fisherTest_.lowerTailedPValue(m, l, n, k);
        }
    }

	std::string genes = "";
	for (auto gene = intersection.begin(); gene != intersection.end(); ++gene) {
		genes += (gene == intersection.begin()) ? *gene : "," + *gene;
	}

    return std::make_tuple(p.convert_to<double>(), expected_k, genes);
}
