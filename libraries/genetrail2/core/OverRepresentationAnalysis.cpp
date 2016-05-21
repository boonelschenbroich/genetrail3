/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

double OverRepresentationAnalysis::numberOfHits(const Category& category) const {
	return static_cast<double>(Category::intersect("null", category, test_set_).size());
}

double OverRepresentationAnalysis::expectedNumberOfHits(const Category& category) const {
	auto l = Category::intersect("null", category, reference_set_).size();
	return (l * n_) / static_cast<double>(m_);
}

double OverRepresentationAnalysis::computePValue_(size_t l, size_t k, bool enriched) const {
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


double OverRepresentationAnalysis::computePValue(const Category& category) const
{
	// GeneTrail 1
	// size_t l = category.size();

	size_t k = Category::intersect("null", category, test_set_).size();
	size_t l = Category::intersect("null", category, reference_set_).size();

	auto expected_k = ((double)l * n_) / ((double)m_);
	bool enriched = expected_k < k;

	return computePValue_(l, k, enriched);
}


double OverRepresentationAnalysis::computeUpperTailedPValue(const Category& category) const
{
	// GeneTrail 1
	// size_t l = category.size();

	size_t k = Category::intersect("null", category, test_set_).size();
	size_t l = Category::intersect("null", category, reference_set_).size();

	return computePValue_(l, k, true);
}


double OverRepresentationAnalysis::computeLowerTailedPValue(const Category& category) const
{
	// GeneTrail 1
	// size_t l = category.size();

	size_t k = Category::intersect("null", category, test_set_).size();
	size_t l = Category::intersect("null", category, reference_set_).size();

	return computePValue_(l, k, false);
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
