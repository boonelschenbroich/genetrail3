#include "Parameters.h"

namespace GeneTrail {

Params::Params()
	: significance(0.05),
	minimum(2),
	maximum(700),
	numPermutations(100000),
	randomSeed(std::random_device{}()),
	adjustSeparately(false),
	pValueMode(PValueMode::RowWise)
	{
	}

}