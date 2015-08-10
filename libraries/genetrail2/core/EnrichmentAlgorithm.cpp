#include "EnrichmentAlgorithm.h"

namespace GeneTrail
{
	EnrichmentAlgorithm::EnrichmentAlgorithm(PValueMode mode) : mode_(mode) {}

	bool EnrichmentAlgorithm::pValuesComputed() const
	{
		return mode_ == PValueMode::RowWise && rowWisePValueIsDirect();
	}
}
