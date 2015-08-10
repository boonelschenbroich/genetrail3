#ifndef GT2_P_VALUE_COMPUTER_H
#define GT2_P_VALUE_COMPUTER_H

#include <memory>

namespace GeneTrail
{
	class PValueComputer
	{
	};

	using PValueComputerPtr = std::unique_ptr<PValueComputer>;

	namespace internal{

	class PValueComputerWorker : PValueComputer
	{
	};
	}

	template<typename Ts...>
	PValueComputerPtr createPValueComputer();
}

#endif // GT2_P_VALUE_COMPUTER_H