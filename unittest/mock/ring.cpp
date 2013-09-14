#include "ring.h"
namespace Mock
{
	std::unique_ptr<Ring> CreateRing()
	{
		std::unique_ptr<Ring> result;
		result.reset(new Ring());
		return result;
	}
}
