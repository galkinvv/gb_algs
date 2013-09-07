#pragma once

#include "settings.h"
#include "types.h"
namespace F4MPI{
namespace PerryOrig
{
	PolynomSet F5ByPerry(PolynomSet F, const F4AlgData* f4options);
}
namespace F5Orig
{
	PolynomSet IncrementalF5(PolynomSet F, const F4AlgData* f4options);
}

	
} //namespace F4MPI

