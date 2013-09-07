#include "gbimpl.h"
#include "f4main.h"
#include "f5main.h"
#include <iterator>
#include <stdexcept>

namespace F4MPI{
	PolynomSet GB(PolynomSet F, const F4AlgData* f4options)
	{
		const decltype(GB) *funcs[] = {F4, F5Orig::IncrementalF5, PerryOrig::F5ByPerry};
		if (std::end(funcs) - std::begin(funcs) > f4options->selectedAlgo)
		{
			return funcs[f4options->selectedAlgo](F, f4options);
		}
		throw std::runtime_error("Unknown algo selected");
	}
}