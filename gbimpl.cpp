#include "gbimpl.h"
#include "f4main.h"
#include "f5main.h"
#include <iterator>
#include <stdexcept>

namespace F4MPI{
#define StaticArraySize(arr) (sizeof(arr)/sizeof(arr[0]))
	/*
	template <class T, int N> int StaticArraySize(T arr[N])
	{
		return N;
	}
*/
	PolynomSet GB(PolynomSet F, const F4AlgData* f4options)
	{
		const decltype(GB) *funcs[] = {F4, F5Orig::IncrementalF5, PerryOrig::F5ByPerry};
		if (StaticArraySize(funcs) > f4options->selectedAlgo)
		{
			return funcs[f4options->selectedAlgo](F, f4options);
		}
		throw std::runtime_error("Unknown algo selected");
	}
}