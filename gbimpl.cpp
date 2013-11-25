#include "gbimpl.h"
#include "f4main.h"
#include "f5main.h"
#include "algs.h"
#include "f4mpi_to_native_bridge.h"
#include "ring_z2_slow.h"
#include "ssg_approx.h"
#include <iterator>
#include <stdexcept>

namespace F4MPI{
	PolynomSet GB(PolynomSet F, const F4AlgData* f4options)
	{
		typedef  RingZ2Slow<CrossRingInfo::MonomialMetadata<CrossRingInfo::MonomialOrder::DegRevLex>, FiniteField<ZField8>> UsedRingZ2;
		const decltype(GB) *funcs[] = {
			F4,
			F5Orig::IncrementalF5,
			PerryOrig::F5ByPerry,
			GBwithRingAlgo<UsedRingZ2, FastZ2SlowBasedRing, ApproxSignatureGroebner>
		};
		if (countof(funcs) > f4options->selectedAlgo)
		{
			return funcs[f4options->selectedAlgo](F, f4options);
		}
		throw std::runtime_error("Unknown algo selected");
	}
}
