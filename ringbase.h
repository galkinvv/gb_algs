#pragma once
#include "utils.h"
#include "cross_ring_info.h"
template <class TMonomialMetadata, class TField>
struct Ring: NoCopy
{
	typedef TField Field;
	typedef TMonomialMetadata MonomialMetadata;

	struct IOData
	{
		typedef CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field> IOPolynomSet;
		
		IOData(const Ring& a_in_ring, const IOPolynomSet& a_in_data):
			in_ring(a_in_ring), in_data(a_in_data)
		{}
	
		const Ring& in_ring;
		const IOPolynomSet& in_data;
		Ring out_ring;
		IOPolynomSet out_data;
	};
};
