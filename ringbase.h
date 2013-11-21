#pragma once
#include "utils.h"
#include "cross_ring_info.h"
template <class TMonomialMetadata, class TField, class DerivedRing>
struct RingBase: NoCopy
{
	typedef TField Field;
	typedef TMonomialMetadata MonomialMetadata;

	struct IOData
	{
		typedef CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field> IOPolynomSet;
		
		IOData(const DerivedRing& a_in_ring, const IOPolynomSet& a_in_data, DerivedRing& a_out_ring):
			in_ring(a_in_ring), in_data(a_in_data), out_ring(a_out_ring)
		{}
	
		const DerivedRing& in_ring;
		const IOPolynomSet& in_data;
		DerivedRing& out_ring;
		std::unique_ptr<const IOPolynomSet> out_data;
	};
};
