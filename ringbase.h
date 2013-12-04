#pragma once
#include "utils.h"
#include "cross_ring_info.h"
template <class TMonomialMetadata, class TField, class DerivedRing>
struct RingBase
{
	RingBase& operator=(const RingBase&) = delete;
	RingBase(const RingBase&) = default;
	typedef TField Field;
	typedef TMonomialMetadata MonomialMetadata;


	RingBase(const MonomialMetadata& monomial_metadata, const Field& field)
		: monomial_metadata_(monomial_metadata)
		,field_(field)
	{}

	struct IOData
	{
		typedef CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field> IOPolynomSet;
		
		IOData(const IOPolynomSet& a_in_data, DerivedRing& a_out_ring):
			in_data(a_in_data), out_ring(a_out_ring)
		{}
	
		const IOPolynomSet& in_data;
		DerivedRing& out_ring;
		std::unique_ptr<const IOPolynomSet> out_data;
	};
	
	protected:
		MonomialMetadata monomial_metadata_;
		Field field_;
};
