#pragma once
#include "utils.h"
#include "cross_ring_info.h"

template <class TMonomialMetadata, class TField>
struct RingBase: private NoCopy
{
	typedef TField Field;
	typedef TMonomialMetadata MonomialMetadata;
	typedef CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field> IOPolynomSet;

	RingBase(const MonomialMetadata& monomial_metadata, const Field& field)
		: monomial_metadata_(monomial_metadata)
		,field_(field)
	{}

	protected:
		MonomialMetadata monomial_metadata_;
		Field field_;
};

template <class DerivedRing>
struct IOData
{
	typedef typename DerivedRing::IOPolynomSet IOPolynomSet;
	
	IOData(const IOPolynomSet& a_in_data, DerivedRing& a_out_ring):
		in_data(a_in_data), out_ring(a_out_ring)
	{}

	const IOPolynomSet& in_data;
	DerivedRing& out_ring;
	unique_deleter_ptr<const IOPolynomSet> out_data;
};
