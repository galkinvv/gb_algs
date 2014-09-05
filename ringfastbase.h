#pragma once
#include "utils.h"

template <class TMonomialMetadata>
struct RingFastBase: private NoCopy
{
	typedef TMonomialMetadata MonomialMetadata;
	typedef CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field> IOPolynomSet;

	RingFastBase(const MonomialMetadata& monomial_metadata)
		: monomial_metadata_(monomial_metadata)
	{}

	protected:
		MonomialMetadata monomial_metadata_;
};
