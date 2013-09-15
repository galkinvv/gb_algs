#pragma once
#include "iopolynomset.h"
namespace F4MPI
{
	template <class TRing, template <class> class TAlgoT>
	PolynomSet GBwithRingAlgo(PolynomSet F, const F4AlgData* /*f4options*/){
		F4MPI::IOPolynomSet io_poly_set_in;
		io_poly_set_in.polys = F;
		io_poly_set_in.type = FieldType::Z;
		io_poly_set_in.field_char = CModular::getMOD();
		io_poly_set_in.mon_order = CMonomial::getOrder();
		auto io_data = TRing::Create(io_poly_set_in);
		TAlgoT<TRing>::Do(*io_data);
		F4MPI::IOPolynomSet io_poly_set_out = TRing::ConvertResult(io_data);
		return io_poly_set_out.polys;
	}
}
	