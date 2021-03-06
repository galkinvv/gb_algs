#pragma once
#include <cassert>
#include <limits>
#include "types.h"
#include "cross_ring_info.h"
#include "cmonomial.h"
#include "finite_field.h"
#include "ringbase.h"
#include "z_ring.h"

namespace F4MPI
{
	template <class Value, class IOPolynomSet>
	void ConvertF4MPIInputData(const PolynomSet& in, IOPolynomSet& out)
	{
		Value coeff_as_value;
		for(auto poly:in) {
			out.BeginPolynomialConstruction(poly.size());
			auto coeff = poly.c_begin();
			for(auto mon = poly.m_begin(); mon != poly.m_end(); ++mon, ++coeff) {
				for (int var_idx = 0; var_idx < CMonomial::theNumberOfVariables; ++var_idx) {
					if (int deg = mon->getDegree(var_idx)) {
						out.AddVariable(CrossRingInfo::PerVariableData::FromDI(deg, var_idx));
					}
				}
				out.Field().Import(unsigned_cast(coeff->toint()), coeff_as_value);
				out.MonomialAdditionDone(coeff_as_value);
			}
		}
	}

	template <class IOPolynomSet, class VariableMapping>
	void CreateF4MPIResult(const IOPolynomSet& in, const VariableMapping& new2old, PolynomSet& out)
	{
		for(auto poly_in_ring:in) {
			F4MPI::CPolynomial poly;
			for(auto mon_in_ring: poly_in_ring) {
				std::vector<F4MPI::CMonomialBase::Deg> mon(CMonomial::theNumberOfVariables);
				for (auto var:mon_in_ring) {
					for (auto old_var:new2old[var.index])
					{
						mon[old_var.index] += var.degree * old_var.degree;
					}
				}
				int exported_coef = in.Field().template Export<int>(mon_in_ring.coef());
				poly.addTerm(F4MPI::CModular(exported_coef), F4MPI::CMonomial(mon));
			}
			out.push_back(poly);
		}
	}

	template <class TRing, template <class Metadata> class TFastRing, template <class Ring, template <class Metadata> class FastRing> class TAlgoT>
	PolynomSet GBwithRingAlgo(PolynomSet F, const F4AlgData* /*f4options*/){
		typedef typename TRing::IOPolynomSet IOPolynomSet;
		typedef typename TRing::Field Field;
		typedef typename TRing::MonomialMetadata MonomialMetadata;
		
		//to handle other orders and felds add template call
		assert(CMonomial::getOrder() == CMonomial::degrevlexOrder);
		assert(MonomialMetadata::order == CrossRingInfo::MonomialOrder::DegRevLex);

		MonomialMetadata monomial_metadata;
		monomial_metadata.var_count = CMonomial::theNumberOfVariables;
		Field field = Field::CreateZpFieldWithChar(unsigned_cast(CModular::getMOD()));
		TRing out_ring{monomial_metadata, field};
		IOPolynomSet io_poly_set_in {monomial_metadata, field};
		ConvertF4MPIInputData<typename Field::Value>(F, io_poly_set_in);
		IOData<TRing> io_data {io_poly_set_in, out_ring};
		TAlgoT<TRing, TFastRing>::Do(io_data);
		PolynomSet result;
		CreateF4MPIResult(*io_data.out_data, *io_data.out_ring.VarMapping(), result);
		return result;
	}
}
