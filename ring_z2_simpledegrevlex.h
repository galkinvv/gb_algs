#pragma once
#include <map>
#include <set>
#include <vector>
#include <memory>
#include "algs.h"
#include "ringbase.h"
#include "z_ring.h"
#include "finite_field.h"
#include "utils.h"

class RingZ2SimpleDegrevlexBase
{
public:
	struct InPolysSetWithOrigMetadata;
	struct OutPolysSetForVariyingMetadata;

	unique_deleter_ptr<OutPolysSetForVariyingMetadata> PrepareEmptyResult();

protected:
	typedef CrossRingInfo::MonomialMetadata<CrossRingInfo::MonomialOrder::DegRevLex> ImplementedOrder;
	struct ImplementedField: FiniteField<ZPlusRing8> {
		ImplementedField():
			FiniteField<ZPlusRing8>(FiniteField<ZPlusRing8>::CreateZpFieldWithChar(2u))
		{}
	};

	bool ReconstructAndInsertNormalizedImpl(const InPolysSetWithOrigMetadata& reconstruction_basis,
	                                      Enumerator<CrossRingInfo::PerVariableData> top_info,
	                                      Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input_polys_mons,
	                                      const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& result);

	struct NewIndices{
		int new_var_index, new_poly_index;
	};
	
	NewIndices ExtendRingWithMonomialToHelpReconstructImpl(const unique_deleter_ptr<InPolysSetWithOrigMetadata>& reconstruction_basis, Enumerator<CrossRingInfo::PerVariableData> info);

	unique_deleter_ptr<InPolysSetWithOrigMetadata> PrepareForReconstructionImpl(Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input);

	//should return only newly added vars
	int VarMappingImplReturningOldVarCount(std::vector<int>& new_monomial_vars) const;

	void ConvertResultToFixedMetadataImpl(const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& constructed_result, CrossRingInfo::MonomialListListWithCoef<ImplementedOrder, ImplementedField>& basic_result);

	explicit RingZ2SimpleDegrevlexBase(int var_count);

	DECLARE_PIMPL;
};

//Z_2 ring with degrevlex oredr on variables
template <class MonomialMetadata, class Field>
struct RingZ2SimpleDegrevlex final: public RingBase<MonomialMetadata, Field>, public RingZ2SimpleDegrevlexBase {
	typedef  RingBase<MonomialMetadata, Field> Base;

	RingZ2SimpleDegrevlex(const MonomialMetadata& monomial_metadata, const Field& field)
		: Base(monomial_metadata, field)
		, RingZ2SimpleDegrevlexBase(monomial_metadata.var_count) {
		assert(field.template ExportZpModulus<int>() == 2);
		assert(MonomialMetadata::order == CrossRingInfo::MonomialOrder::DegRevLex);
	}

	bool ReconstructAndInsertNormalized(const unique_deleter_ptr<InPolysSetWithOrigMetadata>& reconstruction_basis, const std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>>& info, const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& result) {
		auto poly_enumerator = FullRangeEnumerator(*info);
		typedef decltype(poly_enumerator.GetAndMove()) PolynomialAsRange;
		auto mon_enumerator = ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullRangeEnumerator<PolynomialAsRange>)>(poly_enumerator);
		typedef decltype(mon_enumerator.GetAndMove().GetAndMove()) MonomialAsRange;

		auto var_enumerator =
		    ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM((
		                ConverterEnumeratorCFunc<
		                STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullNonSizedRangeEnumerator<MonomialAsRange>),
		                MonomialAsRange>
		            ))>(mon_enumerator);


		return ReconstructAndInsertNormalizedImpl(
		           *reconstruction_basis,
		           FullNonSizedRangeEnumerator(info->TopInfo()),
		           var_enumerator,
		           result
		       );
	}

	std::unique_ptr<const CrossRingInfo::AddedVarInfo<MonomialMetadata>>  ExtendRingWithMonomialToHelpReconstruct(const unique_deleter_ptr<InPolysSetWithOrigMetadata>& reconstruction_basis, const std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>>& info){
		auto top_monomial = info->TopInfo();
		const NewIndices added_indices = ExtendRingWithMonomialToHelpReconstructImpl(reconstruction_basis, FullNonSizedRangeEnumerator(top_monomial));
		auto result = as_unique_ptr(new CrossRingInfo::AddedVarInfo<MonomialMetadata>(Base::monomial_metadata_, added_indices.new_var_index, added_indices.new_poly_index));
		for(auto var : top_monomial)
		{
			result->AddVariable(var);
		}
		return MoveToResultType(result);
	}

	unique_deleter_ptr<InPolysSetWithOrigMetadata> PrepareForReconstruction(const CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field>& input) {
		auto poly_enumerator = FullRangeEnumerator(input);
		typedef decltype(poly_enumerator.GetAndMove()) PolynomialAsRange;
		auto nochecking_mon_enumerator = ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullRangeEnumerator<PolynomialAsRange>)>(poly_enumerator);

		auto polynomial_coef_nonzero_checker = ConverterOfInnerEnumerator(MonomialCoefNonZeroEnsurer(input.Field()));
		
		auto mon_enumerator = ConverterEnumerator(nochecking_mon_enumerator, polynomial_coef_nonzero_checker);
		typedef decltype(mon_enumerator.GetAndMove().GetAndMove()) MonomialAsRange;

		auto var_enumerator =
		    ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM((
		                ConverterEnumeratorCFunc<
		                STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullNonSizedRangeEnumerator<MonomialAsRange>),
		                MonomialAsRange>
		            ))>(mon_enumerator);


		return PrepareForReconstructionImpl(var_enumerator);
	}

	std::unique_ptr<const CrossRingInfo::VariableMapping<MonomialMetadata>> VarMapping()const {
		std::vector<int> new_monomil_vars;
		int old_var_count = VarMappingImplReturningOldVarCount(new_monomil_vars);
		assert(0 == (new_monomil_vars.size() % old_var_count));
		int new_var_count = new_monomil_vars.size() / old_var_count;
		auto result_ptr = as_unique_ptr(new CrossRingInfo::VariableMapping<MonomialMetadata>(Base::monomial_metadata_, new_var_count));
		//copy old variables as is
		for (int old_var_idx = 0; old_var_idx < old_var_count; ++old_var_idx)
		{
			result_ptr->AddVariable(CrossRingInfo::PerVariableData::FromDI(1, old_var_idx));
			result_ptr->MonomialAdditionDone();
		}
		//add new variables to the end
		for (int new_var_idx = 0; new_var_idx < new_var_count; ++new_var_idx)
		{
			for (int old_var_idx = 0; old_var_idx < old_var_count; ++old_var_idx)
			{
				int deg = new_monomil_vars[(new_var_idx * old_var_count) + old_var_idx];
				if (deg > 0)
				{
					result_ptr->AddVariable(CrossRingInfo::PerVariableData::FromDI(deg, old_var_idx));
				}
			}
			result_ptr->MonomialAdditionDone();
		}
		return MoveToResultType(result_ptr);
	}
	
	unique_deleter_ptr<const typename Base::IOPolynomSet> ConvertResultToFixedMetadata(const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& constructed_result) {
		ImplementedOrder implemented_order;
		implemented_order.var_count = Base::monomial_metadata_.var_count;
		ImplementedField implemented_field;
		CrossRingInfo::MonomialListListWithCoef<ImplementedOrder, ImplementedField> basic_result {implemented_order, implemented_field};
		ConvertResultToFixedMetadataImpl(constructed_result, basic_result);
		std::unique_ptr<typename Base::IOPolynomSet> result_ptr {new typename Base::IOPolynomSet{Base::monomial_metadata_, Base::field_}};
		//convert from ImplementedOrder, ImplementedField to actual
		for(auto poly:basic_result) {
			result_ptr->BeginPolynomialConstruction(slow_distance(poly.begin(), poly.end()));
			for(auto mon:poly) {
				for(auto var:mon) {
					result_ptr->AddVariable(CrossRingInfo::PerVariableData::FromDI(var.degree, var.index));
				}
				result_ptr->MonomialAdditionDone(mon.coef());
			}
		}
		return MoveToResultType(result_ptr.release());
	}
};

