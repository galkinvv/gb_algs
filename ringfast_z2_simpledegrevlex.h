#pragma once
#include <map>
#include <set>
#include <vector>
#include <memory>
#include "cross_ring_info.h"
#include "algs.h"
#include "simple_mon.h"
#include "ringfastbase.h"
#include "z_ring.h"
#include "finite_field.h"
#include "utils.h"

class  RingFastZ2SimpleDegrevlexBase
{
public:
	struct LPolyImpl;
	struct MultLPolysQueueImpl;
	struct LPolysResultImpl;
    typedef unique_created_ptr<LPolyImpl> LPoly;
    typedef unique_created_ptr<MultLPolysQueueImpl> MultLPolysQueue;
    typedef unique_created_ptr<LPolysResultImpl> LPolysResult;

	RingFastZ2SimpleDegrevlexBase();

	bool QueueEmpty(const MultLPolysQueue& queue);

	LPoly DequeueSigSmallest(MultLPolysQueue& queue);

	LPolysResult FillWithTrivialSyzygiesOfNonMultElements(const MultLPolysQueue& queue);
	void ReduceCheckingSignatures(LPoly& poly, LPolysResult& reducers);

	bool IsZero(const LPoly& poly);
	void Normalize(LPoly& poly);

	void ExtendQueueBySpairPartsAndFilterUnneeded(const LPolysResult& left_parts, const LPoly& right_part, MultLPolysQueue& queue);
	void InsertInResult(LPoly&& poly, LPolysResult& result);
protected:
	MultLPolysQueue PutInQueueExtendLabeledPolysImpl(Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input);
	void AddLabeledPolyBeforeImpl(int new_var_index, int new_poly_index_in_rec_basis, Enumerator<CrossRingInfo::PerVariableData> monomial, LPolysResult& reducers, const LPoly& poly_before);
	Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> FieldAgnosticReconstructionInfoPolysImpl(const LPoly& poly);
	Enumerator<CrossRingInfo::PerVariableData> FieldAgnosticReconstructionInfoTopImpl(const LPoly& poly);
private:
	DECLARE_PIMPL;
};

template <class MonomialMetadata>
struct RingFastZ2SimpleDegrevlex final: public RingFastZ2SimpleDegrevlexBase
{
	RingFastZ2SimpleDegrevlex(const MonomialMetadata& metadata)
		:metadata_(metadata)
	{
		assert(MonomialMetadata::order == CrossRingInfo::MonomialOrder::DegRevLex);
	}

	template <class Field>
	MultLPolysQueue PutInQueueExtendLabeledPolys(const CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field>& input)
	{
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

		return PutInQueueExtendLabeledPolysImpl(var_enumerator);
	}

	std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>> FieldAgnosticReconstructionInfo(const LPoly& poly)
	{
		auto result = as_unique_ptr(new CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>(metadata_));
		for (auto top_info_var : FieldAgnosticReconstructionInfoTopImpl(poly))
		{
			result->AddVariable(top_info_var);
		}
		result->TopInfoAdditionDone();
		
		for (auto info_poly : FieldAgnosticReconstructionInfoPolysImpl(poly))
		{
			result->BeginPolynomialConstruction(info_poly.size());
			for (auto info_mon : info_poly)
			{
				for (auto info_var : info_mon)
				{
					result->AddVariable(info_var);
				}
				result->MonomialAdditionDone();
			}
		}

		return MoveToResultType(result);
	}

	void AddLabeledPolyBefore(const std::unique_ptr<const CrossRingInfo::AddedVarInfo<MonomialMetadata>>& added_info, LPolysResult& reducers, const LPoly& poly_before)
	{
		AddLabeledPolyBeforeImpl(added_info->AddedMonomialIndex(), added_info->AddedPolynomialIndex(), FullNonSizedRangeEnumerator(*added_info), reducers, poly_before);
	}

	const MonomialMetadata& metadata_;
};
