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

class RingZ2SlowBase
{
public:
	struct InPolysSetWithOrigMetadata;
	struct OutPolysSetForVariyingMetadata;

	unique_deleter_ptr<OutPolysSetForVariyingMetadata> PrepareEmptyResult();

protected:
	typedef CrossRingInfo::MonomialMetadata<CrossRingInfo::MonomialOrder::DegRevLex> ImplementedOrder;
	struct ImplementedField: FiniteField<ZPlusRing8> {
		ImplementedField():
			FiniteField<ZPlusRing8>(FiniteField<ZPlusRing8>::CreateZpFieldWithChar(2))
		{}
	};

	bool ConstructAndInsertNormalizedImpl(const unique_deleter_ptr<const InPolysSetWithOrigMetadata>& prepared_input,
	                                      Enumerator<CrossRingInfo::PerVariableData> top_info,
	                                      Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input_polys_mons,
	                                      const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& result);

	int ExtendRingWithMonomialToHelpReconstructImpl(Enumerator<CrossRingInfo::PerVariableData> info);

	unique_deleter_ptr<const InPolysSetWithOrigMetadata> PrepareForReconstructionImpl(Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input);

	void ConvertResultToFixedMetadataImpl(const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& constructed_result, CrossRingInfo::MonomialListListWithCoef<ImplementedOrder, ImplementedField>& basic_result);

	struct Monomial : std::map<char,int> {
		friend bool operator<(const Monomial&, const Monomial&); //undefined
	};

	struct Polynomial : std::vector<Monomial> {};
	explicit RingZ2SlowBase(int var_count);

	DECLARE_PIMPL;
public:
	struct PolysSet: private std::vector<Polynomial> {
		friend class RingZ2SlowBase;
	};

};

//Z_2 ring with degrevlex oredr on variables
template <class MonomialMetadata, class Field>
struct RingZ2Slow: public RingBase<MonomialMetadata, Field, RingZ2Slow<MonomialMetadata, Field>>, public RingZ2SlowBase {
	typedef  RingBase<MonomialMetadata, Field, RingZ2Slow<MonomialMetadata, Field>> Base;
	RingZ2Slow(const RingZ2Slow& copy_from) = default;

	RingZ2Slow(const MonomialMetadata& monomial_metadata, const Field& field)
		: Base(monomial_metadata, field)
		, RingZ2SlowBase(monomial_metadata.var_count) {
		assert(field.template ExportZpModulus<int>() == 2);
		assert(MonomialMetadata::order == CrossRingInfo::MonomialOrder::DegRevLex);
	}

	RingZ2Slow& operator=(const RingZ2Slow& copy_from) = delete;

	bool ConstructAndInsertNormalized(const unique_deleter_ptr<const InPolysSetWithOrigMetadata>& prepared_input, const std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>>& info, const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& result) {
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


		return ConstructAndInsertNormalizedImpl(
		           prepared_input,
		           FullNonSizedRangeEnumerator(info->TopInfo()),
		           var_enumerator,
		           result
		       );
	}

	std::unique_ptr<const CrossRingInfo::AddedVarInfo<MonomialMetadata>>  ExtendRingWithMonomialToHelpReconstruct(const std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>>& info){
		auto top_monomial = info->TopInfo();
		int added_var_index = ExtendRingWithMonomialToHelpReconstructImpl(FullNonSizedRangeEnumerator(top_monomial));
		auto result = as_unique_ptr(new CrossRingInfo::AddedVarInfo<MonomialMetadata>(Base::monomial_metadata_, added_var_index));
		for(auto var : top_monomial)
		{
			result->AddVariable(var);
		}
		return MoveToResultType(result);
	}

	unique_deleter_ptr<const InPolysSetWithOrigMetadata> PrepareForReconstruction(const CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field>& input) {
		auto poly_enumerator = FullRangeEnumerator(input);
		typedef decltype(poly_enumerator.GetAndMove()) PolynomialAsRange;
		auto mon_enumerator = ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullRangeEnumerator<PolynomialAsRange>)>(poly_enumerator);
		typedef decltype(mon_enumerator.GetAndMove().GetAndMove()) MonomialAsRange;

		//TODO: add checking for coefficient is equal to one
		auto nochecking_var_enumerator =
		    ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM((
		                ConverterEnumeratorCFunc<
		                STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullNonSizedRangeEnumerator<MonomialAsRange>),
		                MonomialAsRange>
		            ))>(mon_enumerator);



		return PrepareForReconstructionImpl(nochecking_var_enumerator);
	}

	std::unique_ptr<const CrossRingInfo::VariableMapping<MonomialMetadata>> VarMapping()const {
		auto result_ptr =  as_unique_ptr(new CrossRingInfo::VariableMapping<MonomialMetadata>(Base::monomial_metadata_, 0));
		//TODO
		return MoveToResultType(result_ptr);
	}
	
	std::unique_ptr<const typename Base::IOData::IOPolynomSet> ConvertResultToFixedMetadata(const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& constructed_result) {
		ImplementedOrder implemented_order;
		implemented_order.var_count = Base::monomial_metadata_.var_count;
		ImplementedField implemented_field;
		CrossRingInfo::MonomialListListWithCoef<ImplementedOrder, ImplementedField> basic_result {implemented_order, implemented_field};
		ConvertResultToFixedMetadataImpl(constructed_result, basic_result);
		std::unique_ptr<typename Base::IOData::IOPolynomSet> result_ptr {new typename Base::IOData::IOPolynomSet{Base::monomial_metadata_, Base::field_}};
		//convert from ImplementedOrder, ImplementedField to actual
		for(auto poly:basic_result) {
			result_ptr->BeginPolynomialConstruction(slow_distance(poly.begin(), poly.end()));
			for(auto mon:poly) {
				for(auto var:mon) {
					result_ptr->AddVariable(CrossRingInfo::PerVariableData(var.degree, var.index));
				}
				result_ptr->MonomialAdditionDone(mon.coef());
			}
		}
		return MoveToResultType(result_ptr);
	}
};

class  FastZ2SlowBasedRingBase:NoCopy
{
public:
	struct LPoly {
		DECLARE_PIMPL;
	};
	struct MultLPolysQueue {
		DECLARE_PIMPL;
	};
	struct LPolysResult {
		DECLARE_PIMPL;
	};

	bool QueueEmpty(const MultLPolysQueue& queue);

	LPoly DequeueSigSmallest(MultLPolysQueue& queue);

	LPolysResult FillWithTrivialSyzygiesOfNonMultElements(const MultLPolysQueue& queue);
	void ReduceCheckingSignatures(LPoly& poly, LPolysResult& reducers);

	bool IsZero(const LPoly& poly);
	void Normalize(LPoly& poly);

	void ExtendQueueBySpairPartsAndFilterUnneeded(const LPolysResult& left_parts, const LPoly& right_part, MultLPolysQueue& queue);
	void InsertInResult(const LPoly& poly, LPolysResult& result);
protected:
	struct Monomial : std::map<int,int> {
		friend bool operator<(const Monomial&, const Monomial&); //undefined
	};
	struct Signature {
		double sig_index;
		Monomial sig_mon;
	};
	virtual bool MonomialLess(const Monomial& m1, const Monomial& m2) const = 0;
	
	MultLPolysQueue PutInQueueExtendLabeledPolysImpl(Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input);
	void AddLabeledPolyBeforeImpl(int new_var_index, Enumerator<CrossRingInfo::PerVariableData> monomial, LPolysResult& reducers, const LPoly& poly_before);
	Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> FieldAgnosticReconstructionInfoPolysImpl(const LPoly& poly);
	Enumerator<CrossRingInfo::PerVariableData> FieldAgnosticReconstructionInfoTopImpl(const LPoly& poly);
private:
	DECLARE_PIMPL;
};

template <class MonomialMetadata>
class FastZ2SlowBasedRing: public FastZ2SlowBasedRingBase
{
public:

	FastZ2SlowBasedRing(const MonomialMetadata& metadata)
		:metadata_(metadata)
	{}

	template <class Field>
	MultLPolysQueue PutInQueueExtendLabeledPolys(const CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field>& input)
	{
		auto poly_enumerator = FullRangeEnumerator(input);
		typedef decltype(poly_enumerator.GetAndMove()) PolynomialAsRange;
		auto mon_enumerator = ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullRangeEnumerator<PolynomialAsRange>)>(poly_enumerator);
		typedef decltype(mon_enumerator.GetAndMove().GetAndMove()) MonomialAsRange;

		//TODO: add checking for coefficient is not equal to zero
		auto nochecking_var_enumerator =
		    ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM((
		                ConverterEnumeratorCFunc<
		                STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullNonSizedRangeEnumerator<MonomialAsRange>),
		                MonomialAsRange>
		            ))>(mon_enumerator);

		return PutInQueueExtendLabeledPolysImpl(nochecking_var_enumerator);
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
		AddLabeledPolyBeforeImpl(added_info->AddedIndex(), FullNonSizedRangeEnumerator(*added_info), reducers, poly_before);
	}

private:
	virtual bool MonomialLess(const Monomial& m1, const Monomial& m2) const
	{
		//TODO
		return false;
	}
	const MonomialMetadata& metadata_;
};
