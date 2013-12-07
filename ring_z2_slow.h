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
	struct InPolysSetWithOrigMetadata:NoCopy
	{
		
	};
	
	struct OutPolysSetForVariyingMetadata:NoCopy
	{
		
	};

	std::unique_ptr<OutPolysSetForVariyingMetadata> PrepareEmptyResult();
	
  protected:
	typedef CrossRingInfo::MonomialMetadata<CrossRingInfo::MonomialOrder::DegRevLex> ImplementedOrder;
	struct ImplementedField: FiniteField<ZRing8>
	{
		ImplementedField():
			FiniteField<ZRing8>(2)
		{}
	};

	bool ConstructAndInsertNormalizedImpl(const std::unique_ptr<const InPolysSetWithOrigMetadata>& prepared_input, 
		Enumerator<CrossRingInfo::PerVariableData> top_info,
		Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input_polys_mons, 
		const std::unique_ptr<OutPolysSetForVariyingMetadata>& result);

	void ExtendWithMonomialImpl(Enumerator<CrossRingInfo::PerVariableData> info);

	std::unique_ptr<const InPolysSetWithOrigMetadata> PrepareForReconstructionImpl(Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input);

	void ConvertResultToFixedMetadataImpl(const std::unique_ptr<OutPolysSetForVariyingMetadata>& constructed_result, CrossRingInfo::MonomialListListWithCoef<ImplementedOrder, ImplementedField>& basic_result);

	struct Monomial : std::map<char,int>
	{
		friend bool operator<(const Monomial&, const Monomial&); //undefined
	};
	
	struct Polynomial : std::vector<Monomial>{};
	explicit RingZ2SlowBase(int var_count);
	RingZ2SlowBase(const RingZ2SlowBase&);
	~RingZ2SlowBase();
	struct Impl;
	std::unique_ptr<Impl> impl_;

public:
	struct PolysSet: private std::vector<Polynomial>
	{
		friend class RingZ2SlowBase;
	};

};

//Z_2 ring with degrevlex oredr on variables
template <class MonomialMetadata, class Field>
struct RingZ2Slow: public RingBase<MonomialMetadata, Field, RingZ2Slow<MonomialMetadata, Field>>, public RingZ2SlowBase
{
	typedef  RingBase<MonomialMetadata, Field, RingZ2Slow<MonomialMetadata, Field>> Base;
	RingZ2Slow(const RingZ2Slow& copy_from) = default;
	
	RingZ2Slow(const MonomialMetadata& monomial_metadata, const Field& field)
		: Base(monomial_metadata, field)
		, RingZ2SlowBase(monomial_metadata.var_count)
	{
		assert(field.IsFiniteZpFieldWithChar(2));
		assert(MonomialMetadata::order == CrossRingInfo::MonomialOrder::DegRevLex);
	}
	
	RingZ2Slow& operator=(const RingZ2Slow& copy_from) = delete;
	
	bool ConstructAndInsertNormalized(const std::unique_ptr<const InPolysSetWithOrigMetadata>& prepared_input, const std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>>& info, const std::unique_ptr<OutPolysSetForVariyingMetadata>& result)
	{			
		auto poly_enumerator = FullRangeEnumerator(*info);
		typedef decltype(poly_enumerator.GetAndMove()) PolynomialAsRange;
		auto mon_enumerator = ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullRangeEnumerator<PolynomialAsRange>)>(poly_enumerator);
		typedef decltype(mon_enumerator.GetAndMove().GetAndMove()) MonomialAsRange;

		auto var_enumerator = 
			ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM((
						ConverterEnumeratorCFunc<
						STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullRangeEnumerator<MonomialAsRange>),
						MonomialAsRange>
						))>(mon_enumerator);
		

		return ConstructAndInsertNormalizedImpl(
			prepared_input,
			FullRangeEnumerator(info->TopInfo()),
			var_enumerator,
			result
		);
	}
		
	void ExtendWithMonomial(const std::unique_ptr<const CrossRingInfo::SingleMonomial<MonomialMetadata>>& info)
	{
		ExtendWithMonomialImpl(FullRangeEnumerator(*info));
	}
		
	std::unique_ptr<const InPolysSetWithOrigMetadata> PrepareForReconstruction(const CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field>& input)
	{
		auto poly_enumerator = FullRangeEnumerator(input);
		typedef decltype(poly_enumerator.GetAndMove()) PolynomialAsRange;
		auto mon_enumerator = ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullRangeEnumerator<PolynomialAsRange>)>(poly_enumerator);
		typedef decltype(mon_enumerator.GetAndMove().GetAndMove()) MonomialAsRange;

		//TODO: add checking for coefficient is equal to one
		auto nochecking_var_enumerator = 
			ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM((
						ConverterEnumeratorCFunc<
						STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullRangeEnumerator<MonomialAsRange>),
						MonomialAsRange>
						))>(mon_enumerator);



		return PrepareForReconstructionImpl(nochecking_var_enumerator);
	}
	
	void ConvertResultToFixedMetadata(const std::unique_ptr<OutPolysSetForVariyingMetadata>& constructed_result, std::unique_ptr<const typename Base::IOData::IOPolynomSet>& final_result)
	{
		ImplementedOrder implemented_order;
		implemented_order.var_count = Base::monomial_metadata_.var_count;
		ImplementedField implemented_field;
		CrossRingInfo::MonomialListListWithCoef<ImplementedOrder, ImplementedField> basic_result{implemented_order, implemented_field};
		ConvertResultToFixedMetadataImpl(constructed_result, basic_result);
		auto result_ptr = new typename Base::IOData::IOPolynomSet{Base::monomial_metadata_, Base::field_};
		//convert from ImplementedOrder, ImplementedField to actual
		for(auto poly:basic_result)
		{
			result_ptr->BeginPolynomialConstruction(slow_distance(poly.begin(), poly.end()));
			for(auto mon:poly)
			{
				for(auto var:mon)
				{
					result_ptr->AddVariable(CrossRingInfo::PerVariableData(var.degree, var.index));
				}
				result_ptr->MonomialAdditionDone(mon.coef());
			}
		}
		final_result.reset(result_ptr);
	}
};

struct  FastZ2SlowBasedRingBase:NoCopy
{
	struct Monomial : std::map<char,int>
	{
		friend bool operator<(const Monomial&, const Monomial&); //undefined
	};
	struct FastPoly : std::vector<Monomial> {};
	struct LPolyImpl
	{
		FastPoly value;
		std::vector<FastPoly> reconstruction_info;
		double sig_index;
		Monomial sig_mon;
	};
	struct MultLPoly
	{
		LPolyImpl poly;
		Monomial mul_by;
	};	
};

template <class MonomialMetadata>
class FastZ2SlowBasedRing: private FastZ2SlowBasedRingBase
{
	
public:
	class LPoly:LPolyImpl
	{
		friend class FastZ2SlowBasedRing;
	};
	class MultLPolysQueue:std::vector<MultLPoly>
	{
		friend class FastZ2SlowBasedRing;
	};
	class LPolysResult:std::vector<LPolyImpl>
	{
		friend class FastZ2SlowBasedRing;
	};

	bool QueueEmpty(const MultLPolysQueue& queue)
	{
		return queue.empty();
	}

	FastZ2SlowBasedRing(const MonomialMetadata&){}

	LPoly DequeueSigSmallest(MultLPolysQueue& queue);

	template <class Field>
	MultLPolysQueue PutInQueueExtendLabeledPolys(const CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field>& input);
	//{return MultLPolysQueue();}

	LPolysResult FillWithTrivialSyzygiesOfNonMultElements(const MultLPolysQueue& queue);
	void ReduceCheckingSignatures(LPoly& poly, LPolysResult& reducers);
	
	std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>> FieldAgnosticReconstructionInfo(const LPoly& poly);

	bool IsZero(const LPoly& poly);
	
	void Normalize(LPoly& poly);
	
	std::unique_ptr<const CrossRingInfo::SingleMonomial<MonomialMetadata>> ExtendRingWithMonomialToHelpReconstruct(const LPoly& poly, LPolysResult& reducers);
	void ExtendQueueBySpairPartsAndFilterUnneeded(const LPolysResult& left_parts, const LPoly& right_part, MultLPolysQueue& queue);
	void InsertInResult(const LPoly& poly, LPolysResult& result);
};
