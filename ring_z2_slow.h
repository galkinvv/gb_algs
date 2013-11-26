#pragma once
#include <map>
#include <set>
#include <vector>
#include <memory>
#include "algs.h"
#include "ringbase.h"
#include "z_field.h"
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
	bool ConstructAndInsertNormalizedImpl(const std::unique_ptr<const InPolysSetWithOrigMetadata>& prepared_input, 
		const Enumerator<CrossRingInfo::PerVariableData>& top_info,  
		const Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>>& input_polys_mons, 
		const std::unique_ptr<OutPolysSetForVariyingMetadata>& result);

	//Z_2 ring with degrevlex oredr on variables
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

template <class MonomialMetadata, class Field>
struct RingZ2Slow: public RingBase<MonomialMetadata, Field, RingZ2Slow<MonomialMetadata, Field>>, public RingZ2SlowBase
{
	typedef  RingBase<MonomialMetadata, Field, RingZ2Slow<MonomialMetadata, Field>> Base;
	RingZ2Slow(const RingZ2Slow& copy_from)
		:RingZ2SlowBase(copy_from)
	{}
	
	RingZ2Slow(const MonomialMetadata& monomial_metadata, const Field& field):
		RingZ2SlowBase(monomial_metadata.var_count)
	{
		assert(field.IsFiniteZpFieldWithChar(2));
		assert(MonomialMetadata::order ==  CrossRingInfo::MonomialOrder::DegRevLex);
	}
	
	RingZ2Slow& operator=(const RingZ2Slow& copy_from) = delete;
	
	bool ConstructAndInsertNormalized(const std::unique_ptr<const InPolysSetWithOrigMetadata>& prepared_input, const std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>>& info, const std::unique_ptr<OutPolysSetForVariyingMetadata>& result)
	{			
		auto poly_enumerator = FullRangeEnumerator(*info);
		typedef typename decltype(poly_enumerator)::value_type IterablePolySet;
		auto mon_enumerator = ConverterEnumeratorCFunc<FUNCTION_WITHTYPE_AS_TEMPLATE_PARAM(FullRangeEnumerator<IterablePolySet>)>(poly_enumerator);
		typedef typename decltype(mon_enumerator.GetAndMove())::value_type IterablePoly;

		auto var_enumerator = 
			ConverterEnumeratorCFunc<FUNCTION_WITHTYPE_AS_TEMPLATE_PARAM((
						ConverterEnumeratorCFunc<
						FUNCTION_WITHTYPE_AS_TEMPLATE_PARAM(FullRangeEnumerator<IterablePoly>),
						decltype(FullRangeEnumerator<IterablePoly>(std::declval<IterablePoly>()))>
						))>(mon_enumerator);
		

		return ConstructAndInsertNormalizedImpl(
			prepared_input,
			FullRangeEnumerator(info->TopInfo()),
			mon_enumerator,
			result
		);
	}
		
	void ExtendWithMonomial(const std::unique_ptr<const CrossRingInfo::SingleMonomial<MonomialMetadata>>& info)
	{
	}
		
	std::unique_ptr<const InPolysSetWithOrigMetadata> PrepareForReconstruction(const CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field>& input)
	{
		return nullptr;
	}
	
	void ConvertResultToFoxedMetadata(const std::unique_ptr<OutPolysSetForVariyingMetadata>& constructed_result, std::unique_ptr<const typename Base::IOData::IOPolynomSet>& final_result)
	{
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
