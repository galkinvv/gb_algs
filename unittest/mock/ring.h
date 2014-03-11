#pragma once
#include <map>
#include <set>
#include <vector>
#include <string>
#include <memory>
#include "mock_base.h"
#include "ringbase.h"
#include "algs.h"
#include "utils.h"
#include "cross_ring_info.h"

namespace Mock
{
	template <class MonomialMetadata, class Field>
	class Ring: public RingBase<MonomialMetadata, Field, Ring<MonomialMetadata, Field>>
	{
		typedef RingBase<MonomialMetadata, Field, Ring<MonomialMetadata, Field>> Base;
	public:

		struct InPolysSetWithOrigMetadata:NoCopy{};
		struct OutPolysSetForVariyingMetadata:NoCopy{};
		Ring(const Ring& copy_from){IgnoreIfUnused(copy_from);}
		Ring(const MonomialMetadata& monomial_metadata, const Field& field):
			Base(monomial_metadata, field)
		{}
		Ring& operator=(const Ring& copy_from) = delete;

		bool ConstructAndInsertNormalized(const std::unique_ptr<InPolysSetWithOrigMetadata>& reconstruction_basis, const std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>>& info, const std::unique_ptr<OutPolysSetForVariyingMetadata>& result)
		{
			IgnoreIfUnused(reconstruction_basis, info, result);
			return true;
		}

		std::unique_ptr<const CrossRingInfo::AddedVarInfo<MonomialMetadata>>  ExtendRingWithMonomialToHelpReconstruct(const std::unique_ptr<InPolysSetWithOrigMetadata>& reconstruction_basis, const std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>>& info)
		{
			IgnoreIfUnused(reconstruction_basis, info);
			return nullptr;
		}

		std::unique_ptr<InPolysSetWithOrigMetadata> PrepareForReconstruction(const CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field>& input)
		{
			IgnoreIfUnused(input);
			return nullptr;
		}
		
		std::unique_ptr<OutPolysSetForVariyingMetadata> PrepareEmptyResult()
		{
			return nullptr;
		}
		
		std::unique_ptr<const typename Base::IOData::IOPolynomSet> ConvertResultToFixedMetadata(const std::unique_ptr<OutPolysSetForVariyingMetadata>& constructed_result)
		{
			IgnoreIfUnused(constructed_result);
			return MoveToResultType(new typename Base::IOData::IOPolynomSet(this->monomial_metadata_, this->field_));
		}
		
		std::unique_ptr<const CrossRingInfo::VariableMapping<MonomialMetadata>> VarMapping()const
		{
			auto mapping = new CrossRingInfo::VariableMapping<MonomialMetadata>(this->monomial_metadata_, this->monomial_metadata_.var_count);
			for (int i =0; i < this->monomial_metadata_.var_count; ++i)
			{
				mapping->AddVariable(CrossRingInfo::PerVariableData::FromDI(1, i));
				mapping->MonomialAdditionDone();
			}
			return MoveToResultType(mapping);
		}
	};

	template <class MonomialMetadata>
	class FastRingWithTracking : NoCopy
	{
	public:
		struct LPoly{};
		struct MultLPolysQueue{};
		struct LPolysResult{};

		FastRingWithTracking(const MonomialMetadata& metadata)
		{
			IgnoreIfUnused(metadata);
		}
		bool QueueEmpty(const MultLPolysQueue& queue)
		{
			IgnoreIfUnused(queue);
			return true;
		}
		LPoly DequeueSigSmallest(MultLPolysQueue& queue)
		{
			IgnoreIfUnused(queue);
			return LPoly();
		}
		
		template <class Field>
		MultLPolysQueue PutInQueueExtendLabeledPolys(const CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field>& input)
		{
			IgnoreIfUnused(input);
			return MultLPolysQueue();
		}
		LPolysResult FillWithTrivialSyzygiesOfNonMultElements(const MultLPolysQueue& queue)
		{
			IgnoreIfUnused(queue);
			return LPolysResult();
		}
		void ReduceCheckingSignatures(LPoly& poly, LPolysResult& reducers)
		{
			IgnoreIfUnused(poly, reducers);
		}

		std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>>  FieldAgnosticReconstructionInfo(const LPoly& poly)
		{
			IgnoreIfUnused(poly);
			return nullptr;
		}
		bool IsZero(const LPoly& poly)
		{
			IgnoreIfUnused(poly);
			return true;
		}
		void Normalize(LPoly& poly)
		{
			IgnoreIfUnused(poly);
		}
		void AddLabeledPolyBefore(const std::unique_ptr<const CrossRingInfo::AddedVarInfo<MonomialMetadata>>& added_info, LPolysResult& reducers, const LPoly& poly_before)
		{
			IgnoreIfUnused(added_info, reducers, poly_before);
		}
		void ExtendQueueBySpairPartsAndFilterUnneeded(const LPolysResult& left_parts, const LPoly& right_part, MultLPolysQueue& queue)
		{
			IgnoreIfUnused(left_parts, right_part, queue);
		}
		void InsertInResult(LPoly&& poly, LPolysResult& result)
		{
			IgnoreIfUnused(poly, result);
		}
	};
}
