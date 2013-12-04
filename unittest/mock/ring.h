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

		bool ConstructAndInsertNormalized(const std::unique_ptr<const InPolysSetWithOrigMetadata>& prepared_input, const std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>>& info, const std::unique_ptr<OutPolysSetForVariyingMetadata>& result)
		{
			IgnoreIfUnused(prepared_input, info, result);
			return true;
		}

		void ExtendWithMonomial(const std::unique_ptr<const CrossRingInfo::SingleMonomial<MonomialMetadata>>& info)
		{
			IgnoreIfUnused(info);
		}

		std::unique_ptr<const InPolysSetWithOrigMetadata> PrepareForReconstruction(const CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field>& input)
		{
			IgnoreIfUnused(input);
			return nullptr;
		}
		
		std::unique_ptr<OutPolysSetForVariyingMetadata> PrepareEmptyResult()
		{
			return nullptr;
		}
		
		void ConvertResultToFixedMetadata(const std::unique_ptr<OutPolysSetForVariyingMetadata>& constructed_result, std::unique_ptr<const typename Base::IOData::IOPolynomSet>& final_result)
		{
			IgnoreIfUnused(constructed_result, final_result);
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
		std::unique_ptr<const CrossRingInfo::SingleMonomial<MonomialMetadata>>  ExtendRingWithMonomialToHelpReconstruct(const LPoly& poly, LPolysResult& reducers)
		{
			IgnoreIfUnused(poly, reducers);
			return nullptr;
		}
		void ExtendQueueBySpairPartsAndFilterUnneeded(const LPolysResult& left_parts, const LPoly& right_part, MultLPolysQueue& queue)
		{
			IgnoreIfUnused(left_parts, right_part, queue);
		}
		void InsertInResult(const LPoly& poly, LPolysResult& result)
		{
			IgnoreIfUnused(poly, result);
		}
	};
}
