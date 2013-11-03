#pragma once
#include <map>
#include <set>
#include <vector>
#include <string>
#include <memory>
#include "mock_base.h"
#include "ringbase.h"
#include "algs.h"
#include "cross_ring_info.h"

namespace Mock
{
	template <class MonomialMetadata>
	class Ring: NoCopy
	{

		Ring(){}
		friend class IOData<Ring>;

	public:

		struct PolysSet:NoCopy{};
		
		void CopyTo(Ring& /*to*/)const{}

		bool ConstructAndInsertNormalized(const PolysSet& /*in*/, const std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>>&/*info*/, PolysSet& /*out*/)
		{
			//TODO
			return true;
		}

		void ExtendWithMonomial(const std::unique_ptr<const CrossRingInfo::SingleMonomial<MonomialMetadata>>&/*info*/)
		{
		}

		std::unique_ptr<const CrossRingInfo::MonomialListList<MonomialMetadata>> GetCrossRingInfoForInput(const PolysSet& /*in*/)const
		{
			return nullptr;
		}
		static std::unique_ptr<IOData<Ring>> Create(const std::string& /*in*/);

		static std::string ConvertResult(std::unique_ptr<IOData<Ring>>& /*result*/)
		{
			return "";
		}
	};

	template <class MonomialMetadata>
	class FastRingWithTracking : NoCopy
	{
	public:
		struct LPoly{};
		struct MultLPolysQueue{};
		struct	LPolysResult{};

		FastRingWithTracking()
		{}
		bool QueueEmpty(const MultLPolysQueue& /*queue*/)
		{
			return true;
		}
		LPoly DequeueSigSmallest(MultLPolysQueue& /*queue*/)
		{
			return LPoly();
		}
		MultLPolysQueue PutInQueueExtendLabeledPolys(const std::unique_ptr<const CrossRingInfo::MonomialListList<MonomialMetadata>>& /*in*/)
		{
			return MultLPolysQueue();
		}
		LPolysResult FillWithTrivialSyzygiesOfNonMultElements(const MultLPolysQueue& /*queue*/)
		{
			return LPolysResult();
		}
		void ReduceCheckingSignatures(LPoly& /*poly*/, LPolysResult& /*reducers*/)
		{
		}
		std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>>  FieldAgnosticReconstructionInfo(const LPoly& /*poly*/)
		{
			return nullptr;
		}
		bool IsZero(const LPoly& /*poly*/)
		{
			return true;
		}
		void Normalize(LPoly& /*poly*/)
		{
		}
		std::unique_ptr<const CrossRingInfo::SingleMonomial<MonomialMetadata>>  ExtendRingWithMonomialToHelpReconstruct(const LPoly& /*poly*/, LPolysResult& /*reducers*/)
		{
			return nullptr;
		}
		void ExtendQueueBySpairPartsAndFilterUnneeded(const LPolysResult& /*left_parts*/, const LPoly& /*right_part*/, MultLPolysQueue& /*queue*/)
		{
		}
		void InsertInResult(const LPoly& /*poly*/, LPolysResult& /*result*/)
		{
		}
	};

	typedef CrossRingInfo::MonimailMetaData<CrossRingInfo::MonomialOrder::DegRevLex> MockMonOrder;
	struct RingIOData: IOData<Ring<MockMonOrder>> {};
	template<>
	inline std::unique_ptr<IOData<Ring<MockMonOrder>>> Ring<MockMonOrder>::Create(const std::string& /*in*/)
	{
		auto* data = new RingIOData();
		return std::unique_ptr<IOData<Ring<MockMonOrder>>>(data);
	}
}
