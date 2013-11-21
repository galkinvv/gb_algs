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

namespace F4MPI
{
	class IOPolynomSet;
}

class RingZ2SlowBase
{
  public:
	struct InPolysSetWithOrigMetadata:NoCopy
	{
		
	};
	
	struct OutPolysSetForVariyingMetadata:NoCopy
	{
		
	};

	std::unique_ptr<OutPolysSetForVariyingMetadata> PrepareEmptyResult()
	{
		return nullptr;
	}
	
  protected:
	//Z_2 ring with degrevlex oredr on variables
	struct Monomial : std::map<char,int>
	{
		friend bool operator<(const Monomial&, const Monomial&); //undefined
	};
	
	struct Polynomial : std::vector<Monomial>{};
	RingZ2Slow();
	~RingZ2Slow();
	friend class IOData<RingZ2Slow>;
	struct Impl;
	std::unique_ptr<Impl> impl_;
	struct ReconstructionInfoImpl: std::vector<std::vector<Monomial>>
	{
		Monomial top;
	};
public:
	class FastAssociatedLabeledRingWithTracking;
	class ReconstructionInfo: ReconstructionInfoImpl
	{
		friend class RingZ2Slow;
		friend class FastAssociatedLabeledRingWithTracking;
	};

	void CopyTo(RingZ2Slow& to)const;

	struct PolysSet: private std::vector<Polynomial>
	{
		friend class RingZ2Slow;
	};

	class FastAssociatedLabeledRingWithTracking : NoCopy
	{
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
		
	public:
		class LPoly:LPolyImpl
		{
			friend class FastAssociatedLabeledRingWithTracking;
		};
		class MultLPolysQueue:std::vector<MultLPoly>
		{
			friend class FastAssociatedLabeledRingWithTracking;
		};
		class LPolysResult:std::vector<LPolyImpl>
		{
			friend class FastAssociatedLabeledRingWithTracking;
		};

		FastAssociatedLabeledRingWithTracking(const RingZ2Slow&)
		{}
		bool QueueEmpty(const MultLPolysQueue& queue)
		{
			return queue.empty();
		}
		LPoly DequeueSigSmallest(MultLPolysQueue& queue);
		void PutInQueueExtendLabeledPolys(const PolysSet& in, MultLPolysQueue& queue);
		void FillWithTrivialSyzygiesOfNonMultElements(const MultLPolysQueue& queue, LPolysResult& to_fill);
		void ReduceCheckingSignatures(LPoly& poly, LPolysResult& reducers);
		
		ReconstructionInfo FieldAgnosticReconstructionInfo(const LPoly& poly);

		bool IsZero(const LPoly& poly);
		
		void Normalize(LPoly& poly);
		
		void InsertInResult(const LPoly& poly, LPolysResult& result);
		void ExtendRingWithMonomialToHelpReconstruct(const LPoly& poly, LPolysResult& reducers);
		void ExtendQueueBySpairPartsAndFilterUnneeded(const LPolysResult& left_parts, const LPoly& right_part, MultLPolysQueue& queue);
	};
	
	static std::unique_ptr<IOData<RingZ2Slow>> Create(const F4MPI::IOPolynomSet& in);

	static F4MPI::IOPolynomSet ConvertResult(std::unique_ptr<IOData<RingZ2Slow>>& result);
};

template <class MonomialMetadata, class Field>
class RingZ2Slow: public RingBase<MonomialMetadata, Field>, public RingZ2SlowBase
{
		RingZ2Slow(const RingZ2Slow& copy_from)
			:RingZ2SlowBase(copy_from)
		{
			
		}
		
		RingZ2Slow(const MonomialMetadata& monomial_metadata, const Field& field)
		{
			
		}
		
		RingZ2Slow& operator=(const RingZ2Slow& copy_from) = delete;
		
		bool ConstructAndInsertNormalized(const std::unique_ptr<const InPolysSetWithOrigMetadata>& prepared_input, const std::unique_ptr<const CrossRingInfo::MonomialListListWithTopInfo<MonomialMetadata>>& info, const std::unique_ptr<OutPolysSetForVariyingMetadata>& result)
		{
			return true;
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

