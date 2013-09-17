#pragma once
#include <map>
#include <set>
#include <vector>
#include <memory>
#include "algs.h"
#include "ringbase.h"
namespace F4MPI
{
	class IOPolynomSet;
}
class RingZ2Slow: NoCopy
{
	//Z_2 ring with degrevlex oredr on variables
	typedef std::map<char,int> Monomial;
	typedef std::vector<Monomial> Polynomial;

	RingZ2Slow(){}
	friend class IOData<RingZ2Slow>;

	struct ReconstructionInfoImpl: std::vector<Monomial>{
		Monomial top;
	};
public:
	class FastAssociatedLabeledRingWithTracking;
	class ReconstructionInfo: ReconstructionInfoImpl{
		friend class RingZ2Slow;
		friend class FastAssociatedLabeledRingWithTracking;
	};

	void CopyTo(RingZ2Slow& to)const{}

	struct PolysSet: private std::vector<Polynomial>
	{
		friend class RingZ2Slow;
	};

	class FastAssociatedLabeledRingWithTracking : NoCopy
	{
		typedef std::vector<Monomial> FastPoly;
		struct LPolyImpl
		{
			FastPoly value;
			FastPoly reconstruction_info;
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
		class MultLPolysQueue:std::vector<MultLPoly>{
			friend class FastAssociatedLabeledRingWithTracking;
		};
		class LPolysResult:std::vector<LPoly>{};

		FastAssociatedLabeledRingWithTracking(const RingZ2Slow&)
		{}
		bool QueueEmpty(const MultLPolysQueue& queue)
		{
			return queue.empty();
		}
		LPoly DequeueSigSmallest(MultLPolysQueue& queue);
		void PutInQueueExtendLabeledPolys(const PolysSet& in, MultLPolysQueue& queue);
		void FillWithTrivialSyzygiesOfNonMultElements(const MultLPolysQueue& queue, LPolysResult& to_fill)
		{
			//TODO
		}
		void ReduceCheckingSignatures(LPoly& poly, LPolysResult& reducers)
		{
			//TODO
		}
		ReconstructionInfo FieldAgnosticReconstructionInfo(const LPoly& poly)
		{
			ReconstructionInfo result;
			static_cast<std::vector<Monomial>&>(result) = poly.reconstruction_info;
			//TODO: fill top
			return result;
		}
		bool IsZero(const LPoly& poly)
		{
			return poly.value.empty();
		}
		void Normalize(LPoly& poly)
		{
			//TODO
		}
		void ExtendRingWithMonomialToHelpReconstruct(const LPoly& poly, LPolysResult& reducers)
		{
			//TODO
		}
		void ExtendQueueBySpairPartsAndFilterUnneeded(const LPolysResult& left_parts, const LPoly& right_part, MultLPolysQueue& queue)
		{
			//TODO
		}
		void InsertInResult(const LPoly& poly, LPolysResult& result)
		{
			//TODO
		}

	};
	bool ConstructAndInsertNormalized(const PolysSet& in, const ReconstructionInfo& info, PolysSet& out)
	{
		//TODO
	}

	static std::unique_ptr<IOData<RingZ2Slow>> Create(const F4MPI::IOPolynomSet& in);

	static F4MPI::IOPolynomSet ConvertResult(std::unique_ptr<IOData<RingZ2Slow>>& result);
};
