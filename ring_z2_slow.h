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
	struct Monomial : std::map<char,int>{};
	friend bool operator<(const Monomial&, const Monomial&); //undefined
	struct Polynomial : std::vector<Monomial>{};
	friend bool operator<(const Polynomial&, const Polynomial&); //undefined

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
	class ReconstructionInfo: ReconstructionInfoImpl{
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
	bool ConstructAndInsertNormalized(const PolysSet& in, const ReconstructionInfo& info, PolysSet& out)
	{
		//TODO
		return true;
	}

	static std::unique_ptr<IOData<RingZ2Slow>> Create(const F4MPI::IOPolynomSet& in);

	static F4MPI::IOPolynomSet ConvertResult(std::unique_ptr<IOData<RingZ2Slow>>& result);
};
