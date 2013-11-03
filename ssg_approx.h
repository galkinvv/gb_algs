#pragma once
#include <memory>
#include "ringbase.h"
#include "cross_ring_info.h"

template<class TRing, class TFastRing>
class ApproxSignatureGroebner
{
	typedef typename TRing::PolysSet PolysSet;
	ApproxSignatureGroebner(IOData<TRing>& io_data, TFastRing& fast_ring):
		in_ring_(io_data.in_ring), in_(io_data.in), out_ring_(io_data.out_ring), out_(io_data.out), fast_ring_(fast_ring)
	{}

	const TRing& in_ring_;
	const PolysSet& in_;
	TRing& out_ring_;
	PolysSet& out_;
	TFastRing& fast_ring_;
	void CalcGB()
	{
		auto input_info = in_ring_.GetCrossRingInfoForInput(in_);
		auto B = fast_ring_.PutInQueueExtendLabeledPolys(input_info);
		auto R = fast_ring_.FillWithTrivialSyzygiesOfNonMultElements(B);
		
		while(!fast_ring_.QueueEmpty(B))
		{
			auto labeled_poly_to_reduce = fast_ring_.DequeueSigSmallest(B);
			for(;;)
			{
				fast_ring_.ReduceCheckingSignatures(labeled_poly_to_reduce, R);
				auto reconstruction_info = fast_ring_.FieldAgnosticReconstructionInfo(labeled_poly_to_reduce);
				if (out_ring_.ConstructAndInsertNormalized(in_, reconstruction_info, out_))
				{
					//normalized polynomial was reconstructed in base ring
					break;
				}
				auto added_monomial = fast_ring_.ExtendRingWithMonomialToHelpReconstruct(labeled_poly_to_reduce, R);
				out_ring_.ExtendWithMonomial(added_monomial);
			}
			if (!fast_ring_.IsZero(labeled_poly_to_reduce))
			{
				fast_ring_.Normalize(labeled_poly_to_reduce);
				fast_ring_.ExtendQueueBySpairPartsAndFilterUnneeded(R, labeled_poly_to_reduce, B);
			}
			fast_ring_.InsertInResult(labeled_poly_to_reduce, R);
		}
	}
public:
	static void Do(IOData<TRing>& io_data)
	{
		io_data.in_ring.CopyTo(io_data.out_ring);
		TFastRing fast_ring;
		ApproxSignatureGroebner calculator(io_data, fast_ring);
		calculator.CalcGB();
	}
};

