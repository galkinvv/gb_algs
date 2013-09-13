#pragma once
namespace
{
	template<class TRing>
	class ApproxSigantureGroebner
	{
		typedef typename TRing::PolysSet PolysSet;
		typedef typename TRing::FastAssociatedLabeledRingWithTracking FastRing;
		ApproxSigantureGroebner(const TRing& in_ring, const PolysSet& in, TRing& out_ring, PolysSet& out, FastRing& fast_ring):
			in_ring_(in_ring), in_(in), out_ring_(out_ring), out_(out), fast_ring_(fast_ring)
		{}

		typename FastRing::LPolysResult R_;
		typename FastRing::MultLPolysQueue B_;
		const TRing& in_ring_;
		const PolysSet& in_;
		TRing& out_ring_;
		PolysSet& out_;
		FastRing& fast_ring_;
		void CalcGB()
		{
			fast_ring_.PutInQueueExtendLabeledPolys(in_, B_);
			fast_ring_.FillWithTrivialSyzygiesOfNonMultElements(B_, R_);
			while(!fast_ring_.QueueEmpty(B_))
			{
				auto labeled_poly_to_reduce = fast_ring_.DequeueSigSmallest(B_);
				bool got_zero = false;
				for(;;)
				{
					fast_ring_.ReduceCheckingSignatures(labeled_poly_to_reduce, R_);
					auto reconstruction_info = fast_ring_.FieldAgnosticReconstructionInfo(labeled_poly_to_reduce);
					if (out_ring_.ConstructAndInsertNormalized(in_, reconstruction_info, out_))
					{
						//normalized polynomial was reconstructed in base ring
						break;
					}
					fast_ring_.ExtendRingWithMonomialToHelpReconstruct(labeled_poly_to_reduce, R_);
				}
				if (!fast_ring_.IsZero(labeled_poly_to_reduce))
				{
					fast_ring_.Normalize(labeled_poly_to_reduce);
					fast_ring_.ExtendQueueBySpairPartsAndFilterUnneeded(R_, labeled_poly_to_reduce, B_);
				}
				fast_ring_.InsertInResult(labeled_poly_to_reduce, R_);
			}
		}
	public:
		static void Do(const TRing& in_ring, const PolysSet& in, TRing& out_ring, PolysSet& out)
		{
			in_ring.CopyTo(out_ring);
			FastRing fast_ring(out_ring);
			ApproxSigantureGroebner calculator(in_ring, in, out_ring, out, fast_ring);
			calculator.CalcGB();
		}
	};
}
