#pragma once
#include <memory>
#include "ringbase.h"
#include "cross_ring_info.h"

template<class TRing, template <class Metadata> class TFastRing>
class ApproxSignatureGroebner
{
	typedef typename TRing::IOData IOData;
	typedef typename IOData::IOPolynomialSet IOPolysSet;
	template <class TExactFastRing>
	struct ApproxSignatureGroebnerCalculator
	{
		ApproxSignatureGroebnerCalculator(IOData& io_data, TExactFastRing& fast_ring):
			in_ring_(io_data.in_ring), in_(io_data.in), out_ring_(io_data.out_ring), out_(io_data.out), fast_ring_(fast_ring)
		{}

		const TRing& in_ring_;
		const IOPolysSet& in_;
		TRing& out_ring_;
		IOPolysSet& out_;
		TExactFastRing& fast_ring_;
		void CalcGB()
		{
			auto B = fast_ring_.PutInQueueExtendLabeledPolys(in_);
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
	};
public:
	static void Do(typename TRing::IOData& io_data)
	{
		io_data.in_ring.CopyTo(io_data.out_ring);
		typedef TFastRing<typename TRing::MonomialMetadata> TExactFastRing;
		TExactFastRing fast_ring {io_data.in.metadata()};
		ApproxSignatureGroebnerCalculator<TExactFastRing> calculator(io_data, fast_ring);
		calculator.CalcGB();
	}
};

