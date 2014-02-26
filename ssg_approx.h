#pragma once
#include <memory>
#include "ringbase.h"
#include "cross_ring_info.h"

template<class TRing, template <class Metadata> class TFastRing>
class ApproxSignatureGroebner
{
	typedef typename TRing::IOData IOData;
	typedef typename IOData::IOPolynomSet IOPolysSet;
	template <class TExactFastRing>
	struct ApproxSignatureGroebnerCalculator
	{
		ApproxSignatureGroebnerCalculator(IOData& io_data, TExactFastRing& fast_ring):
			in_(io_data.in_data), out_ring_(io_data.out_ring), out_(io_data.out_data), fast_ring_(fast_ring)
		{}

		const IOPolysSet& in_;
		TRing& out_ring_;
		std::unique_ptr<const IOPolysSet>& out_;
		TExactFastRing& fast_ring_;
		void CalcGB()
		{
			auto B = fast_ring_.PutInQueueExtendLabeledPolys(in_);
			auto R = fast_ring_.FillWithTrivialSyzygiesOfNonMultElements(B);
			auto reconstruction_basis = out_ring_.PrepareForReconstruction(in_);
			auto reconstructed_result = out_ring_.PrepareEmptyResult();
			
			while(!fast_ring_.QueueEmpty(B))
			{
				auto labeled_poly_to_reduce = fast_ring_.DequeueSigSmallest(B);
				for(;;)
				{
					fast_ring_.ReduceCheckingSignatures(labeled_poly_to_reduce, R);
					auto reconstruction_info = fast_ring_.FieldAgnosticReconstructionInfo(labeled_poly_to_reduce);
					if (out_ring_.ConstructAndInsertNormalized(reconstruction_basis, reconstruction_info, reconstructed_result))
					{
						//normalized polynomial was reconstructed in base ring
						break;
					}
					auto added_info = out_ring_.ExtendRingWithMonomialToHelpReconstruct(reconstruction_basis, reconstruction_info);
					fast_ring_.AddLabeledPolyBefore(added_info, R, labeled_poly_to_reduce);
				}
				if (!fast_ring_.IsZero(labeled_poly_to_reduce))
				{
					fast_ring_.Normalize(labeled_poly_to_reduce);
					fast_ring_.ExtendQueueBySpairPartsAndFilterUnneeded(R, labeled_poly_to_reduce, B);
				}
				fast_ring_.InsertInResult(labeled_poly_to_reduce, R);
			}
			out_ = out_ring_.ConvertResultToFixedMetadata(reconstructed_result);
		}
	};
public:
	static void Do(typename TRing::IOData& io_data)
	{
		typedef TFastRing<typename TRing::MonomialMetadata> TExactFastRing;
		TExactFastRing fast_ring {io_data.in_data.Metadata()};
		ApproxSignatureGroebnerCalculator<TExactFastRing> calculator {io_data, fast_ring};
		calculator.CalcGB();
	}
};

