#include <map>
#include <set>
#include <vector>
#include <memory>
#include "mock_base.h"
namespace Mock
{
	class NoCopy
	{
		void operator=(const NoCopy&);
		NoCopy(const NoCopy&);
	protected:
		NoCopy(){}
	};

	class Ring: NoCopy
	{
		typedef std::map<char,int> Monomial;
		typedef std::vector<Monomial> Polynomial;
		
		Ring(){}
		friend std::unique_ptr<Ring> CreateRing();
		
		struct ReconstructionInfoImpl: std::vector<Monomial>{
			Monomial top;
		};
	public:
		class FastAssociatedLabeledRingWithTracking;
		class ReconstructionInfo: ReconstructionInfoImpl{
			friend class Ring;
			friend class FastAssociatedLabeledRingWithTracking;
		};

		void CopyTo(Ring& to)const{}

		struct PolysSet: private std::vector<Polynomial>
		{
			FRIEND_FOR_TEST
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
			
			FastAssociatedLabeledRingWithTracking(const Ring&)
			{}
			bool QueueEmpty(const MultLPolysQueue& queue)
			{
				return queue.empty();
			}
			LPoly DequeueSigSmallest(MultLPolysQueue& queue);
			void PutInQueueExtendLabeledPolys(const PolysSet& in, MultLPolysQueue& queue)
			{
				//TODO
			}
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
	};
	std::unique_ptr<Ring> CreateRing();	
}