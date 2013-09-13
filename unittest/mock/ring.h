#include <map>
#include <set>
#include <vector>
#include <memory>
#include <gtest/gtest.h>
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
		struct Monomial: std::map<char,int>{};
		struct Polynomial: std::vector<Monomial>{};
		
		Ring(){}
		friend std::unique_ptr<Ring> CreateRing();
		
	public:
		class FastAssociatedLabeledRingWithTracking;
		class ReconstructionInfo:std::vector<Monomial>{
			Monomial top;
			friend class Ring;
			friend class FastAssociatedLabeledRingWithTracking;
		};

		void CopyTo(Ring& to)const{}

		struct PolysSet: private std::vector<Polynomial>
		{
			void expect_empty()
			{
				EXPECT_EQ(this->size(), 0);
			}
		};
		class FastAssociatedLabeledRingWithTracking : NoCopy
		{
			class FastPoly: std::vector<Monomial>
			{
				friend class FastAssociatedLabeledRingWithTracking;
			};
			class LPoly
			{
				FastPoly value;
				FastPoly reconstruction_info;
				double sig_index;
				Monomial sig_mon;
				friend class FastAssociatedLabeledRingWithTracking;
			};
		private:
			struct MultLPoly
			{
				LPoly poly;
				Monomial mul_by;
			};
		public:
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
			LPoly DequeueSigSmallest(MultLPolysQueue& queue)
			{
				//TODO
			}
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

	std::unique_ptr<Ring> CreateRing()
	{
		std::unique_ptr<Ring> result;
		result.reset(new Ring());
		return result;
	}
	
	
}