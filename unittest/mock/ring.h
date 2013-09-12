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

	class BaseRing
	{
	protected:
		typedef std::map<char,int> Monomial;
		typedef std::vector<Monomial> Polynomial;
	};

	struct Ring:BaseRing, NoCopy
	{
		struct FastAssociatedLabeledRingWithTracking : NoCopy
		{
			FastAssociatedLabeledRingWithTracking(FastAssociatedLabeledRingWithTracking &&)
			{}
			static FastAssociatedLabeledRingWithTracking Create()
			{
				return FastAssociatedLabeledRingWithTracking();
			}
		private:
			FastAssociatedLabeledRingWithTracking (){}
		};
		
		struct PolysSet: private std::vector<Polynomial>
		{
			void expect_empty()
			{
				EXPECT_EQ(this->size(), 0);
			}
		};

		friend std::unique_ptr<Ring> CreateRing();
	private:
		Ring(){}
	};

	std::unique_ptr<Ring> CreateRing()
	{
		std::unique_ptr<Ring> result;
		result.reset(new Ring());
		return result;
	}
	
	
}