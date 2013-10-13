#include <gtest/gtest.h>
#include "cross_ring_info.h"

typedef CrossRingInfo::MonimailMetaData<CrossRingInfo::MonomialOrder::DegRevLex> DegRevLex;
typedef CrossRingInfo::PerVariableData V;
bool VarDegEqual(const V& v0, const V& v1)
{
	return (v0.index  == v1.index) && (v0.degree == v1.degree);
}
template <class Iterator, class Container, class Comparator>
void ExpectIterEqualToContainer(const Iterator& begin, const Iterator& end, const Container& container, const Comparator& equal)
{
	SCOPED_TRACE(container.size());
	Iterator curr_it = begin;
	auto container_cur = container.begin();
	for(;;)
	{
		bool it_ends = !(curr_it != end);
		bool cont_ends = !(container_cur != container.end());
		EXPECT_EQ(it_ends, cont_ends);
		if (it_ends || cont_ends)
		{
			break;
		}
		EXPECT_PRED2(equal, *curr_it, *container_cur);
		++curr_it;
		++container_cur;
	}
}

class CrossRingInfoTest: public ::testing::Test
{
public:
	DegRevLex order_;
	CrossRingInfoTest()
	{
		order_.var_count = 7;
	}
};
typedef CrossRingInfoTest CrossRingInfoDeathTest;

TEST_F(CrossRingInfoDeathTest, AddToEmpty)
{
	CrossRingInfo::BasisElementReconstructionInfo<DegRevLex> poly_rec_info(order_, 0);
	EXPECT_DEATH(poly_rec_info.AddVariable(V(2,0)), "");
}

TEST_F(CrossRingInfoDeathTest, AddVarTooBig)
{
	CrossRingInfo::BasisElementReconstructionInfo<DegRevLex> poly_rec_info(order_, 1);
	poly_rec_info.AddVariable(V(2,6));
	EXPECT_DEATH(poly_rec_info.AddVariable(V(2,7)), "");
}

TEST_F(CrossRingInfoDeathTest, DoubleAddVar)
{
	CrossRingInfo::BasisElementReconstructionInfo<DegRevLex> poly_rec_info(order_, 1);
	poly_rec_info.AddVariable(V(2,6));
	EXPECT_DEATH(poly_rec_info.AddVariable(V(1,6)), "");
}

TEST_F(CrossRingInfoDeathTest, AddVarPastEnd)
{
	CrossRingInfo::BasisElementReconstructionInfo<DegRevLex> poly_rec_info(order_, 1);
	poly_rec_info.MonomialAdditionDone();
	EXPECT_DEATH(poly_rec_info.AddVariable(V(1,6)), "");
}

TEST_F(CrossRingInfoTest, Empty)
{
	CrossRingInfo::BasisElementReconstructionInfo<DegRevLex> poly_rec_info(order_, 0);
	EXPECT_FALSE(poly_rec_info.begin() != poly_rec_info.end());
	CrossRingInfo::InputElementsConstructionInfo<DegRevLex> input_rec_info(order_);
	EXPECT_FALSE(input_rec_info.begin() != input_rec_info.end());
}
TEST_F(CrossRingInfoTest, 3x3)
{
	CrossRingInfo::BasisElementReconstructionInfo<DegRevLex> poly_rec_info(order_, 3);
	poly_rec_info.AddVariable(V(1,6));
	poly_rec_info.AddVariable(V(99,5));
	poly_rec_info.AddVariable(V(2,0));
	poly_rec_info.MonomialAdditionDone();
	poly_rec_info.AddVariable(V(1,1));
	poly_rec_info.AddVariable(V(2,3));
	poly_rec_info.MonomialAdditionDone();
	poly_rec_info.AddVariable(V(1,6));
	poly_rec_info.AddVariable(V(99,5));
	poly_rec_info.AddVariable(V(2,4));
	poly_rec_info.MonomialAdditionDone();
	auto monomial_it = poly_rec_info.begin();
	EXPECT_NE(monomial_it , poly_rec_info.end());
	ExpectIterEqualToContainer(monomial_it ->begin(), monomial_it ->end(),  ilist({V(2,0), V(99,5), V(1,6)}), VarDegEqual);
	++monomial_it ;
	EXPECT_NE(monomial_it , poly_rec_info.end());
	ExpectIterEqualToContainer(monomial_it ->begin(), monomial_it ->end(),  ilist({V(1,1), V(2,3)}), VarDegEqual);
	++monomial_it ;
	EXPECT_NE(monomial_it , poly_rec_info.end());
	ExpectIterEqualToContainer(monomial_it ->begin(), monomial_it ->end(),  ilist({V(2,4), V(99,5), V(1,6)}), VarDegEqual);
	++monomial_it ;
	EXPECT_FALSE(monomial_it  != poly_rec_info.end());
}
