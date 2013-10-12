#include <gtest/gtest.h>
#include "cross_ring_info.h"

typedef CrossRingInfo::MonimailMetaData<CrossRingInfo::MonomialOrder::DegRevLex> DegRevLex;
typedef CrossRingInfo::PerVariableData V;
bool VarDegEqual(const V& v0, const V& v1)
{
	return (v0.index  == v1.index) && (v0.degree == v1.degree);
}
template <class Iterator, class Container, class Comparator>
bool IterEqualToContainer(const Iterator& begin, const Iterator& end, const Container& container, const Comparator& cmp)
{
	Iterator curr_it = begin;
	auto container_cur = container.begin();
	for(;;)
	{
		bool it_ends = !(curr_it != end);
		bool cont_ends = !(container_cur != container.end());
		if (it_ends && cont_ends)
		{
			return true;
		}
		if (it_ends != cont_ends)
		{
			return false;
		}
		if (!cmp(*curr_it, *container_cur))
		{
			return false;
		}
		++curr_it;
		++container_cur;
	}
}

TEST(CrossRingInfo, Empty)
{
	DegRevLex order;
	order.var_count = 7;
	CrossRingInfo::BasisElementReconstructionInfo<DegRevLex> poly_rec_info(order, 0);
	EXPECT_FALSE(poly_rec_info.begin() != poly_rec_info.end());
	CrossRingInfo::InputElementsConstructionInfo<DegRevLex> input_rec_info(order);
	EXPECT_FALSE(input_rec_info.begin() != input_rec_info.end());
}
TEST(CrossRingInfo, 3x3)
{
	DegRevLex order;
	order.var_count = 7;
	CrossRingInfo::BasisElementReconstructionInfo<DegRevLex> poly_rec_info(order, 3);
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
	auto b = poly_rec_info.begin();
	EXPECT_NE(b, poly_rec_info.end());
	auto expected1 = {V(2,0), V(1,6), V(99,5)};
	EXPECT_TRUE(IterEqualToContainer(b->begin(), b->end(),  expected1, VarDegEqual));
	++b;
	EXPECT_NE(b, poly_rec_info.end());
	++b;
	EXPECT_NE(b, poly_rec_info.end());
	++b;
	EXPECT_FALSE(b != poly_rec_info.end());
	CrossRingInfo::InputElementsConstructionInfo<DegRevLex> input_rec_info(order);
	EXPECT_FALSE(input_rec_info.begin() != input_rec_info.end());
}
