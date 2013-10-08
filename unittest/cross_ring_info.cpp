#include <gtest/gtest.h>
#include "cross_ring_info.h"

typedef CrossRingInfo::MonimailMetaData<CrossRingInfo::MonomialOrder::DegRevLex> DegRevLex;
typedef CrossRingInfo::PerVariableData V;
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
	poly_rec_info.AddVariable(V(99,2));
	poly_rec_info.AddVariable(V(2,3));
	poly_rec_info.MonomialAdditionDone();
	poly_rec_info.AddVariable(V(1,4));
	poly_rec_info.AddVariable(V(99,5));
	poly_rec_info.AddVariable(V(2,6));
	poly_rec_info.MonomialAdditionDone();
	auto b = poly_rec_info.begin();
	EXPECT_NE(b, poly_rec_info.end());
	++b;
	EXPECT_NE(b, poly_rec_info.end());
	++b;
	EXPECT_NE(b, poly_rec_info.end());
	++b;
	EXPECT_FALSE(b != poly_rec_info.end());
	CrossRingInfo::InputElementsConstructionInfo<DegRevLex> input_rec_info(order);
	EXPECT_FALSE(input_rec_info.begin() != input_rec_info.end());
}
