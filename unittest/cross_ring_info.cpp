#include <gtest/gtest.h>
#include "cross_ring_info.h"

typedef CrossRingInfo::MonimailMetaData<CrossRingInfo::MonomialOrder::DegRevLex> DegRevLex;

TEST(CrossRingInfo, Empty)
{
	DegRevLex order;
	order.var_count = 7;
	CrossRingInfo::BasisElementReconstructionInfo<DegRevLex> poly_rec_info(order, 0);
	EXPECT_FALSE(poly_rec_info.begin() != poly_rec_info.end());
	CrossRingInfo::InputElementConstructionInfo<DegRevLex> input_rec_info(order);
}
