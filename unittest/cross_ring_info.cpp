#include <gtest/gtest.h>
#include "cross_ring_info.h"

typedef CrossRingInfo::MonimailMetaData<CrossRingInfo::MonomialOrder::DegRevLex> DegRevLex;
typedef CrossRingInfo::PerVariableData V;
bool VarDegEqual(const V& v0, const V& v1)
{
	return (v0.index  == v1.index) && (v0.degree == v1.degree);
}

template <class ContainerWithoutSize, class Container, class Comparator>
void ExpectEqualToContainer(const ContainerWithoutSize& cont_no_size, const Container& container, const Comparator& equal)
{
	SCOPED_TRACE(container.size());
	auto container_cur = container.begin();
	for(auto item:cont_no_size)
	{
		ASSERT_NE(container_cur, container.end());
		EXPECT_PRED2(equal, item, *container_cur);
		++container_cur;
	}
	EXPECT_EQ(container_cur, container.end());
}

template <class SubComparator>
struct ContainerEqualExpect
{
	explicit ContainerEqualExpect(const SubComparator& equal)
		:equal_(equal)
	{}
	
	template <class ContainerWithoutSize, class Container>
	bool operator()(const ContainerWithoutSize& cont_no_size, const Container& container)const
	{
		ExpectEqualToContainer(cont_no_size, container, equal_);
		return true;
	}
private:
	const SubComparator& equal_;
};

template <class SubComparator> ContainerEqualExpect<SubComparator> CreateContainerComparator(const SubComparator& equal)
{
	return ContainerEqualExpect<SubComparator>(equal);
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
	ExpectEqualToContainer(*monomial_it,  ilist<V>({}), VarDegEqual);
	++monomial_it ;
	EXPECT_NE(monomial_it , poly_rec_info.end());
	ExpectEqualToContainer(*monomial_it,  ilist({V(1,1), V(2,3)}), VarDegEqual);
	++monomial_it ;
	EXPECT_NE(monomial_it , poly_rec_info.end());
	ExpectEqualToContainer(*monomial_it,  ilist({V(2,4), V(99,5), V(1,6)}), VarDegEqual);
	++monomial_it ;
	EXPECT_FALSE(monomial_it  != poly_rec_info.end());
	CrossRingInfo::InputElementsConstructionInfo<DegRevLex> input_rec_info(order_);
	input_rec_info.BeginPolynomialConstruction(1);
	input_rec_info.AddVariable(V(4,1));
	input_rec_info.AddVariable(V(1,3));
	input_rec_info.MonomialAdditionDone();
	input_rec_info.BeginPolynomialConstruction(3);
	input_rec_info.AddVariable(V(1,1));
	input_rec_info.AddVariable(V(2,3));
	input_rec_info.MonomialAdditionDone();
	input_rec_info.AddVariable(V(1,6));
	input_rec_info.AddVariable(V(99,5));
	input_rec_info.AddVariable(V(2,4));
	input_rec_info.MonomialAdditionDone();
	input_rec_info.MonomialAdditionDone();
	input_rec_info.BeginPolynomialConstruction(0);
	

	ExpectEqualToContainer(input_rec_info,
		ilist({
			ilist({
				ilist({V(2,4), V(99,5), V(1,6)})
			})
		}), 
		CreateContainerComparator(CreateContainerComparator(VarDegEqual))
	);
	for (auto input_poly:input_rec_info)
	{
		for(auto mon:input_poly)
		{
			ExpectEqualToContainer(mon,  ilist({V(2,4), V(99,5), V(1,6)}), VarDegEqual);
		}
	}
	EXPECT_FALSE(input_rec_info.begin() != input_rec_info.end());

}
