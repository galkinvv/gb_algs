#include <gtest/gtest.h>
#include "cross_ring_info.h"

typedef CrossRingInfo::MonimailMetaData<CrossRingInfo::MonomialOrder::DegRevLex> DegRevLex;
typedef CrossRingInfo::PerVariableData V;
bool VarDegEqual(const V& v0, const V& v1)
{
	return (v0.index  == v1.index) && (v0.degree == v1.degree);
}

template <class Container1, class Container2, class Comparator>
void ExpectEqualToContainer(const Container1& container1, const Container2& container2, const Comparator& equal)
{
}

template <class SubComparator>
struct ContainerEqualExpect
{
	explicit ContainerEqualExpect(const SubComparator& equal)
		:equal_(equal)
	{}
	
	template <class Container1, class Container2>
	void ExpectEqual(const Container1& container1, const Container2& container2)const
	{
		auto container2_cur = container2.begin();
		for(auto container1_cur:container1)
		{
			ASSERT_NE(container2_cur, container2.end());
			EXPECT_PRED2(equal_, container1_cur, *container2_cur);
			++container2_cur;
		}
		EXPECT_EQ(container2_cur, container2.end());
	}

	template <class Container1, class Container2>
	bool operator()(const Container1& container1, const Container2& container2)const
	{
		ExpectEqual(container1, container2);
		return true;
	}
private:
	const SubComparator& equal_;
};

template <class SubComparator> ContainerEqualExpect<SubComparator> ExpecterContainerEqual(const SubComparator& equal)
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
	ExpecterContainerEqual(VarDegEqual)(*monomial_it,  ilist<V>({}));
	++monomial_it ;
	EXPECT_NE(monomial_it , poly_rec_info.end());
	ExpecterContainerEqual(VarDegEqual)(*monomial_it,  ilist({V(1,1), V(2,3)}));
	++monomial_it ;
	EXPECT_NE(monomial_it , poly_rec_info.end());
	ExpecterContainerEqual(VarDegEqual)(*monomial_it,  ilist({V(2,4), V(99,5), V(1,6)}));
	++monomial_it ;
	EXPECT_FALSE(monomial_it  != poly_rec_info.end());
	CrossRingInfo::InputElementsConstructionInfo<DegRevLex> input_rec_info(order_);
	input_rec_info.BeginPolynomialConstruction(1);
	input_rec_info.AddVariable(V(4,1));
	input_rec_info.AddVariable(V(1,3));
	input_rec_info.MonomialAdditionDone();
	input_rec_info.BeginPolynomialConstruction(0);
	input_rec_info.BeginPolynomialConstruction(3);
	input_rec_info.AddVariable(V(1,1));
	input_rec_info.AddVariable(V(2,3));
	input_rec_info.MonomialAdditionDone();
	input_rec_info.AddVariable(V(1,6));
	input_rec_info.AddVariable(V(99,5));
	input_rec_info.AddVariable(V(2,4));
	input_rec_info.MonomialAdditionDone();
	input_rec_info.MonomialAdditionDone();
	
	EXPECT_NE(input_rec_info.begin(), input_rec_info.end());

	ExpecterContainerEqual(ExpecterContainerEqual(ExpecterContainerEqual(VarDegEqual)))(input_rec_info,
		ilist({
			ilist({
				ilist({V(4,1), V(1,3)}),
			}),
			ilist<std::initializer_list<V>>({}),
			ilist({
				ilist({V(1,1), V(2,3)}),
				ilist({V(2,4), V(99,5), V(1,6)}),
				ilist<V>({}),
			}),
		})
	);
}
