#include <gtest/gtest.h>
#include "cross_ring_info.h"

typedef CrossRingInfo::MonimailMetaData<CrossRingInfo::MonomialOrder::DegRevLex> DegRevLex;
typedef CrossRingInfo::PerVariableData V;
bool VarDegEqual(const V& v0, const V& v1)
{
	return (v0.index  == v1.index) && (v0.degree == v1.degree);
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
	DegRevLex order_ = []{DegRevLex order; order.var_count = 7; return order;}();
	CrossRingInfo::MonomialListList<DegRevLex> basis_info_{order_};
	CrossRingInfo::MonomialListListWithTopInfo<DegRevLex> poly_rec_info_{order_};
	const CrossRingInfo::MonomialListList<DegRevLex>& const_basis_info_ = basis_info_;
	const CrossRingInfo::MonomialListListWithTopInfo<DegRevLex>& const_poly_rec_info_ = poly_rec_info_;
};
typedef CrossRingInfoTest CrossRingInfoDeathTest;

TEST_F(CrossRingInfoDeathTest, AddToEmpty)
{
	EXPECT_DEATH(basis_info_.AddVariable(V(2,0)), "");
}

TEST_F(CrossRingInfoDeathTest, AddToEmptyWithTopInfo)
{
	poly_rec_info_.TopInfoAdditionDone();
	EXPECT_DEATH(poly_rec_info_.AddVariable(V(2,0)), "");
}

TEST_F(CrossRingInfoDeathTest, AddVarTooBig)
{
	poly_rec_info_.AddVariable(V(2,6));
	EXPECT_DEATH(poly_rec_info_.AddVariable(V(2,7)), "");
}

TEST_F(CrossRingInfoDeathTest, DoubleAddVar)
{
	poly_rec_info_.AddVariable(V(2,6));
	EXPECT_DEATH(poly_rec_info_.AddVariable(V(1,6)), "");
}

TEST_F(CrossRingInfoDeathTest, AddVarPastEnd)
{
	basis_info_.BeginPolynomialConstruction(3);
	basis_info_.MonomialAdditionDone();
	basis_info_.MonomialAdditionDone();
	basis_info_.MonomialAdditionDone();
	EXPECT_DEATH(basis_info_.AddVariable(V(1,6)), "");
}

TEST_F(CrossRingInfoDeathTest, AddVarPastEndWithTopInfo)
{
	poly_rec_info_.TopInfoAdditionDone();
	EXPECT_DEATH(poly_rec_info_.MonomialAdditionDone(), "");
}

TEST_F(CrossRingInfoDeathTest, ForgotTopInfoAddMon)
{
	EXPECT_DEATH(poly_rec_info_.MonomialAdditionDone(), "");
}

TEST_F(CrossRingInfoDeathTest, ForgotTopInfoAddPoly)
{
	EXPECT_DEATH(poly_rec_info_.BeginPolynomialConstruction(3), "");
}

TEST_F(CrossRingInfoDeathTest, ExtraTopInfo)
{
	poly_rec_info_.TopInfoAdditionDone();
	poly_rec_info_.BeginPolynomialConstruction(1);
	EXPECT_DEATH(poly_rec_info_.TopInfoAdditionDone(), "");
}

TEST_F(CrossRingInfoTest, Empty)
{
	poly_rec_info_.TopInfoAdditionDone();
	EXPECT_FALSE(const_poly_rec_info_.begin() != const_poly_rec_info_.end());
	EXPECT_FALSE(const_basis_info_.begin() != const_basis_info_.end());
}

TEST_F(CrossRingInfoTest, HasMetaData)
{
	EXPECT_EQ(&const_poly_rec_info_.MetaData(), &order_);
	EXPECT_EQ(&const_basis_info_.MetaData(), &order_);
}

TEST_F(CrossRingInfoTest, TopMon)
{
	poly_rec_info_.AddVariable(V(1,6));
	poly_rec_info_.AddVariable(V(99,5));
	poly_rec_info_.TopInfoAdditionDone();
	ExpecterContainerEqual(VarDegEqual)(const_poly_rec_info_.TopInfo(), ilist({V(99,5), V(1,6)}));
}

TEST_F(CrossRingInfoTest, TopMonEmpty)
{
	poly_rec_info_.TopInfoAdditionDone();
	ExpecterContainerEqual(VarDegEqual)(const_poly_rec_info_.TopInfo(), ilist<V>({}));
}


TEST_F(CrossRingInfoTest, 3x3)
{
	poly_rec_info_.TopInfoAdditionDone();
	poly_rec_info_.BeginPolynomialConstruction(3);
	poly_rec_info_.MonomialAdditionDone();
	poly_rec_info_.AddVariable(V(1,1));
	poly_rec_info_.AddVariable(V(2,3));
	poly_rec_info_.MonomialAdditionDone();
	poly_rec_info_.AddVariable(V(1,6));
	poly_rec_info_.AddVariable(V(99,5));
	poly_rec_info_.AddVariable(V(2,4));
	poly_rec_info_.MonomialAdditionDone();
	auto first_poly = *const_poly_rec_info_.begin();
	auto monomial_it = first_poly.begin();
	EXPECT_NE(monomial_it , first_poly.end());
	ExpecterContainerEqual(VarDegEqual)(*monomial_it,  ilist<V>({}));
	++monomial_it ;
	EXPECT_NE(monomial_it , first_poly.end());
	ExpecterContainerEqual(VarDegEqual)(*monomial_it,  ilist({V(1,1), V(2,3)}));
	++monomial_it ;
	EXPECT_NE(monomial_it , first_poly.end());
	ExpecterContainerEqual(VarDegEqual)(*monomial_it,  ilist({V(2,4), V(99,5), V(1,6)}));
	++monomial_it ;
	EXPECT_FALSE(monomial_it  != first_poly.end());

	basis_info_.BeginPolynomialConstruction(1);
	basis_info_.AddVariable(V(4,1));
	basis_info_.AddVariable(V(1,3));
	basis_info_.MonomialAdditionDone();
	basis_info_.BeginPolynomialConstruction(0);
	basis_info_.BeginPolynomialConstruction(3);
	basis_info_.AddVariable(V(1,1));
	basis_info_.AddVariable(V(2,3));
	basis_info_.MonomialAdditionDone();
	basis_info_.AddVariable(V(1,6));
	basis_info_.AddVariable(V(99,5));
	basis_info_.AddVariable(V(2,4));
	basis_info_.MonomialAdditionDone();
	basis_info_.MonomialAdditionDone();
	
	EXPECT_NE(const_basis_info_.begin(), const_basis_info_.end());

	ExpecterContainerEqual(ExpecterContainerEqual(ExpecterContainerEqual(VarDegEqual)))(const_basis_info_,
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
