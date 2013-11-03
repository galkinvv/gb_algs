#include <gtest/gtest.h>
#include <functional>
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
			auto elements_comparator = [this](decltype(container1_cur) v1, decltype(*container2_cur) v2){return SavingResultCompareWithComparator(v1,v2, equal_);};
			EXPECT_PRED2(elements_comparator, container1_cur, *container2_cur);
			++container2_cur;
		}
		auto iterator_comparator = [this](decltype(container2_cur) v1, decltype(container2.end()) v2){return SavingResultCompare(v1,v2);};
		EXPECT_PRED2(iterator_comparator, container2_cur, container2.end());
		no_asserts_happen_ = true;
	}

	template <class Container1, class Container2>
	bool operator()(const Container1& container1, const Container2& container2)const
	{
		all_equal_ = true;
		no_asserts_happen_ = false;
		ExpectEqual(container1, container2);
		return all_equal_ && no_asserts_happen_;
	}
private:
	template <class V1, class V2>
	bool SavingResultCompare(V1& v1, V2& v2)const
	{
		return SavingResultCompareWithComparator(v1, v2);
	}

	template <class V1, class V2, class Comparator =  std::equal_to<V1>>
	bool SavingResultCompareWithComparator(V1& v1, V2& v2, const Comparator& comp = std::equal_to<V1>())const
	{
		bool result = comp(v1, v2);
		if (!result)
		{
			all_equal_ = false;
		}
		return result;
	}
	
	const SubComparator&  equal_;
	mutable bool all_equal_ = false;
	mutable bool no_asserts_happen_ = false;
};

template <class SubComparator> ContainerEqualExpect<SubComparator> ExpecterContainerEqual(const SubComparator& equal)
{
	return ContainerEqualExpect<SubComparator>(equal);
}

class CrossRingInfoTest: public ::testing::Test
{
public:
	const int initial_var_count_ = 7;
	DegRevLex order_ = [&]{DegRevLex order; order.var_count = initial_var_count_; return order;}();
	CrossRingInfo::MonomialListList<DegRevLex> basis_info_{order_};
	CrossRingInfo::MonomialListListWithTopInfo<DegRevLex> poly_rec_info_{order_};
	CrossRingInfo::SingleMonomial<DegRevLex> single_mon_{order_};
	const CrossRingInfo::MonomialListList<DegRevLex>& const_basis_info_ = basis_info_;
	const CrossRingInfo::MonomialListListWithTopInfo<DegRevLex>& const_poly_rec_info_ = poly_rec_info_;
	const	CrossRingInfo::SingleMonomial<DegRevLex>& const_single_mon_ = single_mon_;
	CrossRingInfoTest()
	{
		++order_.var_count;
	}
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
	poly_rec_info_.AddVariable(V(2,initial_var_count_-1));
	EXPECT_DEATH(poly_rec_info_.AddVariable(V(2,initial_var_count_)), "");
}

TEST_F(CrossRingInfoDeathTest, DoubleAddVar)
{
	poly_rec_info_.AddVariable(V(2,6));
	EXPECT_DEATH(poly_rec_info_.AddVariable(V(1,6)), "");
}

TEST_F(CrossRingInfoDeathTest, DoubleAddVarToSingle)
{
	single_mon_.AddVariable(V(2,0));
	EXPECT_DEATH(single_mon_.AddVariable(V(1,0)), "");
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
	EXPECT_FALSE(const_single_mon_.begin() != const_single_mon_.end());
}

TEST_F(CrossRingInfoTest, KeepsMetaData)
{
	EXPECT_EQ(const_poly_rec_info_.MetaData().var_count, initial_var_count_);
	EXPECT_EQ(CopyValue(const_poly_rec_info_.MetaData().order), CopyValue(DegRevLex::order));

	EXPECT_EQ(const_basis_info_.MetaData().var_count, initial_var_count_);
	EXPECT_EQ(const_single_mon_.MetaData().var_count, initial_var_count_);
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
	
	single_mon_.AddVariable(V(3,3));
	single_mon_.AddVariable(V(6,1));
	single_mon_.AddVariable(V(5,5));
	ExpecterContainerEqual(VarDegEqual)(const_single_mon_, ilist({V(6,1), V(3,3), V(5,5)}));
}
