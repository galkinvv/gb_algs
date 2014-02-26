#include "test_base.h"
#include <functional>
#include "cross_ring_info.h"

typedef CrossRingInfo::MonomialMetadata<CrossRingInfo::MonomialOrder::DegRevLex> DegRevLex;
typedef CrossRingInfo::PerVariableData V;

bool VarDegEqual(const V& v0, const V& v1)
{
	return (v0.index  == v1.index) && (v0.degree == v1.degree);
}

class StrangeFieldValueInitializer
{
};

class StrangeFieldInitializer
{
};

class StrangeFieldValue
{
	void operator=(const StrangeFieldValue&) = delete;
	StrangeFieldValue(const StrangeFieldValue&) = delete;

public:
	StrangeFieldValue(StrangeFieldValue&&) = default;
	StrangeFieldValue(const StrangeFieldValueInitializer&){}
};

inline std::ostream& operator << (std::ostream& s, const StrangeFieldValue& data)
{
	IgnoreIfUnused(data);
	return s << "qu";
}

class StrangeField
{
	void operator=(const StrangeField&) = delete;
	StrangeField(StrangeField&&) = delete;

  public:
	typedef StrangeFieldValue Value;
	StrangeField(const StrangeField&) = default;
	StrangeField(const StrangeFieldInitializer&){}
};

struct IntField
{
	typedef int Value;
	int field_option;
};
class CrossRingInfoTest: public ::testing::Test
{
public:
	const int initial_var_count_ = 7;
	const int initial_field_option_ = 7;
	const int added_poly_index_ = 42;
	StrangeField strange_field_{StrangeFieldInitializer()};
	IntField int_field_ = [&]{IntField int_field; int_field.field_option = initial_field_option_; return int_field;}();
	DegRevLex order_ = [&]{DegRevLex order; order.var_count = initial_var_count_; return order;}();
	CrossRingInfo::MonomialListListWithCoef<DegRevLex, StrangeField> strange_field_basis_info_{order_, strange_field_};
	CrossRingInfo::MonomialListListWithCoef<DegRevLex, IntField> int_field_basis_info_{order_, int_field_};
	CrossRingInfo::MonomialListListWithTopInfo<DegRevLex> poly_rec_info_{order_};
	CrossRingInfo::AddedVarInfo<DegRevLex> single_mon_{order_, order_.var_count, added_poly_index_};
	const CrossRingInfo::MonomialListListWithCoef<DegRevLex, StrangeField>& const_basis_info_ = strange_field_basis_info_;
	const CrossRingInfo::MonomialListListWithTopInfo<DegRevLex>& const_poly_rec_info_ = poly_rec_info_;
	const CrossRingInfo::AddedVarInfo<DegRevLex>& const_single_mon_ = single_mon_;
	CrossRingInfoTest()
	{
		++order_.var_count;
	}
};
typedef CrossRingInfoTest CrossRingInfoDeathTest;

TEST_F(CrossRingInfoDeathTest, AddToEmpty)
{
	EXPECT_DEATH(strange_field_basis_info_.AddVariable(V(2,0)), "");
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
	strange_field_basis_info_.BeginPolynomialConstruction(3);
	strange_field_basis_info_.MonomialAdditionDone(StrangeFieldValueInitializer());
	strange_field_basis_info_.MonomialAdditionDone(StrangeFieldValueInitializer());
	strange_field_basis_info_.MonomialAdditionDone(StrangeFieldValueInitializer());
	EXPECT_DEATH(strange_field_basis_info_.AddVariable(V(1,6)), "");
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

TEST_F(CrossRingInfoDeathTest, MapGetPostEnd)
{
	CrossRingInfo::VariableMapping<DegRevLex> mapping{order_, 1};
	mapping.MonomialAdditionDone();
	EXPECT_DEATH(mapping[2], "");
}

TEST_F(CrossRingInfoTest, Empty)
{
	poly_rec_info_.TopInfoAdditionDone();
	EXPECT_FALSE(const_poly_rec_info_.begin() != const_poly_rec_info_.end());
	EXPECT_FALSE(const_basis_info_.begin() != const_basis_info_.end());
	EXPECT_FALSE(const_single_mon_.begin() != const_single_mon_.end());
}

TEST_F(CrossRingInfoTest, KeepsMetadata)
{
	EXPECT_EQ(const_poly_rec_info_.Metadata().var_count, initial_var_count_);
	EXPECT_EQ(CopyValue(const_poly_rec_info_.Metadata().order), CopyValue(DegRevLex::order));

	EXPECT_EQ(const_basis_info_.Metadata().var_count, initial_var_count_);
	EXPECT_EQ(const_single_mon_.Metadata().var_count, initial_var_count_);
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

	strange_field_basis_info_.BeginPolynomialConstruction(1);
	strange_field_basis_info_.AddVariable(V(4,1));
	strange_field_basis_info_.AddVariable(V(1,3));
	strange_field_basis_info_.MonomialAdditionDone(StrangeFieldValueInitializer());
	strange_field_basis_info_.BeginPolynomialConstruction(0);
	strange_field_basis_info_.BeginPolynomialConstruction(3);
	strange_field_basis_info_.AddVariable(V(1,1));
	strange_field_basis_info_.AddVariable(V(2,3));
	strange_field_basis_info_.MonomialAdditionDone(StrangeFieldValueInitializer());
	strange_field_basis_info_.AddVariable(V(1,6));
	strange_field_basis_info_.AddVariable(V(99,5));
	strange_field_basis_info_.AddVariable(V(2,4));
	strange_field_basis_info_.MonomialAdditionDone(StrangeFieldValueInitializer());
	strange_field_basis_info_.MonomialAdditionDone(StrangeFieldValueInitializer());
	
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

TEST_F(CrossRingInfoTest, SimpleAdditionWithCoef)
{
	strange_field_basis_info_.BeginPolynomialConstruction(1);
	strange_field_basis_info_.MonomialAdditionDone(StrangeFieldValueInitializer());
	int_field_basis_info_.BeginPolynomialConstruction(1);
	int_field_basis_info_.MonomialAdditionDone(42);
}
