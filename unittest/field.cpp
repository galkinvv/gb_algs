#include <limits>
#include <cstdint>
#include "test_base.h"
#include "z_ring.h"
#include "finite_field.h"

template <class Field>
struct CreateFieldForTest
{
	static Field Create()
	{
		return Field();
	}
};

template <class BaseForFiniteField>
struct FiniteFieldCreator
{
	template <uint64_t module>
	struct Module
	{
		typedef FiniteField<BaseForFiniteField> Field;
		static Field Create()
		{
			return Field::CreateZpFieldWithChar(module);
		}
	};
};


template <class Creator>
struct FieldTest :  ::testing::Test
{
	FieldTest()
		:f_(Creator::Create())
	{}
	typename Creator::Field f_;
	typename Creator::Field::Value to_check;
	template <class Integer>
	typename Creator::Field::Value Imp(const Integer& i)
	{
		typename Creator::Field::Value result;
		f_.Import(i, result);
		return result;
	}
};


TYPED_TEST_CASE_P(FieldTest);

TYPED_TEST_P(FieldTest, HasNestedTypes)
{
	EXPECT_GT(sizeof(typename decltype(this->f_)::Value), 0);
	EXPECT_GT(sizeof(typename decltype(this->f_)::DivResult), 0);
}

TYPED_TEST_P(FieldTest, ZeroTests)
{
	//0 == 0
	EXPECT_TRUE(FieldHelpers::IsZero(this->f_, FieldHelpers::Zero(this->f_)));
	EXPECT_TRUE(FieldHelpers::IsEqual(this->f_, FieldHelpers::Zero(this->f_), FieldHelpers::Zero(this->f_)));
	//0-0*0 == 0
	this->f_.Subtract(FieldHelpers::Zero(this->f_), FieldHelpers::Zero(this->f_), FieldHelpers::DivByOne(this->f_, FieldHelpers::Zero(this->f_)), this->to_check);
	EXPECT_TRUE(FieldHelpers::IsZero(this->f_, this->to_check));
	//0-1*0 == 0
	this->f_.Subtract(FieldHelpers::Zero(this->f_), FieldHelpers::One(this->f_), FieldHelpers::DivByOne(this->f_, FieldHelpers::Zero(this->f_)), this->to_check);
	EXPECT_TRUE(FieldHelpers::IsZero(this->f_, this->to_check));
	//0-0*1 == 0
	this->f_.Subtract(FieldHelpers::Zero(this->f_), FieldHelpers::Zero(this->f_), FieldHelpers::DivByOne(this->f_, FieldHelpers::One(this->f_)), this->to_check);
	EXPECT_TRUE(FieldHelpers::IsZero(this->f_, this->to_check));
	//0-0*1 == 0
	this->f_.Subtract(FieldHelpers::Zero(this->f_), FieldHelpers::Zero(this->f_), FieldHelpers::DivByOne(this->f_, FieldHelpers::One(this->f_)), this->to_check);
	EXPECT_TRUE(FieldHelpers::IsZero(this->f_, this->to_check));
	
	//0-0*2 == 0
	this->f_.Subtract(FieldHelpers::Zero(this->f_), FieldHelpers::Zero(this->f_), FieldHelpers::DivByOne(this->f_, this->Imp(2u)), this->to_check);
	EXPECT_TRUE(FieldHelpers::IsZero(this->f_, this->to_check));
	
	//0-0*17 == 0
	//this->f_.Subtract(FieldHelpers::Zero(this->f_), FieldHelpers::Zero(this->f_), FieldHelpers::DivByOne(this->f_, this->Imp(17ull)), this->to_check);
	EXPECT_TRUE(FieldHelpers::IsZero(this->f_, this->to_check));

	//0-1*1 != 0
	this->f_.Subtract(FieldHelpers::Zero(this->f_), FieldHelpers::One(this->f_), FieldHelpers::DivByOne(this->f_, FieldHelpers::One(this->f_)), this->to_check);
	EXPECT_EQ(this->f_.Subtract(FieldHelpers::Zero(this->f_), FieldHelpers::One(this->f_), FieldHelpers::DivByOne(this->f_, FieldHelpers::One(this->f_)), this->to_check), ExactSubtractionResultInfo::NonZero);
	EXPECT_FALSE(FieldHelpers::IsZero(this->f_, this->to_check));
}

TYPED_TEST_P(FieldTest, OneTests)
{
	//0 != 1
	EXPECT_FALSE(FieldHelpers::IsZero(this->f_, FieldHelpers::One(this->f_)));
	EXPECT_FALSE(FieldHelpers::IsEqual(this->f_, FieldHelpers::One(this->f_), FieldHelpers::Zero(this->f_)));
	EXPECT_FALSE(FieldHelpers::IsEqual(this->f_, FieldHelpers::Zero(this->f_), FieldHelpers::One(this->f_)));
	//1 == 1
	EXPECT_TRUE(FieldHelpers::IsEqual(this->f_, FieldHelpers::One(this->f_), FieldHelpers::One(this->f_)));

}

REGISTER_TYPED_TEST_CASE_P(FieldTest, HasNestedTypes, ZeroTests, OneTests);


INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_2, FieldTest, FiniteFieldCreator<ZPlusRing8>::Module<2>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_3, FieldTest, FiniteFieldCreator<ZPlusRing8>::Module<3>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_5, FieldTest, FiniteFieldCreator<ZPlusRing8>::Module<5>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_7, FieldTest, FiniteFieldCreator<ZPlusRing8>::Module<7>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_251, FieldTest, FiniteFieldCreator<ZPlusRing8>::Module<251>);

INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_2, FieldTest, FiniteFieldCreator<ZPlusRing16>::Module<2>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_3, FieldTest, FiniteFieldCreator<ZPlusRing16>::Module<3>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_5, FieldTest, FiniteFieldCreator<ZPlusRing16>::Module<5>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_7, FieldTest, FiniteFieldCreator<ZPlusRing16>::Module<7>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_251, FieldTest, FiniteFieldCreator<ZPlusRing16>::Module<251>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_65521, FieldTest, FiniteFieldCreator<ZPlusRing16>::Module<65521>);

INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_2, FieldTest, FiniteFieldCreator<ZPlusRing32>::Module<2>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_3, FieldTest, FiniteFieldCreator<ZPlusRing32>::Module<3>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_5, FieldTest, FiniteFieldCreator<ZPlusRing32>::Module<5>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_7, FieldTest, FiniteFieldCreator<ZPlusRing32>::Module<7>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_251, FieldTest, FiniteFieldCreator<ZPlusRing32>::Module<251>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_65521, FieldTest, FiniteFieldCreator<ZPlusRing32>::Module<65521>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_4294967291, FieldTest, FiniteFieldCreator<ZPlusRing32>::Module<4294967291>);

TEST(BadConstructionDeathTest, SmallNotPrime)
{
	EXPECT_ASSERT(FiniteFieldCreator<ZPlusRing8>::Module<4>::Create());
}

TEST(BadConstructionDeathTest, BigNotPrime)
{
	EXPECT_ASSERT(FiniteFieldCreator<ZPlusRing32>::Module<4294967295>::Create());
}

TEST(BadConstructionDeathTest, SmallRingBigNum)
{
	EXPECT_ASSERT(FiniteFieldCreator<ZPlusRing8>::Module<257>::Create());
}

TEST(BadConstructionDeathTest, BigRingVeryBigNum)
{
	EXPECT_ASSERT(FiniteFieldCreator<ZPlusRing8>::Module<4294967311ull>::Create());
}
