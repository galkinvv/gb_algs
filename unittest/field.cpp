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
};


TYPED_TEST_CASE_P(FieldTest);

TYPED_TEST_P(FieldTest, UseValueType)
{
	ASSERT_GT(sizeof(typename decltype(this->f_)::Value), 0);
}

TYPED_TEST_P(FieldTest, ZeroIsZero)
{
	ASSERT_TRUE(FieldHelpers::IsZero(this->f_, FieldHelpers::Zero(this->f_)));
}

TYPED_TEST_P(FieldTest, OneIsNotZero)
{
	ASSERT_FALSE(FieldHelpers::IsZero(this->f_, FieldHelpers::One(this->f_)));
}

REGISTER_TYPED_TEST_CASE_P(FieldTest, UseValueType, ZeroIsZero, OneIsNotZero);


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
