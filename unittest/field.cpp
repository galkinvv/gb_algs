#include <limits>
#include <cstdint>
#include "z_ring.h"
#include "finite_field.h"

#define PRIVATE_TYPED_TEST FieldTest
#include "test_base.h"

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

template <class Param>
struct FieldTest :  ::testing::Test
{
	FieldTest()
		:f_(Param::Create())
	{}
	typename Param::Field f_;
	typedef typename Param::Field::Value Value;
	typedef typename Param::Field::DivResult DivResult;
	
	template <class Integer>
	Value Imp(const Integer& i)
	{
		Value result;
		f_.Import(i, result);
		return result;
	}

	Value MinusOne()
	{
		Value result;
		f_.Subtract(FieldHelpers::Zero(f_), FieldHelpers::One(f_), FieldHelpers::DivByOne(f_, FieldHelpers::One(f_)), result);
		return result;
	}

	void NonEmptyNestedTypes()
	{
		EXPECT_GT(sizeof(Value), 0);
		EXPECT_GT(sizeof(DivResult), 0);		
	}

	void ZeroTests()
	{
		Value to_check;
		//0 == 0
		EXPECT_TRUE(FieldHelpers::IsZero(f_, FieldHelpers::Zero(f_)));
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, FieldHelpers::Zero(f_), Imp(0u)));
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, FieldHelpers::Zero(f_), FieldHelpers::Zero(f_)));
		//0-0*0 == 0
		f_.Subtract(FieldHelpers::Zero(f_), FieldHelpers::Zero(f_), FieldHelpers::DivByOne(f_, FieldHelpers::Zero(f_)), to_check);
		EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));
		//0-1*0 == 0
		f_.Subtract(FieldHelpers::Zero(f_), FieldHelpers::One(f_), FieldHelpers::DivByOne(f_, FieldHelpers::Zero(f_)), to_check);
		EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));
		//0-0*1 == 0
		f_.Subtract(FieldHelpers::Zero(f_), FieldHelpers::Zero(f_), FieldHelpers::DivByOne(f_, FieldHelpers::One(f_)), to_check);
		EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));
		//0-0*1 == 0
		f_.Subtract(FieldHelpers::Zero(f_), FieldHelpers::Zero(f_), FieldHelpers::DivByOne(f_, FieldHelpers::One(f_)), to_check);
		EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));
		
		//0-0*2 == 0
		f_.Subtract(FieldHelpers::Zero(f_), FieldHelpers::Zero(f_), FieldHelpers::DivByOne(f_, Imp(2u)), to_check);
		EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));
		
		//0-0*17 == 0
		f_.Subtract(FieldHelpers::Zero(f_), FieldHelpers::Zero(f_), FieldHelpers::DivByOne(f_, Imp(17ull)), to_check);
		EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));

		//0-1*1 != 0
		f_.Subtract(FieldHelpers::Zero(f_), FieldHelpers::One(f_), FieldHelpers::DivByOne(f_, FieldHelpers::One(f_)), to_check);
		EXPECT_EQ(f_.Subtract(FieldHelpers::Zero(f_), FieldHelpers::One(f_), FieldHelpers::DivByOne(f_, FieldHelpers::One(f_)), to_check), ExactSubtractionResultInfo::NonZero);
		EXPECT_FALSE(FieldHelpers::IsZero(f_, to_check));		
	}

	void OneTests()
	{
		Value to_check;
		
		//0 != 1
		EXPECT_FALSE(FieldHelpers::IsZero(f_, FieldHelpers::One(f_)));
		EXPECT_FALSE(FieldHelpers::IsEqual(f_, FieldHelpers::One(f_), FieldHelpers::Zero(f_)));
		EXPECT_FALSE(FieldHelpers::IsEqual(f_, FieldHelpers::Zero(f_), FieldHelpers::One(f_)));
		//1 == 1
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, FieldHelpers::One(f_), Imp(1u)));
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, FieldHelpers::One(f_), FieldHelpers::One(f_)));	
		//1-1*1=0
		f_.Subtract(FieldHelpers::One(f_), FieldHelpers::One(f_), FieldHelpers::DivByOne(f_, FieldHelpers::One(f_)), to_check);
		EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));
		
		//17-1*1=16
		f_.Subtract(Imp(17u), FieldHelpers::One(f_), FieldHelpers::DivByOne(f_, FieldHelpers::One(f_)), to_check);
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, to_check, Imp(16u)));

		//255-1*1=254
		f_.Subtract(Imp(255u), FieldHelpers::One(f_), FieldHelpers::DivByOne(f_, FieldHelpers::One(f_)), to_check);
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, to_check, Imp(254u)));
		
		//17-1*17=0
		f_.Subtract(Imp(17u), FieldHelpers::One(f_), FieldHelpers::DivByOne(f_, Imp(17u)), to_check);
		EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));	

		//2-1*2=0
		f_.Subtract(Imp(2u), FieldHelpers::One(f_), FieldHelpers::DivByOne(f_, Imp(2u)), to_check);
		EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));	

		//255-1*255=0
		f_.Subtract(Imp(255u), FieldHelpers::One(f_), FieldHelpers::DivByOne(f_, Imp(255u)), to_check);
		EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));

		//2-1*1=1
		f_.Subtract(Imp(2u), FieldHelpers::One(f_), FieldHelpers::DivByOne(f_, FieldHelpers::One(f_)), to_check);
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, to_check, FieldHelpers::One(f_)));
	}

	void MinusOneTests()
	{
		EXPECT_FALSE(FieldHelpers::IsZero(f_, MinusOne()));
		Value to_check;
		//1 == -1 if and only if 0 ==2 (field with char two)
		bool oneEqualToMinusOne = FieldHelpers::IsEqual(f_, FieldHelpers::One(f_), MinusOne());
		bool twoIsZero= FieldHelpers::IsZero(f_, Imp(2u));
		EXPECT_EQ(oneEqualToMinusOne, twoIsZero);
		
		// 0-1*(-1/1) ==1
		f_.Subtract(FieldHelpers::Zero(f_), FieldHelpers::One(f_), FieldHelpers::DivByOne(f_, MinusOne()), to_check);
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, to_check, FieldHelpers::One(f_)));	
		
		DivResult minusOneByMinusOne;
		f_.Divide(MinusOne(), MinusOne(), minusOneByMinusOne);

		// 0-(-1)*(-1/-1) ==1
		f_.Subtract(FieldHelpers::Zero(f_), MinusOne(), minusOneByMinusOne, to_check);
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, to_check, FieldHelpers::One(f_)));	

		// 0-1*(-1/-1) == -1
		f_.Subtract(FieldHelpers::Zero(f_), FieldHelpers::One(f_), minusOneByMinusOne, to_check);
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, to_check, MinusOne()));	

		//17-(-1)*1=18
		f_.Subtract(Imp(17u), MinusOne(), FieldHelpers::DivByOne(f_, FieldHelpers::One(f_)), to_check);
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, to_check, Imp(18u)));

		//254-(-1)*1=254
		f_.Subtract(Imp(254u), MinusOne(), FieldHelpers::DivByOne(f_, FieldHelpers::One(f_)), to_check);
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, to_check, Imp(255u)));
	}
	
	void SubtractTests()
	{
		Value to_check;
		//255-25*10==5
		f_.Subtract(Imp(255u), Imp(25u), FieldHelpers::DivByOne(f_, Imp(10u)), to_check);
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, to_check, Imp(5u)));		

		//15-2*8==-1
		f_.Subtract(Imp(15u), Imp(2u), FieldHelpers::DivByOne(f_, Imp(8u)), to_check);
		EXPECT_TRUE(FieldHelpers::IsEqual(f_, to_check, MinusOne()));		
	}
	void DivideTests()
	{
		unsigned dividers []={2,3,5,7,61};
		bool skipped = false;
		for(auto divider:dividers)
		{
			Value imp_divider = Imp(divider);
			if (FieldHelpers::IsZero(f_, imp_divider))
			{
				EXPECT_FALSE(skipped);
				skipped = true;
			}
			else
			{
				DivResult d1, d2, d25, d255;
				f_.Divide(Imp(1u), imp_divider, d1);
				f_.Divide(Imp(2u), imp_divider, d2);
				f_.Divide(Imp(25u), imp_divider, d25);
				f_.Divide(Imp(255u), imp_divider, d255);
				
				Value to_check, premult;

				f_.Subtract(Imp(1u), imp_divider, d1, to_check);
				EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));
				
				f_.Subtract(Imp(2u), imp_divider, d2, to_check);
				EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));
				
				f_.Subtract(Imp(25u), imp_divider, d25, to_check);
				EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));
				
				f_.Subtract(Imp(255u), imp_divider, d255, to_check);
				EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));
				
				f_.Subtract(FieldHelpers::Zero(f_), Imp(5u), d1, premult);
				f_.Subtract(FieldHelpers::Zero(f_), imp_divider, FieldHelpers::DivByOne(f_, premult), to_check);
				EXPECT_TRUE(FieldHelpers::IsEqual(f_, to_check, Imp(5u)));

				f_.Subtract(FieldHelpers::Zero(f_), Imp(7u), d25, premult);
				f_.Subtract(FieldHelpers::Zero(f_), imp_divider, FieldHelpers::DivByOne(f_, premult), to_check);
				EXPECT_TRUE(FieldHelpers::IsEqual(f_, to_check, Imp(175u)));
				
				if (divider == 2)
				{
					//-1 == 2/2
					f_.Subtract(FieldHelpers::Zero(f_), FieldHelpers::One(f_), d2, to_check);
					EXPECT_TRUE(FieldHelpers::IsEqual(f_, to_check, MinusOne()));
				}
				if (divider == 3)
				{
					//85 == 255/3
					f_.Subtract(Imp(85u), FieldHelpers::One(f_), d255, to_check);
					EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));					
				}
				if (divider == 5)
				{
					//5 == 25/5
					f_.Subtract(Imp(5u), FieldHelpers::One(f_), d25, to_check);
					EXPECT_TRUE(FieldHelpers::IsZero(f_, to_check));					
				}
			}
		}
	}

};


TEST_METHOD(NonEmptyNestedTypes)
TEST_METHOD(ZeroTests)
TEST_METHOD(OneTests)
TEST_METHOD(MinusOneTests)
TEST_METHOD(SubtractTests)
TEST_METHOD(DivideTests)
REGISTER_TYPED_TEST_CASE_P(FieldTest, NonEmptyNestedTypes, ZeroTests, OneTests, MinusOneTests, SubtractTests, DivideTests);


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

TEST(ImportTest, ring8)
{
	typedef FiniteFieldCreator<ZPlusRing8>::Module<3> Creator;
	auto field = Creator::Create();
	typename Creator::Field::Value value0, value1;
	field.Import(uint8_t(2), value0);
	field.Import(uint8_t(5), value1);
	EXPECT_TRUE(FieldHelpers::IsEqual(field, value0, value1));
	field.Import(uint8_t(3), value0);
	field.Import(uint8_t(6), value1);
	EXPECT_TRUE(FieldHelpers::IsEqual(field, value0, value1));
	EXPECT_TRUE(FieldHelpers::IsZero(field, value0));
	EXPECT_TRUE(FieldHelpers::IsZero(field, value1));
	field.Import(uint8_t(4), value0);
	field.Import(uint8_t(1), value1);
	EXPECT_TRUE(FieldHelpers::IsEqual(field, value0, value1));
	field.Import(uint8_t(254), value0);
	field.Import(uint8_t(251), value1);
	EXPECT_TRUE(FieldHelpers::IsEqual(field, value0, value1));
	field.Import(uint8_t(252), value0);
	field.Import(uint8_t(255), value1);
	EXPECT_TRUE(FieldHelpers::IsEqual(field, value0, value1));
	field.Import(uint64_t(255), value0);
	field.Import(uint64_t(252), value1);
	EXPECT_TRUE(FieldHelpers::IsEqual(field, value0, value1));
}

TEST(ImportTest, ring32)
{
	typedef FiniteFieldCreator<ZPlusRing32>::Module<257> Creator;
	auto field = Creator::Create();
	typename Creator::Field::Value value0, value1;
	field.Import(uint16_t(2), value0);
	field.Import(uint16_t(2+257), value1);
	EXPECT_TRUE(FieldHelpers::IsEqual(field, value0, value1));
	field.Import(uint64_t(257<<13), value0);
	field.Import(uint32_t(257), value1);
	EXPECT_TRUE(FieldHelpers::IsEqual(field, value0, value1));
	EXPECT_TRUE(FieldHelpers::IsZero(field, value0));
	EXPECT_TRUE(FieldHelpers::IsZero(field, value1));
	field.Import(uint8_t(4), value0);
	field.Import(uint8_t(1), value1);
	EXPECT_FALSE(FieldHelpers::IsEqual(field, value0, value1));
	field.Import(uint8_t(254), value0);
	field.Import(uint8_t(251), value1);
	EXPECT_FALSE(FieldHelpers::IsEqual(field, value0, value1));
	field.Import(uint8_t(252), value0);
	field.Import(uint8_t(255), value1);
	EXPECT_FALSE(FieldHelpers::IsEqual(field, value0, value1));
	field.Import(uint64_t(0xFFFFFFFFu), value0);
	field.Import(uint64_t(0xFFFFFFFFu-257u), value1);
	EXPECT_TRUE(FieldHelpers::IsEqual(field, value0, value1));
}


TEST(BadImportDeathTest, ring8_16)
{
	typedef FiniteFieldCreator<ZPlusRing8>::Module<3> Creator;
	auto field = Creator::Create();
	typename Creator::Field::Value value;
	EXPECT_ASSERT(field.Import(uint16_t(256), value));
}

TEST(BadImportDeathTest, ring8_32)
{
	typedef FiniteFieldCreator<ZPlusRing8>::Module<3> Creator;
	auto field = Creator::Create();
	typename Creator::Field::Value value;
	EXPECT_ASSERT(field.Import(uint32_t(256), value));
}

TEST(BadImportDeathTest, ring8_64)
{
	typedef FiniteFieldCreator<ZPlusRing8>::Module<3> Creator;
	auto field = Creator::Create();
	typename Creator::Field::Value value;
	EXPECT_ASSERT(field.Import(uint64_t(256), value));
}

TEST(BadImportDeathTest, ring8_64_big)
{
	typedef FiniteFieldCreator<ZPlusRing8>::Module<3> Creator;
	auto field = Creator::Create();
	typename Creator::Field::Value value;
	EXPECT_ASSERT(field.Import(0xFFFFFFFFFFFFFFFFull, value));
}

TEST(BadImportDeathTest, ring32_64)
{
	typedef FiniteFieldCreator<ZPlusRing32>::Module<2> Creator;
	auto field = Creator::Create();
	typename Creator::Field::Value value;
	EXPECT_ASSERT(field.Import(0x100000000ull, value));
}

TEST(BadImportDeathTest, ring32_64_big)
{
	typedef FiniteFieldCreator<ZPlusRing32>::Module<2> Creator;
	auto field = Creator::Create();
	typename Creator::Field::Value value;
	EXPECT_ASSERT(field.Import(0xFFFFFFFFFFFFFFFFull, value));
}


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
