#include <limits>
#include "test_base.h"
#include "z_ring.h"
#include "finite_field.h"

template <class Field>
struct CreateFiledForTest
{
	static Field Create()
	{
		return Field();
	}
};

template <class BaseForFiniteField>
struct CreateFiledForTest<FiniteField<BaseForFiniteField>>
{
	static FiniteField<BaseForFiniteField> Create()
	{
		return FiniteField<BaseForFiniteField>(std::numeric_limits<typename BaseForFiniteField::Value>::max());
	}
};

template <class Field>
struct FieldTest :  ::testing::Test
{
	FieldTest()
	:f_(CreateFiledForTest<Field>::Create())
	{}
	Field f_;
};


TYPED_TEST_CASE_P(FieldTest);

TYPED_TEST_P(FieldTest, UseValueType)
{
	ASSERT_GT(sizeof(typename decltype(this->f_)::Value), 0);
}

REGISTER_TYPED_TEST_CASE_P(FieldTest, UseValueType);


INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8, FieldTest, FiniteField<ZRing8>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16, FieldTest, FiniteField<ZRing16>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32, FieldTest, FiniteField<ZRing32>);