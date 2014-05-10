#include <limits>
#include <vector>
#include <cstdint>
#include "z_ring.h"
#include "finite_field.h"
#include "combined_field.h"
#include "sparse_matrix.h"

#define PRIVATE_TYPED_TEST SparseMatrixExactTest
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
struct FiniteFieldParam
{
	template <uint64_t module>
	struct Module
	{
		typedef FiniteField<BaseForFiniteField> Field;
		typedef SparseMatrix::Element<Field> Element;
		typedef std::vector<Element> Result;

		struct Matrix
		{
			Matrix():
				field(Field::CreateZpFieldWithChar(module))
			{}
			Field field;
			std::vector<std::vector<Element>> matrix;
			void AddRow()
			{
				matrix.emplace_back();
			}
			template <class Integer>
			void AddElement(int column, Integer value)
			{
				Element element;
				element.column = column;
				field.Import(value, element.value);
				matrix.back().push_back(element);
			}
			Result RunSolver()const
			{
				Result result;
				struct FactoryData{};
				FactoryData factory_data;
				struct Factory
				{
					ExactFieldAsCombined<Field> CreateFieldExpectedAsSuitable(FactoryData& data, const Field& field)
					{
						IgnoreIfUnused(data);
						return ExactFieldAsCombined<Field>(field);
					}
				};
				SparseMatrix::SolveWithRightSideContainigSingleOne<Factory>(field, matrix, result, factory_data);
				return result;
			}
		};
	};
};

#define EXPECT_FIELD_VALUE_EQ(field_value, int_value) EXPECT_TRUE(FieldHelpers::IsEqual(m.field, (field_value), FieldHelpers::Imp(m.field,  int_value)))

template <class Param>
struct SparseMatrixExactTest :  ::testing::Test
{
	typename Param::Matrix m;
	typedef typename Param::Field::Value Value;
	typedef typename Param::Result Result;
	
	void SolvesOneZeroColumn()
	{
		m.AddRow();
		m.AddElement(0, 1u);
		auto result = m.RunSolver();
		EXPECT_EQ(result.size(), 1);
		EXPECT_EQ(result.front().column, 0);
		EXPECT_FIELD_VALUE_EQ(result.front().value, 1u);
	}
	void SolvesOneBigColumn()
	{
	}
	void NotSolvesZeroZeroColumn()
	{
	}
	void NotSolvesZeroBigColumn()
	{
	}
};


TEST_METHOD(SolvesOneZeroColumn)
TEST_METHOD(SolvesOneBigColumn)
TEST_METHOD(NotSolvesZeroZeroColumn)
TEST_METHOD(NotSolvesZeroBigColumn)
REGISTER_TYPED_TEST_CASE_P(SparseMatrixExactTest, SolvesOneZeroColumn, SolvesOneBigColumn, NotSolvesZeroZeroColumn, NotSolvesZeroBigColumn);


INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_2, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing8>::Module<2>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_3, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing8>::Module<3>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_5, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing8>::Module<5>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_7, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing8>::Module<7>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_251, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing8>::Module<251>);

INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_2, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing16>::Module<2>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_3, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing16>::Module<3>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_5, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing16>::Module<5>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_7, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing16>::Module<7>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_251, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing16>::Module<251>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_65521, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing16>::Module<65521>);

INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_2, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing32>::Module<2>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_3, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing32>::Module<3>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_5, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing32>::Module<5>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_7, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing32>::Module<7>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_251, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing32>::Module<251>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_65521, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing32>::Module<65521>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_4294967291, SparseMatrixExactTest, FiniteFieldParam<ZPlusRing32>::Module<4294967291>);

