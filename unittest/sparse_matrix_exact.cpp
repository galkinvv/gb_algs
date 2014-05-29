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
				if (!FieldHelpers::IsZero(field, element.value))
				{
					matrix.back().push_back(element);
				}
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
			void Clear()
			{
				matrix.clear();
			}
			
			friend std::ostream& operator <<(std::ostream &output, const Matrix &x) {
				return output << "Matrix with filed " << x.field;
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
		m.AddRow();
		m.AddElement(42, 1u);
		auto result = m.RunSolver();
		EXPECT_EQ(result.size(), 1);
		EXPECT_EQ(result.front().column, 42);
		EXPECT_FIELD_VALUE_EQ(result.front().value, 1u);
	}
	void SolvesOneAndZeroInColumn()
	{
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddRow();
		auto result = m.RunSolver();
		EXPECT_EQ(result.size(), 1);
		EXPECT_EQ(result.front().column, 0);
		EXPECT_FIELD_VALUE_EQ(result.front().value, 1u);
	}
	void SolvesOneAndManyZeroInColumn()
	{
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddRow();
		m.AddRow();
		m.AddRow();
		m.AddRow();
		auto result = m.RunSolver();
		EXPECT_EQ(result.size(), 1);
		EXPECT_EQ(result.front().column, 0);
		EXPECT_FIELD_VALUE_EQ(result.front().value, 1u);
	}
	void SolvesTwoColumn()
	{
		m.AddRow();
		static const unsigned col_values[] = {2,3};
		m.AddElement(0, col_values[0]);
		m.AddElement(1, col_values[1]);
		auto result = m.RunSolver();
		EXPECT_GT(result.size(), 0);
		EXPECT_LE(result.size(), 2);
	}
	void SolvesTwoColumnTwoRows()
	{
		m.AddRow();
		static const unsigned col_values[] = {2,3};
		m.AddElement(0, col_values[0]);
		m.AddElement(1, col_values[1]);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		auto result = m.RunSolver();
		EXPECT_GT(result.size(), 0);
		EXPECT_LE(result.size(), 2);
	}
	void SolvesTwoColumnManyRows()
	{
		m.AddRow();
		static const unsigned col_values[] = {2,3};
		m.AddElement(0, col_values[0]);
		m.AddElement(1, col_values[1]);
		m.AddRow();
		m.AddElement(0, 2u);
		m.AddElement(1, 2u);
		m.AddRow();
		m.AddElement(0, 3u);
		m.AddElement(1, 3u);
		m.AddRow();
		m.AddElement(0, 255u);
		m.AddElement(1, 255u);
		auto result = m.RunSolver();
		EXPECT_GT(result.size(), 0);
		EXPECT_LE(result.size(), 2);
	}
	void NotSolvesZeroAndOneInColumn()
	{
		m.AddRow();
		m.AddRow();
		m.AddElement(0, 1u);
		auto result = m.RunSolver();
		EXPECT_EQ(result.size(), 0);
	}
	void NotSolvesTwoOnesTwoColumns()
	{
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		auto result = m.RunSolver();
		EXPECT_EQ(result.size(), 0);
	}
	void NotSolvesTwoOnes()
	{
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		auto result = m.RunSolver();
		EXPECT_EQ(result.size(), 0);
	}
	void NotSolvesZero()
	{
		m.AddRow();
		auto result = m.RunSolver();
		EXPECT_EQ(result.size(), 0);
	}
	void NotSolvesTwoZeroes()
	{
		m.AddRow();
		m.AddRow();
		auto result = m.RunSolver();
		EXPECT_EQ(result.size(), 0);
	}
};


TEST_METHOD(SolvesOneZeroColumn)
TEST_METHOD(SolvesOneBigColumn)
TEST_METHOD(SolvesOneAndZeroInColumn)
TEST_METHOD(SolvesOneAndManyZeroInColumn)
TEST_METHOD(SolvesTwoColumn)
TEST_METHOD(SolvesTwoColumnTwoRows)
TEST_METHOD(SolvesTwoColumnManyRows)
TEST_METHOD(NotSolvesZeroAndOneInColumn)
TEST_METHOD(NotSolvesTwoOnesTwoColumns)
TEST_METHOD(NotSolvesTwoOnes)
TEST_METHOD(NotSolvesZero)
TEST_METHOD(NotSolvesTwoZeroes)
REGISTER_TYPED_TEST_CASE_P(SparseMatrixExactTest, SolvesOneZeroColumn, SolvesOneBigColumn, SolvesOneAndZeroInColumn, SolvesOneAndManyZeroInColumn, SolvesTwoColumn, SolvesTwoColumnTwoRows, SolvesTwoColumnManyRows, 
NotSolvesZeroAndOneInColumn, NotSolvesTwoOnesTwoColumns, NotSolvesTwoOnes, NotSolvesZero, NotSolvesTwoZeroes);


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

static struct
{
	template <class MatWithField>
	bool operator()(const MatWithField& matrix)const
	{
		EXPECT_FUNCTION_BEGIN
		auto result = matrix.RunSolver();
		EXPECT_2(EqualTo, result.empty(), false);
		for(auto& row:matrix.matrix)
		{
			auto negated_sum = FieldHelpers::Zero(matrix.field);
			auto last_subtraction_result = ExactSubtractionResultInfo::Zero;
			for (auto& res_elem:result)
			{
				for (auto& row_elem:row)
				{
					if (res_elem.column == row_elem.column)
					{
						last_subtraction_result = matrix.field.Subtract(negated_sum, res_elem.value, FieldHelpers::DivByOne(matrix.field, row_elem.value), negated_sum);					
					}
				}
			}
			bool is_first = (&row == &matrix.matrix.front());
			if (is_first)
			{
				EXPECT_2(EqualTo, last_subtraction_result, ExactSubtractionResultInfo::NonZero);
				decltype(negated_sum) minus_one;
				matrix.field.Subtract(FieldHelpers::Zero(matrix.field), FieldHelpers::One(matrix.field), FieldHelpers::DivByOne(matrix.field, FieldHelpers::One(matrix.field)), minus_one);
				EXPECT_3(FieldHelpers::IsEqual<decltype(matrix.field)>, matrix.field, negated_sum, minus_one);
			}
			else
			{
				EXPECT_2(EqualTo, last_subtraction_result, ExactSubtractionResultInfo::Zero);
				EXPECT_2(FieldHelpers::IsZero<decltype(matrix.field)>, matrix.field, negated_sum);
			}
		}
		EXPECT_FUNCTION_RETURN
	}
} ExpectGoodSolution;

static struct
{
	template <class MatWithField>
	bool operator()(const MatWithField& matrix, const std::initializer_list<typename decltype(matrix.RunSolver())::value_type>& list)const
	{
		EXPECT_FUNCTION_BEGIN
		auto result = matrix.RunSolver();
		typedef  typename decltype(matrix.RunSolver())::value_type Element;
		//due to treating field as inexact zeroes can be in result. Erase them.
		result.erase(std::remove_if(result.begin(), result.end(), [&](const Element& e){return FieldHelpers::IsZero(matrix.field, e.value);}), result.end());
		EXPECT_2(EqualTo, result.size(), list.size());
		struct LessColumn{
			bool operator()(const Element& e1, const Element& e2)
			{
				return e1.column < e2.column;
			}
		};
		std::sort(result.begin(), result.end(), LessColumn());
		const auto eq_element = [&](const Element& e1, const Element& e2)
		{
			return e1.column == e2.column && FieldHelpers::IsEqual(matrix.field, e1.value, e2.value);
		};
		EXPECT_2(ExpecterContainerEqual(eq_element), result, list);
		EXPECT_1(ExpectGoodSolution, matrix);
		EXPECT_FUNCTION_RETURN
	}
} ExpectKnownSolution;


template <class MatWithField>
void ExpectNoSolution(const MatWithField& matrix)
{
	EXPECT_EQ(matrix.RunSolver().size(), 0);
}

TEST(SparseMatrixExactValues, Z2determined)
{
	typedef FiniteFieldParam<ZPlusRing32>::Module<2> Param;
	Param::Matrix m;
	const auto E = [&](int column, unsigned long long value){return Param::Element::FromCV(column, FieldHelpers::Imp(m.field, value));};
	{
		//zero matrix size 1
		m.Clear();
		m.AddRow();
		ExpectNoSolution(m);
	}
	{
		//zero matrix size 3
		m.Clear();
		m.AddRow();
		m.AddRow();
		m.AddRow();
		ExpectNoSolution(m);
	}
	{
		//identity matrix size 1
		m.Clear();
		m.AddRow();
		m.AddElement(0, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(0, 1u)}));
	}
	{
		//identity matrix size 3
		m.Clear();
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddRow();
		m.AddElement(2, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(0, 1u)}));
	}
	{
		//matrix with big column number
		m.Clear();
		m.AddRow();
		m.AddElement(42, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(42, 1u)}));
	}
	{
		//matrix with single zero column
		m.Clear();
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		ExpectNoSolution(m);
	}
	{
		//matrix with single columnwith ones
		m.Clear();
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		ExpectNoSolution(m);
	}
	{
		//matrix with single columnw eqaul to result, non-obviousely initialized
		m.Clear();
		m.AddRow();
		m.AddElement(0, 3u);
		m.AddRow();
		m.AddElement(0, 2u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(0, 1u)}));
	}
	{
		//matrix with every column needed in sum
		m.Clear();
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddElement(3, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(3, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(2, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(0, 1u), E(1, 1u), E(2, 1u), E(3, 1u)}));
	}
}
