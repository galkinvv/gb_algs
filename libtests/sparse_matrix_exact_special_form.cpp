#include <limits>
#include "z_ring.h"
#include "finite_field.h"
#include "combined_field.h"

#define PRIVATE_TYPED_TEST SparseMatrixSpecialFormTest
#include "sparse_matrix_base.h"

template <class Param>
struct SparseMatrixSpecialFormTest :  ::testing::Test
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
REGISTER_TYPED_TEST_CASE_P(SparseMatrixSpecialFormTest, SolvesOneZeroColumn, SolvesOneBigColumn, SolvesOneAndZeroInColumn, SolvesOneAndManyZeroInColumn, SolvesTwoColumn, SolvesTwoColumnTwoRows, SolvesTwoColumnManyRows, 
NotSolvesZeroAndOneInColumn, NotSolvesTwoOnesTwoColumns, NotSolvesTwoOnes, NotSolvesZero, NotSolvesTwoZeroes);


INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_2, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing8>::Module<2>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_3, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing8>::Module<3>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_5, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing8>::Module<5>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_7, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing8>::Module<7>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField8_251, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing8>::Module<251>);

INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_2, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing16>::Module<2>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_3, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing16>::Module<3>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_5, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing16>::Module<5>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_7, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing16>::Module<7>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_251, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing16>::Module<251>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField16_65521, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing16>::Module<65521>);

INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_2, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing32>::Module<2>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_3, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing32>::Module<3>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_5, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing32>::Module<5>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_7, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing32>::Module<7>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_251, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing32>::Module<251>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_65521, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing32>::Module<65521>);
INSTANTIATE_TYPED_TEST_CASE_P(FiniteField_ZField32_4294967291, SparseMatrixSpecialFormTest, FiniteFieldParam<ZPlusRing32>::Module<4294967291>);

