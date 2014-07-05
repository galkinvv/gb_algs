#include <vector>
#include <cstdint>
#include "sparse_matrix.h"
#include "test_base.h"
	
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
				assert(!matrix.empty()); // has rows
				auto& last_row = matrix.back();
				assert(last_row.empty() || last_row.back().column < column);
				Element element;
				element.column = column;
				field.Import(value, element.value);
				if (!FieldHelpers::IsZero(field, element.value))
				{
					last_row.push_back(element);
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

DECLARE_FUNCTOR_TEMPLATE_T(bool, ExpectGoodSolution, const T& matrix)
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

static struct{//used to expect exact solution
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


static struct{//used to expect partially known solution
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
} ExpectColumnsInSolution;


template <class MatWithField>
void ExpectNoSolution(const MatWithField& matrix)
{
	EXPECT_EQ(matrix.RunSolver().size(), 0);
}
