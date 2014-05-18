#pragma once
#include "finite_field.h"
#include "combined_field.h"
#include "z_ring.h"
#include "utils.h"
#include <algorithm>
#include <vector>
#include <exception>

namespace SparseMatrix
{
template <class Field>
struct Element
{
	static Element FromCV(int column, typename Field::Value value)
	{
		Element result;
		result.column = column;
		result.value = value;
		return result;
	}
	int column;
	typename Field::Value value;
	
		friend std::ostream& operator <<(std::ostream &output, const Element &x) {
			return output << "[" << x.column << "]=" << x.value;
        }

};

namespace{
	
	template <class Field>
	struct RowWithRightPart
	{
		std::vector<Element<Field>> left;
		typename Field::Value right;
	};
	
	template <class Field>
	void SubMatrixRows(const Field& field, const RowWithRightPart<Field>& orig, const RowWithRightPart<Field>& modifier, RowWithRightPart<Field>& result, typename std::vector<Element<Field>>::const_iterator lead_orig, typename std::vector<Element<Field>>::const_iterator lead_modifier)
	{
		typename Field::DivResult row_multiplier;
		assert(lead_orig->column == lead_modifier->column);
		field.Divide(lead_orig->value, lead_modifier->value, row_multiplier);
		field.Subtract(orig.right, modifier.right, row_multiplier, result.right);
		const int result_reserve = orig.left.size() + modifier.left.size() - 2;
		result.left.reserve(result_reserve);
		auto io = orig.left.cbegin();
		auto io_end = orig.left.cend();
		auto im = modifier.left.cbegin();
		auto im_end = modifier.left.cend();
		bool lead_column_found_and_ignored = false;
		while(io!=io_end || im!=im_end)
		{			
			if (im==im_end || io->column < im->column)
			{
				result.left.emplace_back(*io);
				++io;
			}
			else if (io == io_end || io->column > im->column)
			{
				result.left.emplace_back(*im);
				auto sub_result = field.Subtract(FieldHelpers::Zero(field), im->value, row_multiplier, result.left.back().value);
				assert(sub_result == ExactSubtractionResultInfo::NonZero);
				++im;
			}
			else
			{
				assert(io!=io_end && im!=im_end && io->column == im->column);
				if (io != lead_orig)
				{
					typename Field::Value new_value;
					if (ExactSubtractionResultInfo::NonZero == field.Subtract(io->value, im->value, row_multiplier, new_value))
					{
						result.left.emplace_back(*io);
						result.left.back().value = new_value;
					}					
				}
				else
				{
					assert(im == lead_modifier);
					lead_column_found_and_ignored = true;
				}
				++io;
				++im;
			}
		}
		assert(lead_column_found_and_ignored);
		assert(result.left.size() <= result_reserve);
	}
	
	template <class Field> bool columnOfElementLessCompare(const Element<Field>& el, int column)
	{
		return el.column < column;
	}

	//triangulates matrix. lead_columns - filled with column numbers treated as leading.
	//triangulation is in a form that leading element is subtracted from all rows with smaller number
	//i-th number in lead_columns corresponds to i-th row in matrix
	//lead elements of triangulated matrix are NOT set to exact ones - they are arbitrary values

	//throws:
	//incompatible_system_exception when system can't be triangulated (too much non-zero rows for example)
	//cant_detect_zero_equality_exception when combined field detected that n exact part is exactly zero, while approx part isn't
	//unexact_divisor_exception when there is not enough precision in approx part of elements to triangulate
	
	const int kRowIsCompletlyZero = -1; //kRowIsCompletlyZero in lead column means  completely zero row
	struct incompatible_system_exception: std::exception{};
	
	template <class Field>
	void TriangulateMatrix(const Field& field, std::vector<auto_unique_ptr<RowWithRightPart<Field>>>& matrix, std::vector<int>& lead_columns)
	{
		const int kUnitializedLeadColumn = -2; //used only internally
		std::vector<typename std::remove_reference<decltype(matrix.back()->left.cbegin())>::type> min_row_iterators;
		lead_columns.resize(matrix.size(), kUnitializedLeadColumn); //fill with non-initialized
		int used_rows = matrix.size();
		for(;;)
		{			
			min_row_iterators.clear();
			min_row_iterators.reserve(used_rows);
			auto ir=matrix.begin();
			//fill min_row_iterators and check for zero rows
			while (ir != matrix.begin()+used_rows)
			{
				const auto& left_row_part = (**ir).left;
				if (left_row_part.empty())
				{
					// left side is zero
					if (!FieldHelpers::IsZero(field, (**ir).right))
					{
						//zero left side but non-zero right side
						throw incompatible_system_exception();
					}
					//completely zero left and right sides
					--used_rows;
					auto last = matrix.begin()+used_rows;
					if (ir != last)
					{
						ir->swap(*last);
					}
					assert(matrix[used_rows]->left.empty()); //should be empty because we positioned it in a such way
					lead_columns[used_rows] = kRowIsCompletlyZero;
				}
				else
				{					
					//non-empty row
					const auto min_in_row = std::min_element(left_row_part.begin(), left_row_part.end(), [&field](const Element<Field>& v0, const Element<Field>& v1){return field.IsPreciserDivisor(v0.value, v1.value);});					
					assert(min_in_row != left_row_part.end());
					min_row_iterators.push_back(min_in_row);
					++ir;
				}
			}
			if (min_row_iterators.empty())
			{
				//no more rows
				break;
			}
			const auto min_min_row_it = std::min_element(
					min_row_iterators.begin(), 
					min_row_iterators.end(), 
					[&field](decltype(min_row_iterators.front()) v0, decltype(min_row_iterators.front()) v1){return field.IsPreciserDivisor(v0->value, v1->value);}
				);
			
			auto lead_item_in_min_row = *min_min_row_it;
			const int min_row_index = min_min_row_it - min_row_iterators.begin();
			
			const auto found_min_row = matrix.begin()+min_row_index;
			
			--used_rows;
			const auto last = matrix.begin()+used_rows;
			if (found_min_row != last)
			{
				found_min_row->swap(*last);
			}
			assert((matrix[used_rows]->left.cbegin() <= lead_item_in_min_row) && (lead_item_in_min_row < matrix[used_rows]->left.cend())); //should contain the ninimal iterator
			const int lead_column = lead_item_in_min_row->column;
			lead_columns[used_rows] = lead_column;
			const auto& last_row = *last;
			assert((lead_item_in_min_row - last_row->left.begin()) < last_row->left.size());

			//go to subtraction without  normalizing the lead coef in last
			for (ir=matrix.begin(); ir != last; ++ir)
			{
				auto lead_row = std::lower_bound((**ir).left.begin(),(**ir).left.end(), lead_column, columnOfElementLessCompare<Field>);
				if (lead_row == (**ir).left.end() || lead_row->column != lead_column)
				{
					continue;
				}
				auto_unique_ptr<RowWithRightPart<Field>> modified_row;
				SubMatrixRows(field, **ir, *last_row, *modified_row, lead_row, lead_item_in_min_row);
				ir->reset(modified_row.release());
			}			
		}
		//check  anwser corretness: all lead_columns are kRowIsCompletlyZero or correspond to value
		assert(matrix.size() == lead_columns.size());
		for (int ir=0;ir < matrix.size(); ++ir)
		{
			const auto& r = *matrix[ir]; 
			const auto lead_column = lead_columns[ir];
			assert(lead_column != kUnitializedLeadColumn);
			if (lead_column == kRowIsCompletlyZero)
			{
				assert(r.left.empty());
				assert(FieldHelpers::IsZero(field, r.right));
			}
			else
			{
				assert(lead_column >= 0);
				auto lead_it = std::lower_bound(r.left.begin(), r.left.end(), lead_column, columnOfElementLessCompare<Field>);
				assert(lead_it != r.left.end());
				assert(lead_it->column == lead_column);
			}
		}
	}

	//fills combined_matrix with data from matrix, assuming that right side consists if one followed by zeroes
	//expcted  not to throw non-system errors
	template <class CombinedField, class CombinedMatrix, class ElementMatrixContainer>
	void ImportToCombinedMatrix(const CombinedField& combined_field,  CombinedMatrix& combined_matrix, const ElementMatrixContainer& matrix)
	{
		combined_matrix.reserve(matrix.size());
		RandomGenerator rand_functor;
		for (const auto& row: matrix)
		{
			auto& pair_row = emplaced_back(combined_matrix);
			if (combined_matrix.size() == 1)
			{
				//first row is one
				pair_row->right = FieldHelpers::One(combined_field);
			}
			else
			{
				//all other rows - zeroes
				pair_row->right = FieldHelpers::Zero(combined_field);
			}
			pair_row->left.reserve(row.size());
			for(const auto& el:row)
			{
				auto& pair_el = emplaced_back(pair_row->left);
				pair_el.column = el.column;
				combined_field.ExtendWithRandom(el.value, rand_functor, pair_el.value);
				assert(!FieldHelpers::IsZero(combined_field, pair_el.value));
			}
		}
	}
	
	//fills result based on triangulation
	//throws:
	//unexact_divisor_exception when there is not enough precision in approx part of elements to triangulate. This is possible because some lead elements can be not used in inverted form in triangulation

	template <class CombinedField, class CombinedMatrix,  class Field>
	void CreateResultBasedOnTriangulation(const CombinedField& combined_field,  const CombinedMatrix& combined_matrix, const std::vector<int>& lead_columns, const Field& field, std::vector<Element<Field>>& result)
	{
		for (int ir = 0; ir  < lead_columns.size(); ++ir)
		{
			const int lead_column = lead_columns[ir];
			if (lead_column != kRowIsCompletlyZero)
			{
				const auto& r = *combined_matrix[ir];
				const auto lead_element_it = std::lower_bound(r.left.begin(), r.left.end(), lead_column, columnOfElementLessCompare<CombinedField>);
				assert(lead_element_it != r.left.end());
				typename Field::Value lead_value;
				typename Field::Value right_value;
				combined_field.ExtractApprox(lead_element_it->value, lead_value);
				combined_field.ExtractApprox(r.right, right_value);
				//subtract  from ritght side columns corresponding to already foumnd result elements
				for(auto computed_element:result)
				{
					const int other_column = computed_element.column;
					const auto other_element_it = std::lower_bound(r.left.begin(), r.left.end(), other_column, columnOfElementLessCompare<CombinedField>);
					if (other_element_it != r.left.end() && other_element_it->column == other_column)
					{
						typename Field::Value other_value;
						combined_field.ExtractApprox(other_element_it->value, other_value);
						auto other_value_as_frac = FieldHelpers::DivByOne(field, other_value);
						field.Subtract(right_value, computed_element.value, other_value_as_frac, right_value);
					}
				}
				typename Field::DivResult var_value_as_frac;
				field.Divide(right_value, lead_value, var_value_as_frac);
				Element<Field> new_element;
				new_element.column = lead_column;
				new_element.value = FieldHelpers::FracAsValue(field, var_value_as_frac);
				result.emplace_back(new_element);
			}
		}
	}
	
} //namespace{}

//solves matrix equation in field with matrix as left side(Container of rows, where each row is container of elements) assuming that right side consists if one followed by zeroes
//empty result means "no solution"
//expcted  not to throw non-system errors
template <class CombinedFieldFactory, class Field, class ElementMatrixContainer,  class CombinedFieldFactoryData>
void SolveWithRightSideContainigSingleOne(const Field& field, const ElementMatrixContainer& matrix, std::vector<Element<Field>>& result, CombinedFieldFactoryData& combined_creator_data)
{
	CombinedFieldFactory combined_field_factory;
	static const int max_attempts = 2;//maximal number of attmpts to select suitable field
	int attempt = 0;
	result.clear();
	try
	{
		for (;;)
		{
			auto combined_field = combined_field_factory.CreateFieldExpectedAsSuitable(combined_creator_data, field);
			std::vector<auto_unique_ptr<RowWithRightPart<decltype(combined_field)>>> combined_matrix;
			ImportToCombinedMatrix(combined_field, combined_matrix, matrix);
			std::vector<int> lead_columns;
			try
			{
				TriangulateMatrix(combined_field, combined_matrix, lead_columns);
			}
			catch(const cant_detect_zero_equality_exception&)
			{
				++attempt;
				if (attempt >= max_attempts)
				{
					assert(result.empty());
					return;//calculations completed without solution found, because no suitable exact field found
				}
				continue;//non-suitable field chosen, try another
			}
			catch(const incompatible_system_exception&)
			{
				assert(result.empty());
				return;//calculations completed without solution found, because no suitable exact field found
				//system is incompatible 
			}
			assert(lead_columns.size() == combined_matrix.size());
			CreateResultBasedOnTriangulation(combined_field, combined_matrix, lead_columns, field, result);
			assert(!result.empty());
			return;//calculations completed with suitable field, all done
		}
	}
	catch(const unexact_divisor_exception&)
	{
		//not enough data precision to solve system
		result.clear();
		return;
	}
}
}
