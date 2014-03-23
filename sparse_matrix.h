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
	typename Field::Value value;
	int column;
};

namespace{
	
	template <class Field>
	struct RowWithRightPart
	{
		std::vector<Element<Field>> left;
		Element<Field> right;
	};
	
	//triangulates matrix. lead_columns - filled with column numbers treated as leading.
	//i-th number in lead_columns corresponds to i-th row in matrix
	//lead elements of triangulated matrix are set to exact ones

	const int kRowIsCompletlyZero = -1; //kRowIsCompletlyZero in lead column means  completely zero row
	const int kUnitializedLeadColumn = -2; //used only internally
	struct incompatible_system_exception: std::exception{};
	//throws incompatible_system_exception when system can't be triangulated (too much non-zero rows for example)
	
	template <class Field>
	void TriangulateMatrix(const Field& field, std::vector<auto_unique_ptr<RowWithRightPart<Field>>>& matrix, std::vector<int>& lead_columns)
	{
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
					if (!FieldHelpers::IsZero(field, (**ir).right.value))
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
					assert(matrix[used_rows]->left.empty()); //should be empty because we positioned in in a such way
					lead_columns[used_rows] = kRowIsCompletlyZero;
				}
				else
				{					
					//non-empty row
					const auto min_in_row = std::min_element(left_row_part.begin(), left_row_part.end(), [&field](const Element<Field>& v0, const Element<Field>& v1){return field.IsPreciserDivisor(v0.value, v1.value);});					
					assert(min_in_row != left_row_part.end());
					min_row_iterators.push_back(min_in_row);
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
			
			auto item_in_min_row = *min_min_row_it;
			const int min_row_index = min_min_row_it - min_row_iterators.begin();
			
			const auto found_min_row = matrix.begin()+min_row_index;
			
			--used_rows;
			const auto last = matrix.begin()+used_rows;
			if (found_min_row != last)
			{
				found_min_row->swap(*last);
			}
			assert(matrix[used_rows]->left.empty()); //should be empty because we positioned in in a such way
			lead_columns[used_rows] = item_in_min_row - (**last).left.begin();
			assert(lead_columns[used_rows] < (**last).left.size());

			//TODO: normilize last
			for (ir=matrix.begin(); ir != last; ++ir)
			{
				//TODO: substract multiplied last from ir
			}			
		}
		//TODO: check  anwser corretness: all lead_columns are -1 or correspond to exaclyvalue exactly equal to one 
	}
}

//solves matrix equation in field with matrix as left side(Container of rows, where each row is container of elements) assuming that right side consists if one followed by zeroes
//empty result means "no solution"
template <class Field, class ElementMatrixContainer, class CombinedField = ExactFieldAsCombined<Field>>
void SolveWithRightSideContainigSingleOne(const Field& field, const ElementMatrixContainer& matrix, std::vector<Element<Field>>& result, int max_diferent_numbers_in_coefficients)
{
	result.clear();
	//TODO for non-Z2 case:
	//calculate P_0 based on matrix.size() and max_diferent_numbers_in_coefficients
	//calculate N_0 based on P_0. 
	//select field based on P_0 and N_0
	auto combined_field = CombinedField(field);
	std::vector<auto_unique_ptr<RowWithRightPart<CombinedField>>> combined_matrix;
	combined_matrix.reserve(matrix.size());
	RandomGenerator rand_functor;
	for (const auto& row: matrix)
	{
		auto& pair_row = emplaced_back(combined_matrix);
		if (combined_matrix.size() == 1)
		{
			//first row contains one
			combined_field.SetOne(pair_row->right.value);
		}
		else
		{
			//all other rows - zeroes
			combined_field.SetZero(pair_row->right.value);
		}
		pair_row->left.reserve(row.size());
		for(const auto& el:row)
		{
			auto& pair_el = emplaced_back(pair_row->left);
			pair_el.column = el.column;
			combined_field.ExtendWithRandom(el.value, rand_functor, pair_el.value);
		}
	}
	std::vector<int> lead_columns;
	TriangulateMatrix(combined_field, combined_matrix, lead_columns);
	//TODO
}
}
