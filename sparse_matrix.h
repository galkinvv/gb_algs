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
	const int kRowIsCompletlyZero = -1;
	//kRowIsCompletlyZero in lead column means  completely zero row
	struct incompatible_system_exception: std::exception{};
	//throws incompatible_system_exception when system can't be triangulated (too much non-zero rows for example)
	
	template <class Field>
	void TriangulateMatrix(const Field& field, std::vector<auto_unique_ptr<RowWithRightPart<Field>>>& matrix, std::vector<int>& lead_columns)
	{
		std::vector<typename std::remove_reference<decltype(matrix.back()->left.cbegin())>::type> min_row_iterators;
		lead_columns.resize(matrix.size(), -2); //fill with non-initialized
		for (int used_rows = matrix.size(); used_rows >0 ; --used_rows)
		{			
			const auto used_matrix_end = matrix.begin()+used_rows;
			int last_idx = used_rows - 1;
			//check for zero rows
			bool zero_row_found = false;
			for (auto ir=matrix.begin(); ir!=used_matrix_end; ++ir)
			{
				const auto& r = **ir;
				if (r.left.empty())
				{
					if (!FieldHelpers::IsZero(field, r.right.value))
					{
						//zero left side but non-zero right side
						throw incompatible_system_exception();
					}
					//completely zero left and right sides
					auto last = matrix.begin() + last_idx;
					if (ir != last)
					{
						ir->swap(*last);
					}
					assert(matrix[last_idx]->left.empty()); //should be empty because we positioned in in a such way
					lead_columns[last_idx] = kRowIsCompletlyZero;
					zero_row_found = true;
					break;
				}
			}
			if (zero_row_found)
			{
				continue;
			}
			//fill min_row_iterators
			min_row_iterators.clear();
			min_row_iterators.reserve(used_rows);
			for (auto ir = matrix.begin(); ir != used_matrix_end; ++ir)
			{
				const auto& left_row_part = (**ir).left;
				assert(!left_row_part.empty());
				const auto min_in_row = std::min_element(left_row_part.begin(), left_row_part.end(), [&field](const Element<Field>& v0, const Element<Field>& v1){return field.IsPreciserDivisor(v0.value, v1.value);});
				min_row_iterators.push_back(min_in_row);
			}
			//TODO
			//find minimail
			//swap
			//subtract
			//fill lead_columns
		}
	}
}

//solves matrix equation in field with matrix as left side(Container of rows, where each row is container of elements) assuming that right side consists if one followed by zeroes
//empty result means "no solution"
template <class Field, class ElementMatrixContainer, class CombinedField = ExactFieldAsCombined<Field>>
void SolveWithRightSideContainigSingleOne(const Field& field, const ElementMatrixContainer& matrix, std::vector<Element<Field>>& result, int max_diferent_numbers_in_coefficients)
{
	result.clear();
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
