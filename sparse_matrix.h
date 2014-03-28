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
		typename Field::Value right;
	};
	
	//triangulates matrix. lead_columns - filled with column numbers treated as leading.
	//i-th number in lead_columns corresponds to i-th row in matrix
	//lead elements of triangulated matrix are NOT set to exact ones - they are arbitrary values

	const int kRowIsCompletlyZero = -1; //kRowIsCompletlyZero in lead column means  completely zero row
	const int kUnitializedLeadColumn = -2; //used only internally
	struct incompatible_system_exception: std::exception{};
	//throws incompatible_system_exception when system can't be triangulated (too much non-zero rows for example)
	
	template <class Field>
	void SubMatrixRows(const Field& field, const RowWithRightPart<Field>& orig, const RowWithRightPart<Field>& modifier, RowWithRightPart<Field>& result, typename std::vector<Element<Field>>::const_iterator lead_orig, typename std::vector<Element<Field>>::const_iterator lead_modifier)
	{
		typename Field::Value zero;
		field.SetZero(zero);
		typename Field::DivResult row_multiplier;
		assert(lead_orig->column == lead_modifier->column);
		field.Divide(lead_orig->value, lead_modifier->value, row_multiplier);
		result.left.reserve(orig.left.size() + modifier.left.size());
		field.Subtract(orig.right, modifier.right, row_multiplier, result.right);
		auto io = orig.left.begin();
		auto io_end = orig.left.end();
		auto im_end = modifier.left.end();
		//TODO: adapt
		/*
		while(i1!=i1f && i2!=i2f){			
			int dif = p1.getMon(i1).compareTo(p2.getMon(i2));
			CModular c;
			if(dif==0){				
				c = p1.getCoeff(i1) +addCoeff*p2.getCoeff(i2);
				M = p1.getMon(i1);
				++i1;
				++i2;
			}		
			else if(dif>0){
				c = p1.getCoeff(i1);
				M = p1.getMon(i1);				
				++i1;
			}		
			else{
				c = addCoeff*p2.getCoeff(i2);
				M = p2.getMon(i2);				
				++i2;
			}
			if(c!=0){
				Result.pushTermBack(c, M);
			}
		}
		while(i1!=i1f){			
			CModular c;
			c = p1.getCoeff(i1);
			M = p1.getMon(i1);				
			if(c!=0){
				Result.pushTermBack(c, M);
			}
			++i1;
		}

		while(i2!=i2f){			
			CModular c;
			c = addCoeff*p2.getCoeff(i2);
			M = p2.getMon(i2);
			if(c!=0){
				Result.pushTermBack(c, M);
			}
			++i2;
		}	
		*/
	}
	
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
			
			auto lead_item_in_min_row = *min_min_row_it;
			const int min_row_index = min_min_row_it - min_row_iterators.begin();
			
			const auto found_min_row = matrix.begin()+min_row_index;
			
			--used_rows;
			const auto last = matrix.begin()+used_rows;
			if (found_min_row != last)
			{
				found_min_row->swap(*last);
			}
			assert(matrix[used_rows]->left.empty()); //should be empty because we positioned in in a such way
			const int lead_column = lead_item_in_min_row->column;
			lead_columns[used_rows] = lead_column;
			const auto& last_row = *last;
			assert((lead_item_in_min_row - last_row->left.begin()) < last_row->left.size());

			//go to subtraction without  normalizing the  lead coef in last
			for (ir=matrix.begin(); ir != last; ++ir)
			{
				auto lead_row = std::lower_bound((**ir).left.begin(),(**ir).left.end(), lead_column, [](const Element<Field>& el, int column){return el.column < column;});
				if (lead_row == (**ir).left.end() || lead_row->column != lead_column)
				{
					continue;
				}
				auto_unique_ptr<RowWithRightPart<Field>> modified_row;
				SubMatrixRows(field, **ir, *last_row, *modified_row, lead_row, lead_item_in_min_row);
				ir->reset(modified_row.release());
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
			combined_field.SetOne(pair_row->right);
		}
		else
		{
			//all other rows - zeroes
			combined_field.SetZero(pair_row->right);
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
