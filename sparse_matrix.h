#pragma once
#include "finite_field.h"
#include "z_ring.h"
#include "utils.h"
#include <vector>

namespace SparseMatrix
{
template <class Field>
struct Element
{
	typename Field::Value value;
	int column;
};

namespace{
	template <class Field, class ExactField>
	struct PairElement
	{
		typename Field::Value approx;
		typename ExactField::Value exact;
		int column;
	};	
	
	template <class Field, class ExactField>
	struct PairRow
	{
		std::vector<PairElement<Field, ExactField>> left;
		PairElement<Field, ExactField> right;
	};
}

//solves matrix equation in field with matrix as left side(Container of rows, where each row is container of elements) assuming that right side consists if one followed by zeroes
//empty result means "no solution"

template <class Field, class ElementMatrixContainer, class ExactField = FiniteField<ZPlusRing32>>
void SolveWithRightSideContainigSingleOne(const Field& field, const ElementMatrixContainer& matrix, std::vector<Element<Field>>& result, int max_diferent_numbers_in_coefficients)
{
	//calculate P_0 based on matrix.size() and max_diferent_numbers_in_coefficients
	//calculate N_0 based on P_0. 
	//select field based on P_0 and N_0
	auto fast_field = ExactField::CreateZpFieldWithChar(2);
	std::vector<PairRow<Field, ExactField>> pair_matrix;
	pair_matrix.reserve(matrix.size());
	for (const auto& row: matrix)
	{
		auto& pair_row = emplaced_back(pair_matrix);
	}
	result.clear();
	//TODO
}
}
