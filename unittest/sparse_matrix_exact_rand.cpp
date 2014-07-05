#include <limits>
#include "z_ring.h"
#include "finite_field.h"
#include "combined_field.h"
#include "sparse_matrix_base.h"

TEST(SparseMatrixExactRandom, Z2determined)
{
	typedef FiniteFieldParam<ZPlusRing32>::Module<2> Param;
	Param::Matrix m;
	const auto E = [&](int column, unsigned long long value){return Param::Element::FromCV(column, FieldHelpers::Imp(m.field, value));};
	//tests for big randdom matrices based on https://cloud.sagemath.com/projects/16ae0d9e-83b0-46c0-9630-d8d4469ec367/files/solvable%20matrices.sagews
	{//9x9 matrix with det 1
		m.Clear();
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddElement(4, 1u);
		m.AddElement(6, 1u);
		m.AddElement(8, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddElement(4, 1u);
		m.AddElement(5, 1u);
		m.AddElement(7, 1u);
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddElement(3, 1u);
		m.AddElement(4, 1u);
		m.AddElement(5, 1u);
		m.AddElement(6, 1u);
		m.AddElement(8, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddElement(3, 1u);
		m.AddElement(4, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(2, 1u);
		m.AddElement(3, 1u);
		m.AddElement(4, 1u);
		m.AddElement(5, 1u);
		m.AddElement(6, 1u);
		m.AddElement(8, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(3, 1u);
		m.AddElement(6, 1u);
		m.AddElement(8, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(2, 1u);
		m.AddElement(3, 1u);
		m.AddElement(5, 1u);
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddElement(4, 1u);
		m.AddElement(8, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(3, 1u);
		m.AddElement(6, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(0, 1u), E(1, 1u), E(3, 1u), E(4, 1u), E(6, 1u), E(7, 1u)}));
	}
	{//the same matrix with to identical additional colums and rows (on the top)
		m.Clear();
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(10, 1u);
		m.AddRow();
		m.AddElement(2, 1u);
		m.AddElement(5, 1u);
		m.AddElement(7, 1u);
		m.AddElement(9, 1u);
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddElement(3, 1u);
		m.AddElement(5, 1u);
		m.AddElement(6, 1u);
		m.AddElement(8, 1u);
		m.AddRow();
		m.AddElement(2, 1u);
		m.AddElement(4, 1u);
		m.AddElement(5, 1u);
		m.AddElement(6, 1u);
		m.AddElement(7, 1u);
		m.AddElement(9, 1u);
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddElement(3, 1u);
		m.AddElement(4, 1u);
		m.AddElement(5, 1u);
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddElement(3, 1u);
		m.AddElement(4, 1u);
		m.AddElement(5, 1u);
		m.AddElement(6, 1u);
		m.AddElement(7, 1u);
		m.AddElement(9, 1u);
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddElement(4, 1u);
		m.AddElement(7, 1u);
		m.AddElement(9, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(3, 1u);
		m.AddElement(4, 1u);
		m.AddElement(6, 1u);
		m.AddElement(10, 1u);
		m.AddRow();
		m.AddElement(2, 1u);
		m.AddElement(3, 1u);
		m.AddElement(5, 1u);
		m.AddElement(9, 1u);
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddElement(4, 1u);
		m.AddElement(7, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(10, 1u);
		ExpectNoSolution(m);
	}
}
