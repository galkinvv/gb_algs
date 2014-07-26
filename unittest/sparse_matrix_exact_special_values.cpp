#include <limits>
#include "z_ring.h"
#include "finite_field.h"
#include "combined_field.h"
#include "sparse_matrix_base.h"

TEST(SparseMatrixSpecialValues, Z2determined)
{
	typedef FiniteFieldParam<ZPlusRing32>::Module<2> Param;
	Param::Matrix m;
	const auto E = m.ElementConstructor();
	{//zero matrix size 1
		m.Clear();
		m.AddRow();
		ExpectNoSolution(m);
	}
	{//zero matrix size 3
		m.Clear();
		m.AddRow();
		m.AddRow();
		m.AddRow();
		ExpectNoSolution(m);
	}
	{//identity matrix size 1
		m.Clear();
		m.AddRow();
		m.AddElement(0, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(0, 1u)}));
	}
	{//identity matrix size 3
		m.Clear();
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddRow();
		m.AddElement(2, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(0, 1u)}));
	}
	{//matrix with only last column needed
		m.Clear();
		m.AddRow();
		m.AddElement(2, 1u);
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(2, 1u)}));
	}
	{//matrix with big column number
		m.Clear();
		m.AddRow();
		m.AddElement(42, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(42, 1u)}));
	}
	{//matrix with single zero column
		m.Clear();
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		ExpectNoSolution(m);
	}
	{//matrix with single columnwith ones
		m.Clear();
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		ExpectNoSolution(m);
	}
	{//matrix with single columnw eqaul to result, non-obviousely initialized
		m.Clear();
		m.AddRow();
		m.AddElement(0, 3u);
		m.AddRow();
		m.AddElement(0, 2u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(0, 1u)}));
	}
	{//matrix with every column needed in sum
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
	{//matrix with every column needed in sum - other column order
		m.Clear();
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(3, 1u);
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddElement(3, 1u);
		m.AddRow();
		m.AddElement(2, 1u);
		m.AddElement(3, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(0, 1u), E(1, 1u), E(2, 1u), E(3, 1u)}));
	}
	{//matrix with NOT every column needed in sum
		m.Clear();
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(3, 1u);
		m.AddRow();
		m.AddElement(3, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(2, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddElement(3, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(0, 1u), E(2, 1u)}));
	}
	{//matrix with NOT every column needed in sum - other column order
		m.Clear();
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(2, 1u);
		m.AddElement(3, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddElement(3, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddElement(3, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(1, 1u), E(3, 1u)}));
	}
	{//wide matrix with zero columns; wide matrix with non-zero columns would have non-determined result, so is not included in this test
		m.Clear();
		m.AddRow();
		m.AddElement(42, 1u);
		m.AddRow();
		m.AddElement(12, 1u);
		m.AddElement(32, 1u);
		m.AddElement(42, 1u);
		m.AddRow();
		m.AddElement(2, 1u);
		m.AddElement(32, 1u);
		m.AddRow();
		m.AddElement(2, 1u);
		m.AddElement(12, 1u);
		m.AddElement(32, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(2, 1u), E(32, 1u), E(42, 1u)}));
	}
	{//tall matrix
		m.Clear();
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddRow();
		m.AddElement(1, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(2, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(3, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddElement(3, 1u);
		m.AddElement(4, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(3, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddRow();
		m.AddElement(0, 1u);
		m.AddElement(1, 1u);
		m.AddElement(2, 1u);
		m.AddElement(3, 1u);
		m.AddElement(4, 1u);
		EXPECT_PRED2(ExpectKnownSolution, m, ilist({E(0, 1u), E(2, 1u), E(3, 1u), E(4, 1u)}));
	}
}
