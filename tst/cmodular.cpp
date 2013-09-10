#include <gtest/gtest.h>
#include "cmodular.h"
using namespace F4MPI;
TEST(CModular, eq)
{
	CModular m1(1000);
	CModular m2(1000);
	EXPECT_EQ(m1,m2);
}