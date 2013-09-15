#pragma once
#include "types.h"
namespace F4MPI
{
	enum class FieldType
	{
		Z, Approx
	};
	struct IOPolynomSet
	{
		PolynomSet polys;
		FieldType type;
		int field_char;
		CMonomialBase::Order mon_order;
	};
}