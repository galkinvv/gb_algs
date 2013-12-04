#pragma once
#include <cassert>
template <class ZField>
struct FiniteField
{
	typedef typename ZField::Value Value;
	FiniteField(Value mod)
		: mod_(mod)
	{}
	bool IsFiniteZpFieldWithChar(Value mod)const
	{
		return mod == mod_;
	}
	bool IsIdentity(const Value& value)const
	{
		AssertNormalized(value);
		return value == 1;
	}

private:
	void AssertNormalized(const Value& value)const
	{
		assert(value >= 0);
		assert(value < mod_);
	}
	Value mod_;
};
