#pragma once
template <class ZField>
struct FiniteField
{
	typedef typename ZField::Value Value;
	FiniteField(Value mod)
		: mod_(mod)
	{}
	bool IsFiniteZpFieldWithChar(Value mod)
	{
		return mod == mod_;
	}

private:
	Value mod_;
};