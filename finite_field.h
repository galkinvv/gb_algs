#pragma once
template <class ZField>
struct FiniteField
{
	typedef typename ZField::Value Value;
	FiniteField(Value mod)
		: mod_(mod)
	{}
private:
	Value mod_;
};