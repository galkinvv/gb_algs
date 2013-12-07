#pragma once
#include <cassert>
#include "field_base.h"

template <class ZRing>
struct FiniteField
{
	typedef typename ZRing::Value Value;
	typedef Value Frac;
	
	explicit FiniteField(Value mod)
		: mod_(mod)
	{}
	bool IsFiniteZpFieldWithChar(Value mod)const
	{
		return mod == mod_;
	}
	template <class F>
	void SetRandom(F random_functor, Value& value)
	{
		z_.SetRandom(random_functor, value);
		Normalize(value);
	}
	void Divide(const Value& numerator, const Value& denominator, Frac& result)
	{
		
	}
	ExactSubtractionResultInfo Subtract(const Value& from, const Value& what, const Frac& multiplier, Value& result)
	{
		return ExactSubtractionResultInfo::Zero;
	}
	void SetZero(Value& result)
	{
		z_.SetOne(result);
	}
	void SetOne(Value& result)
	{
		z_.SetOne(result);		
	}
private:
	void Normalize(Value& value)
	{
		
	}
	Value mod_;
	ZRing z_;
};
