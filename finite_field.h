#pragma once
#include <cassert>
#include "field_base.h"

template <class ZRing>
struct FiniteField
{
	typedef typename ZRing::Value ModulusValue;
	struct Value: private ModulusValue
	{
		friend class FiniteField;
	};
	typedef Value DivResult;
	
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
	void Divide(const Value& numerator, const Value& denominator, DivResult& result)
	{
		
	}
	ExactSubtractionResultInfo Subtract(const Value& from, const Value& what, const DivResult& multiplier, Value& result)
	{
		return ExactSubtractionResultInfo::Zero;
	}
	template <class Integer>
	void Import(const Integer& value, Value& result)const
	{
		result.Import(value);
		Normalize(result);
	}
	template <class Integer>
	Integer Export(const Value& value)const
	{
		return value.template Export<Integer>();
	}

	void SetZero(Value& result)
	{
		z_.SetOne(result);
	}
	void SetOne(Value& result)
	{
		z_.SetOne(result);		
	}
	static FiniteField CreateZpFieldWithChar(const ModulusValue& mod)
	{
		FiniteField result;
		//mod must be prime
		result.mod_ = mod;
		return result;
	}
private:
	FiniteField(){}

	void Normalize(Value& value)const
	{
		//TODO:add division
	}
	 ModulusValue mod_;
	ZRing z_;
};
