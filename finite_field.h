#pragma once
#include <cassert>
#include "field_base.h"

template <class ZPlusRing>
struct FiniteField
{
	typedef typename ZPlusRing::Value ModulusValue;
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
	void Divide(const Value& divident, const Value& divider, DivResult& result)
	{
		
	}
	ExactSubtractionResultInfo Subtract(const Value& from, const Value& what, const DivResult& multiplier, Value& result)
	{
		return ExactSubtractionResultInfo::Zero;
	}

	template <class Integer>
	void Import(const Integer& value, Value& result)const
	{
		z_.Import(result, value);
		Normalize(result);
	}
	
	template <class Integer>
	Integer Export(const Value& value)const
	{
		return z_.template Export<Integer>(value);
	}

	void SetZero(Value& result)
	{
		z_.SetOne(result);
	}
	void SetOne(Value& result)
	{
		z_.SetOne(result);		
	}
	
	template <class Integer>
	static FiniteField CreateZpFieldWithChar(const Integer& value)
	{
		FiniteField result;
		result.z_.Import(result.mod_, value);
		//result.mod_ must be prime
		return result;
	}
private:
	FiniteField(){}

	void Normalize(Value& value)const
	{
		//TODO:add division
	}
	
	//Finds non-both-zero x and y for equation a*x = b*y using calculation wih only positive numbers
	void ExtendedEuclid(const Value& a, const Value& b, Value& x, Value& y)
	{
		if (z_.Less(a, b))
		{
			ExtendedEuclid(b, a, y, x);
		}
		else if (z_.IsZero(b))
		{
			z_.SetZero(x);
			z_.SetOne(y);
		}
		else
		{            
			typename ZPlusRing::DivResult abDivResult;
			z_.Divide(a, b, abDivResult);
			ExtendedEuclid(b, abDivResult.rem, y, x);
			z_.Add(CopyValue(y), x, abDivResult.qout, y);
		}
	}
	ModulusValue mod_;
	ZPlusRing z_;
};
