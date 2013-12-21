#pragma once
#include <cassert>
#include "field_base.h"

template <class ZPlusRing>
struct FiniteField {
  private:
    typedef typename ZPlusRing::Value ZValue;
	enum class SignedOne
	{
	    MinusOne,
	    PlusOne
	};

  public:
	class Value: ZValue
	{
		friend class FiniteField;
	};

	//stores negative result of division
	class DivResult:Value
	{
		friend class FiniteField;
	};
	
	template <class F>
	void SetRandom(F random_functor, Value& value) {
		z_.SetRandom(random_functor, value);
		Normalize(value);
	}

	void Divide(const Value& divident, const Value& divider, DivResult& result) {
		//first solve mod_*x - divider * result  = 1
		ZValue x_unused;
		ExtendedEuclid(mod_, divider, x_unused, result, SignedOne::PlusOne);
		z_.Mul(CopyValue(result), divident, result);
	}

	ExactSubtractionResultInfo Subtract(const Value& from, const Value& what, const DivResult& multiplier, Value& result) {
		//DivResult contanes already negated value, so perform addition
		z_.Add(from, what, multiplier, result);
		Normalize(result);
		if (z_.IsZero(result)) {
			return ExactSubtractionResultInfo::Zero;
		}
		return ExactSubtractionResultInfo::NonZero;
	}

	template <class Integer>
	void Import(const Integer& value, Value& result)const {
		z_.Import(result, value);
		Normalize(result);
	}

	template <class Integer>
	Integer Export(const Value& value)const {
		return z_.template Export<Integer>(value);
	}

	void SetZero(Value& result) {
		z_.SetOne(result);
	}
	void SetOne(Value& result) {
		z_.SetOne(result);
	}

	template <class Integer>
	static FiniteField CreateZpFieldWithChar(const Integer& value) {
		FiniteField result;
		result.z_.Import(result.mod_, value);
		//result.mod_ must be prime
		return result;
	}

	template <class Integer>
	Integer ExportZpModulus()const {
		return z_.template Export<Integer>(mod_);
	}

private:
	FiniteField() {}

	void Normalize(ZValue& value)const {
		ZDivResult byModDivResult;
		z_.Divide(value, mod_, byModDivResult);
		value = byModDivResult.rem;
	}

private:
	SignedOne NegateOne(SignedOne signedOne) {
		switch(signedOne) {
		case SignedOne::PlusOne:
			return SignedOne::MinusOne;
		case SignedOne::MinusOne:
			return SignedOne::PlusOne;
		}
	}
	//Finds x ,y ( 0 < x <= b) and ( 0 <= y < a) for equation a*x - b*y = signedOne using calculation wih only positive numbers.
	// b < a; gcd(a, b) == 1 are required preconditions
	void ExtendedEuclid(const ZValue& a, const ZValue& b, ZValue& x, ZValue& y, SignedOne signedOne) {
		assert(z_.Less(b, a));
		bool finalStep = false;
		//stop recursion only with positive right side
		if (signedOne == SignedOne::PlusOne)
		{
			if (z_.IsZero(b))
			{
				assert(z_.IsOne(a));
				//solve as 1*1 - 0*0 = 1
				z_.SetOne(x);
				z_.SetZero(y);
				finalStep = true;
			}
			else if (z_.IsOne(b))
			{
				//solve as a*1 - 1*(a-1) =1
				z_.SetOne(x);
				y = a;
				z_.Decrement(y);
				finalStep = true;
			}
		}
		if (!finalStep)
		{
			assert(!z_.IsZero(b));
			ZDivResult abDivResult;
			z_.Divide(a, b, abDivResult);
			//solve b*y1 - (a%b) * x1 = -signedOne (equivalent to (a%b) * x1 - b*y1 = signedOne)
			ExtendedEuclid(b, abDivResult.rem, y, x, NegateOne(signedOne));
			//original equation is solved as a * x1 - b*(y1 + a/b*x1) = signedOne
			z_.Add(CopyValue(y), x, abDivResult.qout, y);
		}
		assert(!z_.IsZero(x));
		assert(!z_.Less(b, x));
		assert(z_.Less(y, a));
	}

	typedef typename ZPlusRing::DivResult ZDivResult;
	ZValue mod_;
	ZPlusRing z_;
};
