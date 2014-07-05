#pragma once
#include <cassert>
#include "field_base.h"
#include "utils.h"


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
		friend std::ostream& operator <<(std::ostream &output, const Value &x) {
			return output << static_cast<ZValue>(x);
        }

		friend class FiniteField;
	};

	class DivResult
	{
		Value negatedValue; 	//stores negative result of division
		friend class FiniteField;
	};
	
	template <class F>
	void SetRandom(F random_functor, Value& value) const{
		z_.SetRandom(random_functor, value);
		Normalize(value);
	}

	void DivByOne(const Value& divident, DivResult& result) const{
		//simple cast is not enough, because DivResult is stored as negated value

		//zero case is specially handled to get rid of normalize call
		if (z_.IsZero(divident))
		{
			z_.SetZero(result.negatedValue);
		}
		else
		{
			z_.SubtractFromNotSmaller(mod_, divident, result.negatedValue);
		}
		assert(IsNormalized(result.negatedValue));
	}

	void Divide(const Value& divident, const Value& divider, DivResult& result) const{
		//first solve mod_*x_unused - divider * divider_inverse  = -1, which is equivalent to inverting in Z_mod field
		ZValue x_unused;
		ZValue negative_divider_inverse;
		ExtendedEuclid(mod_, divider, x_unused, negative_divider_inverse, SignedOne::PlusOne); //SignedOne::PlusOne is used because result should be stored as negative
		/*
		not easy to check because of negation
		 assert([&]{
			ZDivResult expected_minus_one_remainder;
			ZValue expected_minus_one_by_mod;
			ZValue expected_one;
			z_.Mul(divider_inverse, divider, expected_minus_one_by_mod);
			z_.Divide(expected_one_by_mod, mod_, expected_minus_one_remainder);
			z_.
			return z_.IsOne(expected_one_remainder.rem);
			}());
		*/
		z_.MulMod(negative_divider_inverse, divident, mod_, result.negatedValue);
		assert(IsNormalized(result.negatedValue));
	}

	ExactSubtractionResultInfo Subtract(const Value& from, const Value& what, const DivResult& multiplier, Value& result) const{
		//DivResult contanes already negated value, so perform addition
		z_.AddMulMod(from, what, multiplier.negatedValue, mod_, result);
		assert(IsNormalized(result));
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

	void SetZero(Value& result)const{
		z_.SetZero(result);
	}
	
	void SetOne(Value& result) const{
		z_.SetOne(result);
	}

	template <class Integer>
	static FiniteField CreateZpFieldWithChar(const Integer& value) {
		FiniteField result;
		result.z_.Import(result.mod_, value);
		assert(RingHelpers::IsPrime(result.z_, result.mod_));
		//result.mod_ must be prime
		return result;
	}

	template <class Integer>
	Integer ExportZpModulus()const {
		return z_.template Export<Integer>(mod_);
	}
	
	friend std::ostream& operator <<(std::ostream &output, const FiniteField &f) {
		return output << "FiniteField_MOD_" << f.mod_ ;
	}

private:
	FiniteField() {}

	void Normalize(ZValue& value)const {
		ZDivResult byModDivResult;
		z_.Divide(value, mod_, byModDivResult);
		value = byModDivResult.rem;
		assert(IsNormalized(value));
	}

	bool IsNormalized(ZValue& value)const {
		ZValue zero;
		z_.SetZero(zero);
		return z_.Less(value, mod_) && !z_.Less(value, zero);
	}

	SignedOne NegateOne(SignedOne signedOne) const{
		switch(signedOne) {
		case SignedOne::PlusOne:
			return SignedOne::MinusOne;
		case SignedOne::MinusOne:
			return SignedOne::PlusOne;
		default:
			assert(!"bad value for signedOne");
			return SignedOne::PlusOne; //only for return statement checkers or compilation with NDEBUG
		}
	}
	
	//epects that a > b
	//Finds x, y (0 < x <= b+1) and (0 <= y < a) for equation a*x - b*y = signedOne using calculation wih only positive numbers.
	// b < a; gcd(a, b) == 1 are required preconditions
	void ExtendedEuclid(const ZValue& a, const ZValue& b, ZValue& x, ZValue& y, SignedOne signedOne) const{
		assert(z_.Less(b, a));
		bool finalStep = false;
		//stop recursion only with positive right side
		if (signedOne == SignedOne::PlusOne)
		{
			if (z_.IsZero(b))
			{
				assert(z_.IsOne(a)); //fail of this assertion means that a and b has non-one gcd.
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
			z_.AddMul(y, x, abDivResult.quot, y);
		}
		assert([&]{
			if (z_.IsZero(x))
			{
				return true;
			}
			ZValue xMinusOne{x};
			z_.Decrement(xMinusOne);
			return !z_.Less(b ,xMinusOne);
		}());
		assert(z_.Less(y, a));
	}

	typedef typename ZPlusRing::DivResult ZDivResult;
	ZValue mod_;
	ZPlusRing z_;
};
