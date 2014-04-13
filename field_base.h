#pragma once
#include <exception>
#include <gmpxx.h>

enum class ExactSubtractionResultInfo
{
	Zero,
	NonZero
};
enum class ApproxSubtractionResultInfo
{
	MaybeZero,
	Nonzero
};
struct unexact_divisor_exception: std::exception{};
struct cant_detect_zero_equality_exception: std::exception{};

namespace FieldHelpers
{
	template <class Field>
	typename Field::Value Zero(const Field& field)
	{
		typename Field::Value result;
		field.SetZero(result);
		return result;
	}
	
	template <class Field>
	typename Field::Value One(const Field& field)
	{
		typename Field::Value result;
		field.SetOne(result);
		return result;
	}
	
	template <class Field>
	typename Field::DivResult DivByOne(const Field& field, const typename Field::Value& value)
	{
		//TODO:: add "DivideByOne" method
		typename Field::DivResult result;
		field.Divide(value, One(field), result);
		return result;
	}
	
	template <class Field>
	typename Field::Value FracAsValue(const Field& field, const typename Field::DivResult& value)
	{
		typename Field::Value negated_result, result;
		field.Subtract(Zero(field), One(field), value, negated_result);
		field.Subtract(Zero(field), One(field), DivByOne(field, negated_result), result);
		return result;
	}
	
	template <class Field>
	bool IsEqual(const Field& field, const typename Field::Value& value0, const typename Field::Value& value1)
	{
		typename Field::Value unused_result;
		return ExactSubtractionResultInfo::Zero ==  field.Subtract(value0, value1, DivByOne(field, One(field)), unused_result);
	}	

	template <class Field>
	bool IsZero(const Field& field, const typename Field::Value& value)
	{
		return IsEqual(field, value, Zero(field));
	}	

}

namespace RingHelpers{
	template <class Ring>
	bool IsPrime(const Ring& ring, const typename Ring::Value& value)
	{
		auto value_as_mpz = ring.template Export<mpz_class>(value);
		static const int kDefinitelyComposite = 0;
		static const int kSuitableNumberOfChecks = 25;
		return mpz_probab_prime_p(value_as_mpz.get_mpz_t(), kSuitableNumberOfChecks) != kDefinitelyComposite;
	}	
}