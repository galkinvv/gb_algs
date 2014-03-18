#pragma once
#include <exception>
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
	bool IsZero(const Field& field, const typename Field::Value& value)
	{
		typename Field::Value one, zero, unused_result;
		typename Field::DivResult zero_div_result;
		field.SetOne(one);
		field.SetZero(zero);
		field.Divide(zero, one, zero_div_result);
		return ExactSubtractionResultInfo::Zero ==  field.Subtract(value, zero, zero_div_result, unused_result);
	}
}