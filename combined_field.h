#pragma once
template <class ExactField, class ApproxField>
class CombinedField
{
	class Value
	{
		friend class CombinedField;
		typename ExactField::Value exact;
		typename ApproxField::Value approx;
	};
	
	class Frac
	{
		friend class CombinedField;
		typename ExactField::Frac exact;
		typename ApproxField::Frac approx;
	};

	template <class F>
	void ExtendWithRandom(const ApproxField::Value& approx_value, F random_functor, Value& value)
	{
		value.approx = approx_value;
		ef_.SetRandom(random_functor, value.exact);
	}

	void ExtractApprox(const Value& value,  ApproxField::Value& approx_value)
	{
		approx_value = value.approx;
	}

	void Frac(const Value& numerator, const Value& denominator, Frac& result)
	{
		af_.Divide(numerator.approx, denominator.approx, result.approx); //may  throw unexact_divisor_exception
		ef_.Divide(numerator.exact, denominator.exact, result.exact);
	}
	bool IsPreciserDivisor(const Value& numerator, const Value& denominator)
	{
		return af_.IsPreciserDivisor(numerator.approx, denominator.approx);
	}
	
	ExactSubtractionResultInfo Subtract(const Value& from, const Value& what, const Frac& multiplier, Value& result)
	{
		auto approx_result = af_.Subtract(from.approx, what.approx, multiplier.approx, result.approx);
		auto exact_result = af_.Subtract(from.exact, what.exact, multiplier.exact, result.exact);
		if (exact_result == ExactSubtractionResultInfo::NonZero)
		{
			return ExactSubtractionResultInfo::NonZero;
		}
		if (approx_result != ApproxSubtractionResultInfo::MaybeZero)
		{
			throw cant_detect_zero_equality_exception();
		}
		return ExactSubtractionResultInfo::Zero;
	}

	ExactField ef_;
	ApproxField af_;
};