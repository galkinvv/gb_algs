#pragma once
#include "field_base.h"
#include "utils.h"

template <class ExactField, class ApproxField>
class CombinedField
{
	class Value
	{
		friend class CombinedField;
		typename ExactField::Value exact;
		typename ApproxField::Value approx;
	};
	
	class DivResult
	{
		friend class CombinedField;
		typename ExactField::DivResult exact;
		typename ApproxField::DivResult approx;
	};

	template <class F>
	void ExtendWithRandom(const typename ApproxField::Value& approx_value, F random_functor, Value& value)const
	{
		value.approx = approx_value;
		ef_.SetRandom(random_functor, value.exact);
	}

	void ExtractApprox(const Value& value, typename ApproxField::Value& approx_value)const
	{
		approx_value = value.approx;
	}

	void DivByOne(const Value& divident, DivResult& result) const{
		af_.DivByOne(divident.approx, result.approx);
		ef_.DivByOne(divident.exact, result.exact);
	}

	void Divide(const Value& numerator, const Value& denominator, DivResult& result)const
	{
		af_.Divide(numerator.approx, denominator.approx, result.approx); //may  throw unexact_divisor_exception
		ef_.Divide(numerator.exact, denominator.exact, result.exact);
	}
	
	bool IsPreciserDivisor(const Value& divisor0, const Value& divisor1)const
	{
		return af_.IsPreciserDivisor(divisor0.approx, divisor1.approx);
	}
	
	void SetZero(Value& result)const
	{
		af_.SetZero(result.approx);
		ef_.SetZero(result.exact);
	}
	
	void SetOne(Value& result)const
	{
		af_.SetOne(result.approx);
		ef_.SetOne(result.exact);
	}
	
	ExactSubtractionResultInfo Subtract(const Value& from, const Value& what, const DivResult& multiplier, Value& result)const
	{
		auto approx_result = af_.Subtract(from.approx, what.approx, multiplier.approx, result.approx);
		auto exact_result = ef_.Subtract(from.exact, what.exact, multiplier.exact, result.exact);
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

template <class ExactField>
class ExactFieldAsCombined
{
  public:
	ExactFieldAsCombined(const ExactField& other):ef_(other){}

	class Value
	{
		friend class ExactFieldAsCombined;
		typename ExactField::Value exact;
	};
	
	class DivResult
	{
		friend class ExactFieldAsCombined;
		typename ExactField::DivResult exact;
	};

	template <class F>
	void ExtendWithRandom(const typename ExactField::Value& approx_value, F random_functor, Value& value)const
	{
		value.exact = approx_value;
		IgnoreIfUnused(random_functor);
	}
	void ExtractApprox(const Value& value, typename ExactField::Value& approx_value)const
	{
		approx_value = value.exact;
	}

	void DivByOne(const Value& divident, DivResult& result) const{
		ef_.DivByOne(divident.exact, result.exact);
	}

	void Divide(const Value& numerator, const Value& denominator, DivResult& result)const
	{
		ef_.Divide(numerator.exact, denominator.exact, result.exact);
	}
	bool IsPreciserDivisor(const Value& divisor0, const Value& divisor1)const
	{
		//mostly times is uncomparable
		IgnoreIfUnused(divisor0);
		return FieldHelpers::IsZero(ef_, divisor1.exact);
	}
	
	void SetZero(Value& result)const
	{
		ef_.SetZero(result.exact);
	}
	
	void SetOne(Value& result)const
	{
		ef_.SetOne(result.exact);
	}
	
	ExactSubtractionResultInfo Subtract(const Value& from, const Value& what, const DivResult& multiplier, Value& result)const
	{
		return ef_.Subtract(from.exact, what.exact, multiplier.exact, result.exact);
	}
private:
  ExactField ef_;
};
