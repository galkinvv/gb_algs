#pragma once
#include <cstdint>
#include <limits>
#include <cassert>
template <class Integer>
struct ZRing
{
	class Value
	{
	  public:
		Value():
			i()
		{}
		template <class Integer2>
		Value Import(const Integer2& ext)
		{
			i = ext;
			assert(Integer2(i) == ext);
			return *this;
		}
		template <class Integer2>
		Integer2 Export()const
		{
			Integer2 result = i;
			assert(Integer(result) == i);
			return result;
		}
	  private:
		friend struct ZRing;
		Integer i;
	};
	
	struct DivResult
	{
		Value quot;
		Value rem;
	};
	template <class F>
	void SetRandom(F random_functor, Value& value)
	{
		const int single_call_bytes = sizeof(random_functor());
		int random_bytes = 0;
		value.i = 0;
		for(;;)
		{
			value.i ^= Integer(random_functor());
			random_bytes += single_call_bytes;
			if (random_bytes >= sizeof(value.i))
			{
				break;
			}
		}
	}
	
	static_assert(std::is_unsigned<Integer>::value, "Integer must be unsigned for use with ZRing");
	void Divide(const Value& numerator, const Value& denominator, DivResult& result)
	{
		//division never overflows for unsigned integers
		result.quot.i = divident/divider;
		result.rem.i = divident%divider;
	}
	
	void Subtract(const Value& from, const Value& what, const DivResult& multiplier, Value& result)
	{
		assert(std::numeric_limits<Integer>::max() / what.i > multiplier.quot.i); //check for mul overflow
		result.i = from.i - what.i * multiplier.quot.i;
		assert(result.i <= from.i)//check for subtract overflow
	}

	void SetZero(Value& result)const 
	{
		result.i = 0;
	}
	void SetOne(Value& result)const 
	{
		result.i = 1;		
	}
	bool IsZero(const Value& result)const 
	{
		return 0 == result.i;
	}
};
typedef ZRing<std::uint8_t> ZRing8;
typedef ZRing<std::uint16_t> ZRing16;
typedef ZRing<std::uint32_t> ZRing32;
