#pragma once
#include <cstdint>
#include <limits>
#include <cassert>
template <class Integer>
struct ZPlusRing
{
	class Value
	{
	  public:
		Value():
			i()
		{}
	  private:
		friend struct ZPlusRing;
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
	
	static_assert(std::is_unsigned<Integer>::value, "Integer must be unsigned for use with ZPlusRing");
	void Divide(const Value& divident, const Value& divider, DivResult& result)
	{
		//division never overflows for unsigned integers and always is positive
		result.quot.i = divident/divider;
		result.rem.i = divident%divider;
	}
	
	void Add(const Value& from, const Value& what, const Value& multiplier, Value& result)
	{
		assert(std::numeric_limits<Integer>::max() / what.i > multiplier.i); //check for mul overflow
		result.i = from.i + what.i * multiplier.i;
		assert(result.i >= from.i);//check for add overflow
	}

	void SetZero(Value& result)const 
	{
		result.i = 0;
	}
	void SetOne(Value& result)const 
	{
		result.i = 1;		
	}
	
	bool Less(const Value& v1, const Value& v2)const
	{
		return v1.i < v2.i;
	}
	
	bool IsZero(const Value& result)const 
	{
		return 0 == result.i;
	}
	
	template <class Integer2>
	void Import(Value& value, const Integer2& ext)const
	{
		value.i = ext;
		assert(Integer2(value.i) == ext);
	}
	
	template <class Integer2 = Integer>
	Integer2 Export(const Value& value)const
	{
		Integer2 result = value.i;
		assert(Integer(result) == value.i);
		return result;
	}
	
	typedef std::numeric_limits<Integer> numeric_limits;
};
typedef ZPlusRing<std::uint8_t> ZPlusRing8;
typedef ZPlusRing<std::uint16_t> ZPlusRing16;
typedef ZPlusRing<std::uint32_t> ZPlusRing32;
