#pragma once
#include <cstdint>
#include <limits>
#include <cassert>
#include <type_traits>
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
	void SetRandom(F random_functor, Value& value)const
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

	void Divide(const Value& divident, const Value& divider, DivResult& result)const
	{
		//division never overflows for unsigned integers and always is positive
		result.quot.i = divident.i/divider.i;
		result.rem.i = divident.i%divider.i;
	}
	
	void Add(const Value& from, const Value& what, const Value& multiplier, Value& result)const
	{
		Mul(what, multiplier, result);
		result.i = from.i + result.i;
		assert(result.i >= from.i);//check for add overflow
	}

	void Mul(const Value& mult0, const Value& mult1,  Value& result)const
	{
		if (mult0.i)
		{
			assert(numeric_limits::max() / mult0.i >= mult1.i); //check for mul overflow
		}
		result.i =mult0.i * mult1.i;
	}

	void SetZero(Value& result)const 
	{
		result.i = 0;
	}
	
	void SetOne(Value& result)const 
	{
		result.i = 1;		
	}
	
	void Decrement(Value& result)const 
	{
		assert(result.i > 0);
		--result.i;
	}

	bool Less(const Value& v1, const Value& v2)const
	{
		return v1.i < v2.i;
	}
	
	bool IsZero(const Value& result)const 
	{
		return 0 == result.i;
	}
	
	bool IsOne(const Value& result)const 
	{
		return 1 == result.i;
	}
	template <class Integer2>
	void Import(Value& value, const Integer2& ext)const
	{
		assert(ext >= numeric_limits::min());
		assert(ext <= numeric_limits::max());
		value.i = ext;
		assert(Integer2(value.i) == ext);
	}
	
	template <class Integer2 = Integer>
	Integer2 Export(const Value& value)const
	{
		typedef std::numeric_limits<Integer2> numeric_limits2;
		assert(value.i >= numeric_limits2::min());
		assert(value.i <= numeric_limits2::max());
		Integer2 result = value.i;
		assert(Integer(result) == value.i);
		return result;
	}
	
	typedef std::numeric_limits<Integer> numeric_limits;
};
typedef ZPlusRing<std::uint8_t> ZPlusRing8;
typedef ZPlusRing<std::uint16_t> ZPlusRing16;
typedef ZPlusRing<std::uint32_t> ZPlusRing32;
