#pragma once
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <cassert>
#include <type_traits>

#include <utils.h>

//
template <class Integer, class BiggerInteger>
struct ZPlusRing
{
	class Value
	{
	  public:
		Value():
			i()
		{}
		friend std::ostream& operator <<(std::ostream &output, const Value &x) {
			return output << (IntegerForTextIO)x.i;
		}
	  private:
		typedef typename std::conditional<std::is_same<Integer, std::uint8_t>::value, std::uint32_t, Integer>::type IntegerForTextIO;
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
	static_assert(std::is_unsigned<BiggerInteger>::value, "BiggerInteger must be unsigned for use with ZPlusRing");
	static_assert(sizeof(BiggerInteger) >= 2*sizeof(Integer), "BiggerInteger must have at least 2 times more bits than Integer for use with ZPlusRing");

	void Divide(const Value& divident, const Value& divider, DivResult& result)const
	{
		//division never overflows for unsigned integers and always is positive
		result.quot.i = divident.i/divider.i;
		result.rem.i = divident.i%divider.i;
	}
	
	void AddMulMod(const Value& from, const Value& what, const Value& multiplier, const Value& mod, Value& result)const
	{
		MulMod(what, multiplier, mod, result);
		result.i = narrow_cast<Integer>((wide_cast<BiggerInteger>(from.i) + wide_cast<BiggerInteger>(result.i))%wide_cast<BiggerInteger>(mod.i));
	}

	void MulMod(const Value& mult0, const Value& mult1, const Value& mod,  Value& result)const
	{
		result.i = narrow_cast<Integer>((wide_cast<BiggerInteger>(mult0.i) * wide_cast<BiggerInteger>(mult1.i)) % wide_cast<BiggerInteger>(mod.i));
	}

	void AddMul(const Value& from, const Value& what, const Value& multiplier, Value& result)const
 	{
		if (what.i)
		{
			assert(numeric_limits::max() / what.i >= multiplier.i); //check for mul overflow
		}
		result.i =from.i + what.i * multiplier.i;
		assert(result.i >= from.i);//check for add overflow
 	}
 
	void SubtractFromNotSmaller(const Value& from, const Value& what, Value& result)const
 	{
		assert(from.i >= what.i);
		result.i = from.i - what.i;
		assert(result.i <= from.i);
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
		value.i = narrow_cast<Integer>(ext);
	}
	
	template <class Integer2 = Integer>
	Integer2 Export(const Value& value)const
	{
		typedef std::numeric_limits<Integer2> numeric_limits2;
		assert(value.i >= numeric_limits2::min());
		//gmp doesn't define non-zero max numeric limits
		assert(numeric_limits2::max() ==0 || value.i <= numeric_limits2::max());
		Integer2 result = value.i;
		assert(result == value.i);
		return result;
	}
	
	typedef std::numeric_limits<Integer> numeric_limits;
};
typedef ZPlusRing<std::uint_least8_t, std::uint_fast16_t> ZPlusRing8;
typedef ZPlusRing<std::uint_least16_t, std::uint_fast32_t> ZPlusRing16;
typedef ZPlusRing<std::uint_least32_t, std::uint_fast64_t> ZPlusRing32;
