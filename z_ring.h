#pragma once
#include <cstdint>
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
	void SetZero(Value& result)
	{
		result.i = 0;
	}
	void SetOne(Value& result)
	{
		result.i = 1;		
	}
};
typedef ZRing<std::uint8_t> ZRing8;
typedef ZRing<std::uint16_t> ZRing16;
typedef ZRing<std::uint32_t> ZRing32;
