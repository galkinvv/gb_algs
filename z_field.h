#pragma once
#include <cstdint>
template <class Integer>
struct ZField
{
	typedef Integer Value;
	bool IsZField()const
	{
		return true;
	}
	bool IsIdentity(const Value& value)const
	{
		return value == 1;
	}
};
typedef ZField<std::uint8_t> ZField8;
typedef ZField<std::uint16_t> ZField16;
typedef ZField<std::uint32_t> ZField32;
typedef ZField<std::uint64_t> ZField64;
