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
};
typedef ZField<std::uint8_t> ZField8;
typedef ZField<std::uint16_t> ZField16;
typedef ZField<std::uint32_t> ZField32;
typedef ZField<std::uint64_t> ZField64;