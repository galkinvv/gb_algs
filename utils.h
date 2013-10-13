#pragma once
#include <type_traits>
#include <utility>
#include <initializer_list>
class NoCopy
{
	void operator=(const NoCopy&);
	NoCopy(const NoCopy&);
protected:
	NoCopy(){}
};

template <class... Args> void IgnoreIfUnused(Args...){}

template < typename T, size_t N >
constexpr int countof( T ( & /*arr*/ )[ N ] )
{
	static_assert(std::extent< T[ N ] >::value <= size_t(std::numeric_limits<int>::max()), "Static array size exceeds int. countof will not work");
	return int(std::extent< T[ N ] >::value);
}

#define DECLARE_FORWARDING_CONSTRUCTOR(Self, forward_to)\
		template <class... TA> Self(TA&&... args): forward_to(std::forward<TA>(args)...) {}


template <class T>
class PseudoPointer
{
	typedef typename std::remove_reference<T>::type TValue;
	TValue temp_data_;
	
  public:
	DECLARE_FORWARDING_CONSTRUCTOR(PseudoPointer, temp_data_)
	
	TValue* operator->()
	{
		return &temp_data_;
	}
};

template <class T>
const std::initializer_list<T>& ilist(const std::initializer_list<T>&  list)
{
	return list;
}
