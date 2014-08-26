#pragma once
#include <climits>
#include <limits>
#include <type_traits>
#include <utility>
#include <initializer_list>
#include <ostream>
#include <functional>
#include <memory>
#include <random>
#include <cassert>
#include <iostream>


class NoCopy
{
	void operator=(const NoCopy&);
	NoCopy(const NoCopy&);
protected:
	NoCopy(){}
};

template <class... Args> void IgnoreIfUnused(Args&&...){}

template < typename T, size_t N >
constexpr int countof( T ( & /*arr*/ )[ N ] )
{
	static_assert(std::extent< T[ N ] >::value <= size_t(std::numeric_limits<int>::max()), "Static array size exceeds int. countof will not work");
	return int(std::extent< T[ N ] >::value);
}

#define DECLARE_FORWARDING_CONSTRUCTOR(Self, forward_to)\
		template <class... TA> Self(TA&&... args): forward_to(std::forward<TA>(args)...) {}

#define DECLARE_FORWARDING_METHOD(return_type, method, forward_to_expr)\
		template <class... TA> return_type method(TA&&... args){ return (forward_to_expr).method(std::forward<TA>(args)...); }

#define DECLARE_FUNCTOR_TEMPLATE_T(ResultType, Name, ...)\
static struct Name##StructHelper{\
template <class T> ResultType operator()(__VA_ARGS__)const;\
} const Name;\
template <class T> ResultType Name##StructHelper::operator()(__VA_ARGS__)const


#define DECLARE_FUNCTOR_TEMPLATE_T1_T2(ResultType, Name, ...)\
static struct Name##StructHelper{\
template <class T1, class T2> ResultType operator()(__VA_ARGS__)const;\
} const Name;\
template <class T1, class T2> ResultType Name##StructHelper::operator()(__VA_ARGS__)const

#define DECLARE_PARAM_FUNCTOR_TEMPLATE_T(ResultType, Name, ...)\
template<class Param> struct Name##StructHelper{\
Param param_;\
explicit Name##StructHelper(Param&&  param):param_(std::forward<Param>(param)){}\
template <class T> ResultType operator()(__VA_ARGS__)const;\
};\
template <class Param> Name##StructHelper<Param> Name(Param&&  param){return Name##StructHelper<Param>(std::forward<Param>(param));}\
template <class Param> template <class T> ResultType Name##StructHelper<Param>::operator()(__VA_ARGS__)const

#define DECLARE_PARAM_FUNCTOR_TEMPLATE_T1_T2(ResultType, Name, ...)\
template<class Param> struct Name##StructHelper{\
Param param_;\
explicit Name##StructHelper(Param&&  param):param_(std::forward<Param>(param)){}\
template <class T1, class T2> ResultType operator()(__VA_ARGS__)const;\
};\
template <class Param> Name##StructHelper<Param> Name(Param&&  param){return Name##StructHelper<Param>(std::forward<Param>(param));}\
template <class Param> template <class T1, class T2> ResultType Name##StructHelper<Param>::operator()(__VA_ARGS__)const

DECLARE_FUNCTOR_TEMPLATE_T(std::size_t, CallStdHash, const T& v){
	return std::hash<T>()(v);
}

DECLARE_FUNCTOR_TEMPLATE_T1_T2(bool, EqualTo, T1&& v1, T2&& v2){
	return v1 == v2;
}

DECLARE_FUNCTOR_TEMPLATE_T1_T2(bool, NotEqualTo, T1&& v1, T2&& v2){
	return !(v1 == v2);
}

DECLARE_FUNCTOR_TEMPLATE_T1_T2(bool, ReferenceEqualTo, T1&& v1, T2&& v2){
	return std::addressof(v1) == std::addressof(v2);
}

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

template <class T>
typename std::remove_reference<T>::type CopyValue(T value)
{
	return value;
}

template <class Container, class Separator>
struct  ContainerWriterStruct
{
	const Container& data;
	const Separator separator;
	friend std::ostream& operator <<(std::ostream &output, const ContainerWriterStruct &x) {
		for(auto i = x.data.begin(); i != x.data.end(); ++i)
		{
			if (i != x.data.begin())
			{
				output << x.separator;
			}
			output << *i;
		}
		return output;
	}
};

template <class Container, class Separator>
struct  ContainerWriterStructWithSize 
{
	ContainerWriterStruct<Container, Separator> writer;
	const Separator container_name;
	const Separator container_begin;
	const Separator container_end;
	
	friend std::ostream& operator <<(std::ostream &output, const ContainerWriterStructWithSize  &x) {
		return output << x.container_name << "[" << x.writer.data.size() << "]" << x.container_begin << x.writer << x.container_end;
	}	
};

template <class Container, class Separator>
ContainerWriterStruct<Container, Separator> OutputContainer(const Container& data, Separator separator)
{
	return {data, separator};
}

template <class Container>
ContainerWriterStructWithSize<Container, const char*> OutputContainerWithSize(const Container& data, const char* name)
{
	return {{data, ", "}, name, "{", "}"};
}

template <class T>
struct optionalof
{
	optionalof(): has_value_(false), data_()
	{}

	~optionalof()
	{
		DestructIfNeeded();
	}

	optionalof(optionalof&&) = delete;

	optionalof(const optionalof& other)
		:has_value_(other.has_value_), data_()
	{
		if (has_value_)
		{
			new (TypedPtr()) T(other.get());
		}
	}

	optionalof& operator=(const optionalof& other) = delete;

	void operator=(const T& other)
	{
		DestructIfNeeded();
		new (TypedPtr()) T(other);
		has_value_ = true;
	}

	T& get()
	{
		assert(has_value());
		return *TypedPtr();
	}
	
	const T& get() const
	{
		assert(has_value());
		return *TypedPtr();
	}

	bool has_value()const
	{
		return has_value_;
	}
  private:
	T* TypedPtr()
	{
		return reinterpret_cast<T*>(data_);
	}
	const T* TypedPtr()const
	{
		return reinterpret_cast<const T*>(data_);
	}
	void DestructIfNeeded()
	{
		if (has_value_)
		{
			get().~T();
			has_value_ = false;
		}
	}
	char data_[sizeof(T)];
	bool has_value_;
};

template <class T>
struct Enumerator
{
	typedef T value_type;
	struct Impl
	{
		virtual void GoToBegin() = 0;
		virtual T GetAndMove() = 0;
		virtual bool AtEnd() = 0;
		virtual size_t size() = 0;
		virtual ~Impl(){}
	};
	
	template <class Iterator>
	struct NonSizedRangeImpl:Impl
	{
		NonSizedRangeImpl(Iterator a_begin, Iterator a_end)
			:begin(a_begin), current(a_begin), end(a_end)
		{}

		Iterator begin, current, end;
		
		void GoToBegin() override
		{
			current = begin;
		}

		T GetAndMove() override
		{
			T result = *current;
			++current;
			return result;
		}
		bool AtEnd() override
		{
			return !(current != end);
		}

		size_t size() override
		{
			throw std::logic_error("size is unknown for this iterator");
		}		
	};
	
	template <class Iterator>
	struct RangeImpl:Impl
	{
		RangeImpl(Iterator a_begin, Iterator a_end, size_t a_size)
			:begin(a_begin), current(a_begin), end(a_end), orig_size(a_size)
		{}
		
		Iterator begin, current, end;
		size_t orig_size;
		
		void GoToBegin() override
		{
			current = begin;
		}

		T GetAndMove() override
		{
			T result = *current;
			++current;
			return result;
		}
		bool AtEnd() override
		{
			return !(current != end);
		}

		size_t size() override
		{
			return orig_size;
		}		
	};
	

	template <class ConvertFrom>
	struct ConverterImpl:Impl
	{
		ConverterImpl(Enumerator<ConvertFrom> a_orig_enumerator)
			:orig_enumerator(a_orig_enumerator)
		{}
		Enumerator<ConvertFrom> orig_enumerator;
		std::function<T(const ConvertFrom&)> converter;
		
		void GoToBegin() override
		{
			orig_enumerator.GoToBegin();
		}

		T GetAndMove() override
		{
			return converter(orig_enumerator.GetAndMove());
		}
		
		bool AtEnd() override
		{
			return orig_enumerator.AtEnd();
		}
		
		size_t size() override
		{
			return orig_enumerator.size();
		}
	};
	
	template <class Iterator>
	static Enumerator<T> NonSizedRange(Iterator begin, Iterator end)
	{
		auto impl = std::make_shared<NonSizedRangeImpl<Iterator>>(begin, end);
		return Enumerator(impl);
	}

	template <class Iterator>
	static Enumerator<T> Range(Iterator begin, Iterator end, size_t size)
	{
		auto impl = std::make_shared<RangeImpl<Iterator>>(begin, end, size);
		return Enumerator(impl);
	}

	template <class ConvertFrom>
	static Enumerator<T> Converter(Enumerator<ConvertFrom> orig_enumerator,	std::function<T(const ConvertFrom&)> converter)
	{
		auto impl = std::make_shared<ConverterImpl<ConvertFrom>>(orig_enumerator);
		impl->converter = converter;
		return Enumerator(impl);
	}
	
	struct WrapperIterator
	{
		typedef typename std::remove_reference<T>::type TNoReference;
		typedef optionalof<T> TStorage;
		WrapperIterator()
			:enumerator_(nullptr), last_value_()
		{}
		explicit WrapperIterator(Enumerator& enumerator)
			:enumerator_(&enumerator), last_value_()
		{
			++(*this);
		}
		T operator*()const
		{
			return last_value_.get();
		}
		TNoReference* operator->()const
		{
			return &last_value_.get();
		}
		bool operator!=(const WrapperIterator& other)const
		{
			return enumerator_ != other.enumerator_;
		}
		void operator++()
		{
			assert(enumerator_);
			if(enumerator_->AtEnd())
			{
				enumerator_ = nullptr;
			}
			else
			{
				last_value_ = enumerator_->GetAndMove();
			}
		}
	private:
		Enumerator* enumerator_;
		TStorage last_value_;
	};
	
	WrapperIterator begin()
	{
		return WrapperIterator(*this);
	}

	WrapperIterator end()
	{
		return WrapperIterator();
	}
	
	void GoToBegin()
	{
		impl_->GoToBegin();
	}

	T GetAndMove()
	{
		return impl_->GetAndMove();
	}
	
	bool AtEnd()
	{
		return impl_->AtEnd();
	}
  
	size_t size()
	{
		return impl_->size();
	}

	Enumerator(){}
	
	void operator=(const Enumerator& other)
	{
		impl_ = other.impl_;
	}
  private:
	Enumerator(const std::shared_ptr<Impl>& impl)
		:impl_(impl)
	{}
	
	std::shared_ptr<Impl> impl_;
};

template <class TRange>
Enumerator<decltype(*std::begin(std::declval<TRange>()))> FullNonSizedRangeEnumerator(const TRange& range)
{
	return Enumerator<decltype(*std::begin(std::declval<TRange>()))>::NonSizedRange(range.begin(), range.end());
}

template <class TRange>
Enumerator<decltype(*std::begin(std::declval<TRange>()))> FullRangeEnumerator(const TRange& range)
{
	return Enumerator<decltype(*std::begin(std::declval<TRange>()))>::Range(range.begin(), range.end(), range.size());
}

template <class Func, class ConvertFrom>
static auto ConverterEnumerator(Enumerator<ConvertFrom> orig_enumerator, Func converter) -> Enumerator<decltype(converter(orig_enumerator.GetAndMove()))>
{
	return Enumerator<decltype(converter(orig_enumerator.GetAndMove()))>::template Converter<ConvertFrom>(orig_enumerator, converter);
}

template <class Func>
struct ConverterOfInnerEnumeratorImpl
{
	template <class ValueConvertedByFunc>
	auto operator()(Enumerator<ValueConvertedByFunc> orig_enumerator) const ->decltype(ConverterEnumerator(orig_enumerator, std::declval<Func>()))
	{
		return ConverterEnumerator(orig_enumerator, func_);
	}
	Func func_;
};

template <class Func>
static ConverterOfInnerEnumeratorImpl<Func> ConverterOfInnerEnumerator(Func converter)
{
	return ConverterOfInnerEnumeratorImpl<Func>{converter};
}

template <class Func, Func converter, class ConvertFrom>
static auto ConverterEnumeratorCFunc(Enumerator<ConvertFrom> orig_enumerator) -> decltype(ConverterEnumerator(orig_enumerator, converter))
{
	return ConverterEnumerator(orig_enumerator, converter);
}

#define STATIC_WITHTYPE_AS_TEMPLATE_PARAM(f) decltype(&(f)), (&(f))

template <class Iterator>
int slow_distance(Iterator from, Iterator to)
{
	int result = 0;
	while(from != to)
	{
		++from;
		++result;
	}
	return result;
}

template <class T>
struct auto_unique_ptr: std::unique_ptr<T>
{
	auto_unique_ptr():
		std::unique_ptr<T>(new T())
	{}
	
	auto_unique_ptr& operator=(auto_unique_ptr&& other) = default;
	auto_unique_ptr(auto_unique_ptr&& other):
		std::unique_ptr<T>(std::move(other))
	{}
};

//incapsulates reference to a deleter function (typically destructor) in  pointer so that the object can be destroyed where only a forward declaration is available
template <class T>
struct unique_deleter_ptr: std::unique_ptr<T, void(*)(T*)>
{
	typedef typename std::remove_cv<T>::type NonConstT;
	typedef std::unique_ptr<T, void(*)(T*)> BasePtr;
	unique_deleter_ptr():
		BasePtr(nullptr, nullptr)
	{}
	explicit unique_deleter_ptr(T* ptr):
		BasePtr(ptr, &unique_deleter_ptr::Delete)
	{}
	unique_deleter_ptr(T* ptr, void(deleter)(T*)):
		BasePtr(ptr, deleter)
	{}
	unique_deleter_ptr(unique_deleter_ptr<NonConstT>&& other):
		BasePtr(other.release(), reinterpret_cast<typename BasePtr::deleter_type>(other.get_deleter()))
	{
	}
	//fix dunplicate specialization with enable_if
	unique_deleter_ptr(unique_deleter_ptr&& other):
		BasePtr(other)
	{
	}
	unique_deleter_ptr& operator=(unique_deleter_ptr&& other) = default;

	static void Delete(T* ptr)
	{
		std::default_delete<T> cheking_incomplete_type_deleter;
		cheking_incomplete_type_deleter(ptr);
	}
};

template <class T>
unique_deleter_ptr<T> as_deleter_ptr(T* ptr)
{
	return unique_deleter_ptr<T>(ptr, unique_deleter_ptr<T>::Delete);
}

template <class T>
std::unique_ptr<T> as_unique_ptr(T* ptr)
{
	return std::unique_ptr<T>(ptr);
}

template <class T>
struct ImplicitlyConvertible 
{
	explicit ImplicitlyConvertible(T value):
		value_(std::forward<T>(value))
	{}
	template <class T2>
	operator T2()
	{
		return T2(std::forward<T>(value_));
	}
private:
  T value_;
};

template <class T>
ImplicitlyConvertible<typename std::remove_reference<T>::type&&> MoveToResultType(T&& value)
{
	return ImplicitlyConvertible<typename std::remove_reference<T>::type&&>(std::move(value));
}

template <class T>
struct auto_pimpl_ptr: unique_deleter_ptr<T>
{
	auto_pimpl_ptr():
		unique_deleter_ptr<T>(new T{})
	{}
};

//pimpl declaration that requires pimpl_ setting
#define DECLARE_PIMPL  struct Impl; unique_deleter_ptr<Impl> impl_
//pimpl declaration that allows constructing of on objetc only in places where impl is defined and always value-initializes impl
#define DECLARE_AUTO_PIMPL  struct Impl; auto_pimpl_ptr<Impl> impl_

template <class Enumerator>
struct EnumeratorTraits
{
	typedef  decltype(*(std::declval<Enumerator>().begin())) ValueType;
};

template <class T>
T Initialized()
{
	//value_initiaalized
	return T{};
}

template <class T, class MemberType, class... Args>
T Initialized(MemberType T::*member, MemberType value, Args... args)
{
	auto result = Initialized<T>(args...);
	result.*member = value;
	return result;
}

template <class Container, class... Args>
auto emplaced_back(Container& container, Args... args) -> decltype(std::declval<Container>().back())
{
	container.emplace_back(std::forward<Args>(args)...);
	return container.back();
}

 inline std::ostream& PrintTo(std::ostream& stream)
{
	return stream;
}

template <class T, class... Args>
std::ostream& PrintTo(std::ostream& stream, T&& arg0, Args&&... args)
{
	return PrintTo(stream << std::forward<T>(arg0), std::forward<Args>(args)...);
}

template <class... Args>
void L(Args&&... args)
{
	PrintTo(std::clog, std::forward<Args>(args)...) << std::endl;
}

class RandomGenerator
{
  public:
	int operator()()
	{
		return dis(gen);
	}
  private:
	std::mt19937 gen{std::random_device()()};
	std::uniform_int_distribution<> dis;
};

template <class NarrowInteger, class WideInteger> NarrowInteger narrow_cast(const WideInteger& i)
{
	static_assert(
		(std::is_unsigned<NarrowInteger>::value && std::is_unsigned<WideInteger>::value) ||
		(std::is_signed<NarrowInteger>::value && std::is_signed<WideInteger>::value)
		, "Integers must be both unsigned or signed");
	NarrowInteger result = i;
	assert(WideInteger(result) == i);
	return result;
}

template <class  WideInteger, class NarrowInteger> WideInteger wide_cast(const NarrowInteger& i)
{
	static_assert(
		(std::is_unsigned<NarrowInteger>::value && std::is_unsigned<WideInteger>::value) ||
		(std::is_signed<NarrowInteger>::value && std::is_signed<WideInteger>::value)
		, "Integers must be both unsigned or signed");
	static_assert(sizeof(NarrowInteger) <= sizeof(WideInteger), "wide_cast can't cast to smaller size"); //can add support for casting to wider signed type
	return WideInteger(i);
}

template <class ValueContainer, decltype(ValueContainer::value) expected_value>
struct StaticAsserter{
	static bool const value = ValueContainer::value == expected_value;
	static_assert(value, "StaticAsserter reports that it's type arguments is bad, see below for other assertion explaining why");
};

//need a macro to be able to cast to non-publicly visible bases
#define TO_BASE_CAST(TBase, child) \
	([](decltype(child)&& child_in_lambda)->TBase{ \
		static_assert(StaticAsserter<std::is_base_of<typename std::remove_reference<TBase>::type, typename std::remove_reference<decltype(child)>::type>, true>::value, "TO_BASE_CAST can cast only to base");\
		return static_cast<TBase>(child_in_lambda); \
	}(child))

template <class Signed> typename std::make_unsigned<Signed>::type unsigned_cast(const Signed& i)
{
	typedef typename std::make_unsigned<Signed>::type Unsigned;
	assert(i >= Signed(0));
	return Unsigned(i);
}

template <class Unsigned> typename std::make_signed<Unsigned>::type signed_cast(const Unsigned& i)
{
	typedef typename std::make_unsigned<Unsigned>::type Signed;
	assert(i <= Unsigned(std::numeric_limits<Signed>::max()));
	auto result = Signed(i);
	assert(result >=0);
	return result;
}

#define BIT_SIZEOF(type_or_value) (sizeof(type_or_value) * CHAR_BIT)

template <class Unsigned> Unsigned RotateBitsRight(const Unsigned& v, int shift)
{
	static_assert(std::is_unsigned<Unsigned>::value, "Can be used only with Unsigned type");
	assert(shift >= 0);
	static const int kBitsInValue = BIT_SIZEOF(Unsigned);
	int normalized_shift = shift % kBitsInValue; //the division is a quite heavy operation in a bit rotation function implementation, but kBitsInValue in all typical saces would be power of 2.
	return (v >> normalized_shift) ^ (v << (kBitsInValue - normalized_shift));
}
//std::hash for pairs
namespace std
{
	template <class T1, class T2>
	struct hash<std::pair<T1, T2>>
	{
		std::size_t operator()(const std::pair<T1, T2>& value)const
		{
			static const int kHalfSizeShift = BIT_SIZEOF(std::size_t)/2; //bit count in std::size_t
			return CallStdHash(value.first) ^ RotateBitsRight(CallStdHash(value.second), kHalfSizeShift);
		}
	};
}

//hash function for small collections
DECLARE_FUNCTOR_TEMPLATE_T(std::size_t, SmallCollectionHash, const T& collection){
	std::size_t result = 0;
	if (collection.size() != 0)
	{
		const int rotation_per_value = std::max(1, signed_cast(BIT_SIZEOF(result) / collection.size()));
		for(auto i : collection)
		{
			result = RotateBitsRight(result, rotation_per_value) ^ CallStdHash(i);
		}
	}
	return result;
}
