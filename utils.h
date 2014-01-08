#pragma once
#include <limits>
#include <type_traits>
#include <utility>
#include <initializer_list>
#include <ostream>
#include <functional>
#include <memory>
#include <cassert>

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
	const Separator& separator;
};

template <class Container, class Separator>
ContainerWriterStruct<Container, Separator> OutputContainer(const Container& data, const Separator& separator)
{
	return {data, separator};
}

template <class Container, class Separator>
std::ostream& operator<<(std::ostream& s, const ContainerWriterStruct<Container, Separator>& writer)
{
	for(auto i = writer.data.begin(); i != writer.data.end(); ++i)
	{
		if (i != writer.data.begin())
		{
			s << writer.separator;
		}
		s << *i;
	}
	return s;
}

template <class T>
std::ostream& operator<<(std::ostream& s, const std::initializer_list<T>& ilist)
{
	return s << "{" << OutputContainer(ilist, ", ") << "}";
}

template <class T>
struct optionalof
{
	optionalof(): has_value_(false)
	{}

	~optionalof()
	{
		DestructIfNeeded();
	}

	optionalof(optionalof&&) = delete;

	optionalof(const optionalof& other)
		:has_value_(other.has_value_)
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
		virtual T GetAndMove() = 0;
		virtual bool AtEnd() = 0;
		virtual ~Impl(){}
	};
	
	template <class Iterator>
	struct RangeImpl:Impl
	{
		RangeImpl(Iterator a_current, Iterator a_end)
			:current(a_current), end(a_end)
		{}
		Iterator current, end;
		
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
	};
	
	template <class ConvertFrom>
	struct ConverterImpl:Impl
	{
		ConverterImpl(Enumerator<ConvertFrom> a_orig_enumerator)
			:orig_enumerator(a_orig_enumerator)
		{}
		Enumerator<ConvertFrom> orig_enumerator;
		std::function<T(const ConvertFrom&)> converter;
		
		T GetAndMove() override
		{
			return converter(orig_enumerator.GetAndMove());
		}
		
		bool AtEnd() override
		{
			return orig_enumerator.AtEnd();
		}
	};
	
	template <class Iterator>
	static Enumerator<T> Range(Iterator begin, Iterator end)
	{
		auto impl = std::make_shared<RangeImpl<Iterator>>(begin, end);
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
	
	T GetAndMove()
	{
		return impl_->GetAndMove();
	}
	
	bool AtEnd()
	{
		return impl_->AtEnd();
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
Enumerator<decltype(*std::begin(std::declval<TRange>()))> FullRangeEnumerator(const TRange& range)
{
	return Enumerator<decltype(*std::begin(std::declval<TRange>()))>::Range(range.begin(), range.end());
}

template <class Func, class ConvertFrom>
static auto ConverterEnumerator(Enumerator<ConvertFrom> orig_enumerator, Func converter) -> Enumerator<decltype(converter(orig_enumerator.GetAndMove()))>
{
	return Enumerator<decltype(converter(orig_enumerator.GetAndMove()))>::template Converter<ConvertFrom>(orig_enumerator, converter);
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
	static void Delete(T* ptr)
	{
		std::default_delete<T> cheking_incomplete_type_deleter;
		cheking_incomplete_type_deleter(ptr);
	}
};

template <class T>
unique_deleter_ptr<T> create_deleter_ptr(T* ptr)
{
	return unique_deleter_ptr<T>(ptr, unique_deleter_ptr<T>::Delete);
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

#define DECLARE_PIMPL  struct Impl; unique_deleter_ptr<Impl> impl_