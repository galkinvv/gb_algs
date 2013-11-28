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
		WrapperIterator()
			:enumerator_(nullptr), last_value_()
		{}
		explicit WrapperIterator(Enumerator& enumerator)
			:enumerator_(enumerator), last_value_()
		{
			++(*this);
		}
		T operator*()const
		{
			return last_value_;
		}
		T* operator->()const
		{
			return &last_value_;
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
		T last_value_;
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
  private:
	Enumerator(const std::shared_ptr<Impl>& impl)
		:impl_(impl)
	{}
	
	const std::shared_ptr<Impl> impl_;
};

template <class TRange>
Enumerator<typename std::add_lvalue_reference<typename std::add_const<decltype(*std::begin(std::declval<TRange>()))>::type>::type> FullRangeEnumerator(const TRange& range)
{
	return Enumerator<typename std::add_lvalue_reference<typename std::add_const<decltype(*std::begin(std::declval<TRange>()))>::type>::type>::Range(range.begin(), range.end());
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
