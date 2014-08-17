#pragma once
#include <gtest/gtest.h>
#include "utils.h"

#ifdef PRIVATE_TEST
namespace UnitTest{
	struct PRIVATE_TEST;
}
using UnitTest::PRIVATE_TEST;
#define TEST_METHOD(method) \
TEST_F(PRIVATE_TEST, method) {this->method();}
#endif

#ifdef PRIVATE_TYPED_TEST
namespace UnitTest{
	template <class Param>
	struct PRIVATE_TYPED_TEST;
}
using UnitTest::PRIVATE_TYPED_TEST;
#define TEST_METHOD(method) \
TYPED_TEST_P(PRIVATE_TYPED_TEST, method) {PRIVATE_TYPED_TEST<TypeParam>::method();}
TYPED_TEST_CASE_P(PRIVATE_TYPED_TEST);
#endif

#define EXPECT_ASSERT(expr) EXPECT_DEATH((expr), "")

//EXPECT_FUNCTION_NO_FAILURES is local variable name used in functions EXPECT_FUNCTION_*
#define EXPECT_FUNCTION_BEGIN bool EXPECT_FUNCTION_NO_FAILURES = true;
#define EXPECT_FUNCTION_RETURN return EXPECT_FUNCTION_NO_FAILURES;

template <class Comparator>
struct SavingValueExpector
{
	explicit SavingValueExpector(const Comparator& comp, bool& no_failures):
		no_failures_(no_failures),
		comp_(comp)
	{}
	template <class... Args>
	bool operator()(Args... args)const
	{
		bool result = comp_(std::forward<Args>(args)...);
		if (!result)
		{
			no_failures_ = false;
		}
		return result;
	}
	bool& no_failures_;
	const Comparator& comp_;
};

template <class Comparator> SavingValueExpector<Comparator> SaveValue(const Comparator& equal, bool& no_failures)
{
	return SavingValueExpector<Comparator>(equal, no_failures);
}

#define EXPECT_1(pred, arg0) EXPECT_PRED1(SaveValue(pred, EXPECT_FUNCTION_NO_FAILURES), arg0)

#define EXPECT_2(pred, arg0, arg1) EXPECT_PRED2(SaveValue(pred, EXPECT_FUNCTION_NO_FAILURES), arg0, arg1)

#define EXPECT_3(pred, arg0, arg1, arg2) EXPECT_PRED3(SaveValue(pred, EXPECT_FUNCTION_NO_FAILURES), arg0, arg1, arg2)

#define ASSERT_2(pred, arg0, arg1) do\
{bool no_failures_as_assert = true; EXPECT_PRED2(SaveValue(pred, no_failures_as_assert), arg0, arg1);\
if (!no_failures_as_assert){return false;}\
}while(false)

DECLARE_PARAM_FUNCTOR_TEMPLATE_T1_T2(bool, ContainerEqualExpect, const T1& container1, const T2& container2)
{
		EXPECT_FUNCTION_BEGIN
		auto container2_cur = container2.begin();
		for(auto container1_cur:container1)
		{
			ASSERT_2(NotEqualTo, container2_cur, container2.end());
			EXPECT_2(param_, container1_cur, *container2_cur);
			++container2_cur;
		}
		EXPECT_2(EqualTo, container2_cur, container2.end());
		EXPECT_FUNCTION_RETURN;
}

template <class T> 
std::ostream& operator <<(std::ostream &output, const std::vector<T> &x) {
	return output << OutputContainerWithSize(x, "vector");
}

template <class T>
std::ostream& operator<<(std::ostream& output, const std::initializer_list<T>& x)
{
	return output << OutputContainerWithSize(x, "initializer_list");
}

