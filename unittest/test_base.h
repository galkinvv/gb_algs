#pragma once
#include <gtest/gtest.h>

#ifdef PRIVATE_TEST
namespace UnitTest{
	struct PRIVATE_TEST;
}
using UnitTest::PRIVATE_TEST;
#define TEST_METHOD(method) \
TEST_F(PRIVATE_TEST, method) {this->method();}
#endif


template <class SubComparator>
struct ContainerEqualExpect
{
	explicit ContainerEqualExpect(const SubComparator& equal)
		:equal_(equal)
	{}
	
	template <class Container1, class Container2>
	void ExpectEqual(const Container1& container1, const Container2& container2)const
	{
		auto container2_cur = container2.begin();
		for(auto container1_cur:container1)
		{
			ASSERT_NE(container2_cur, container2.end());
			auto elements_comparator = [this](decltype(container1_cur) v1, decltype(*container2_cur) v2){return SavingResultCompareWithComparator(v1,v2, equal_);};
			EXPECT_PRED2(elements_comparator, container1_cur, *container2_cur);
			++container2_cur;
		}
		auto iterator_comparator = [this](decltype(container2_cur) v1, decltype(container2.end()) v2){return SavingResultCompare(v1,v2);};
		EXPECT_PRED2(iterator_comparator, container2_cur, container2.end());
		no_asserts_happen_ = true;
	}

	template <class Container1, class Container2>
	bool operator()(const Container1& container1, const Container2& container2)const
	{
		all_equal_ = true;
		no_asserts_happen_ = false;
		ExpectEqual(container1, container2);
		return all_equal_ && no_asserts_happen_;
	}
private:
	template <class V1, class V2>
	bool SavingResultCompare(V1& v1, V2& v2)const
	{
		return SavingResultCompareWithComparator(v1, v2);
	}

	template <class V1, class V2, class Comparator =  std::equal_to<V1>>
	bool SavingResultCompareWithComparator(V1& v1, V2& v2, const Comparator& comp = std::equal_to<V1>())const
	{
		bool result = comp(v1, v2);
		if (!result)
		{
			all_equal_ = false;
		}
		return result;
	}
	
	const SubComparator&  equal_;
	mutable bool all_equal_ = false;
	mutable bool no_asserts_happen_ = false;
};

template <class SubComparator> ContainerEqualExpect<SubComparator> ExpecterContainerEqual(const SubComparator& equal)
{
	return ContainerEqualExpect<SubComparator>(equal);
}
