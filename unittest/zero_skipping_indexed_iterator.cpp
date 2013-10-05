#include <gtest/gtest.h>
#include <vector>
#include <array>
#include <deque>
#include "zero_skipping_indexed_iterator.h"

template <class Container>
struct ZeroSkippingIndexedIteratorTest :  ::testing::Test
{};

struct Item
{
	Item(int a_value, int a_index)
		:value(a_value), index(a_index)
	{}
	int value, index;
};

TYPED_TEST_CASE_P(ZeroSkippingIndexedIteratorTest);

TYPED_TEST_P(ZeroSkippingIndexedIteratorTest, SkipMidAndLastZeroes) {
  typedef ZeroSkippingIndexedIterator<TypeParam, Item> It;
  TypeParam data;
  data.push_back(1);
  data.push_back(0);
  data.push_back(1);
  data.push_back(3);
  data.push_back(0);
  auto current = It::BeginOf(data.begin(), data. end());
  auto end = It::EndOf(data.begin(), data. end());
  std::array<Item, 3> expected{Item(1,0), Item(1,2), Item(3,3)};
  auto ecurrent = std::begin(expected);
  for(;current != end; ++current, ++ecurrent)
  {
	  EXPECT_EQ((*current).value, ecurrent->value);
	  EXPECT_EQ((*current).index, ecurrent->index);
	  ASSERT_NE(ecurrent, std::end(expected));
  }
  EXPECT_EQ(ecurrent, std::end(expected));
}

TYPED_TEST_P(ZeroSkippingIndexedIteratorTest, SkipFirstZeroes) {
  typedef ZeroSkippingIndexedIterator<TypeParam, Item> It;
  TypeParam data;
  data.push_back(0);
  data.push_back(0);
  data.push_back(1);
  data.push_back(3);
  data.push_back(5);
  auto current = It::BeginOf(data.begin(), data. end());
  auto end = It::EndOf(data.begin(), data. end());
  std::array<Item, 3> expected{Item(1,2), Item(3,3), Item(5,4)};
  auto ecurrent = std::begin(expected);
  for(;current != end; ++current, ++ecurrent)
  {
	  EXPECT_EQ((*current).value, ecurrent->value);
	  EXPECT_EQ((*current).index, ecurrent->index);
	  ASSERT_NE(ecurrent, std::end(expected));
  }
  EXPECT_EQ(ecurrent, std::end(expected));
}

TYPED_TEST_P(ZeroSkippingIndexedIteratorTest, Empty) {
  typedef ZeroSkippingIndexedIterator<TypeParam, Item> It;
  TypeParam data;
  auto current = It::BeginOf(data.begin(), data. end());
  auto end = It::EndOf(data.begin(), data. end());
  EXPECT_FALSE(current != end);
}

TYPED_TEST_P(ZeroSkippingIndexedIteratorTest, SingleNonZero) {
  typedef ZeroSkippingIndexedIterator<TypeParam, Item> It;
  TypeParam data;
  data.push_back(7);
  auto current = It::BeginOf(data.begin(), data. end());
  auto end = It::EndOf(data.begin(), data. end());
  EXPECT_TRUE(current != end);  
  EXPECT_EQ((*current).value, 7);
  EXPECT_EQ((*current).index, 0);
  ++current;
  EXPECT_FALSE(current != end);
}

TYPED_TEST_P(ZeroSkippingIndexedIteratorTest, AllZeroes) {
  typedef ZeroSkippingIndexedIterator<TypeParam, Item> It;
  TypeParam data;
  data.push_back(0);
  data.push_back(0);
  data.push_back(0);
  data.push_back(0);
  data.push_back(0);
  auto current = It::BeginOf(data.begin(), data. end());
  auto end = It::EndOf(data.begin(), data. end());
  EXPECT_FALSE(current != end);
}

REGISTER_TYPED_TEST_CASE_P(ZeroSkippingIndexedIteratorTest, SkipMidAndLastZeroes, AllZeroes, SingleNonZero, Empty, SkipFirstZeroes);

INSTANTIATE_TYPED_TEST_CASE_P(ContVector, ZeroSkippingIndexedIteratorTest, std::vector<int>);
INSTANTIATE_TYPED_TEST_CASE_P(ContDeque, ZeroSkippingIndexedIteratorTest, std::deque<int>);