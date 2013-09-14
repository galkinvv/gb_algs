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
