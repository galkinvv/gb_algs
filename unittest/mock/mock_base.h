#pragma once
#ifdef PRIVATE_TEST
#define FRIEND_FOR_TEST friend class ::UnitTest::PRIVATE_TEST;
#else
#define FRIEND_FOR_TEST
#endif