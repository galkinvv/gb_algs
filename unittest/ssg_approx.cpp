#define PRIVATE_TEST ApproxSigantureGroebnerUT
#include "test_base.h"
#include "ssg_approx.h"
#include "mock/ring.h"

struct ApproxSigantureGroebnerUT: ::testing::Test{
	void ZeroPolys()
	{
		typedef Mock::Ring R;
		std::unique_ptr<const R> ring = Mock::CreateRing();
		std::unique_ptr<R> out_ring = Mock::CreateRing();
		const typename R::PolysSet empty_in;
		typename R::PolysSet out;
		ApproxSigantureGroebner<R>::Do(*ring, empty_in, *out_ring, out);
		EXPECT_EQ(out.size(), 0);
	}
};
TEST_METHOD(ZeroPolys)
