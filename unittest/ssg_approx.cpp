#include <gtest/gtest.h>
#include "ssg_approx.h"
#include "mock/ring.h"

typedef Mock::Ring R;

TEST(ApproxSigantureGroebner, ZeroPolys)
{
	std::unique_ptr<const R> ring = Mock::CreateRing();
	std::unique_ptr<R> out_ring = Mock::CreateRing();
	const typename R::PolysSet empty_in;
	typename R::PolysSet out;
	ApproxSigantureGroebner<R>(*ring, empty_in, *out_ring, out);
	out.expect_empty();
}