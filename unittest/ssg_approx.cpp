#define PRIVATE_TEST ApproxSignatureGroebnerUT
#include "test_base.h"
#include "ssg_approx.h"
#include "mock/ring.h"

struct ApproxSignatureGroebnerUT: ::testing::Test{
	void ZeroPolys()
	{
		typedef Mock::Ring R;
		auto io_data = R::Create("");
		ApproxSignatureGroebner<R>::Do(*io_data);
		std::string result = R::ConvertResult(io_data);
		EXPECT_EQ(result, "");
	}
};
TEST_METHOD(ZeroPolys)
