#define PRIVATE_TEST ApproxSignatureGroebnerUT
#include "test_base.h"
#include "ssg_approx.h"
#include "z_field.h"
#include "mock/ring.h"

struct PRIVATE_TEST: ::testing::Test{
	void ZeroPolys()
	{
		typedef Mock::Ring<CrossRingInfo::MonomialMetadata<CrossRingInfo::MonomialOrder::DegRevLex>, ZField8> R;
		typedef R::Field Field;
		typedef R::MonomialMetadata MonomialMetadata;
		Field field;
		MonomialMetadata monomial_metadata;
		monomial_metadata.var_count = 1;
		R out_ring{monomial_metadata, field};
		R::IOData::IOPolynomSet io_poly_set_in {monomial_metadata, field};
		R::IOData io_data {io_poly_set_in, out_ring};
		ApproxSignatureGroebner<R, Mock::FastRingWithTracking>::Do(io_data);
		EXPECT_EQ(io_data.out_data->begin(), io_data.out_data->end());
	}
};
TEST_METHOD(ZeroPolys)
