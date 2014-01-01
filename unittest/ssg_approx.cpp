#define PRIVATE_TEST ApproxSignatureGroebnerUT
#include "test_base.h"
#include "ssg_approx.h"
#include "finite_field.h"
#include "z_ring.h"
#include "mock/ring.h"

struct PRIVATE_TEST: ::testing::Test{
	typedef Mock::Ring<CrossRingInfo::MonomialMetadata<CrossRingInfo::MonomialOrder::DegRevLex>, FiniteField<ZPlusRing8>> R;
	typedef R::Field Field;
	typedef R::MonomialMetadata MonomialMetadata;

	static Field CreateField()
	{
		return Field::CreateZpFieldWithChar(2);
	}
	void ZeroPolys()
	{
		auto field = CreateField();
		MonomialMetadata monomial_metadata;
		monomial_metadata.var_count = 1;
		R out_ring{monomial_metadata, field};
		R::IOData::IOPolynomSet io_poly_set_in {monomial_metadata, field};
		R::IOData io_data {io_poly_set_in, out_ring};
		ApproxSignatureGroebner<R, Mock::FastRingWithTracking>::Do(io_data);
		EXPECT_EQ(io_data.out_data->begin(), io_data.out_data->end());
		auto mapping = io_data.out_ring.MonMapping();
		int out_var_count = io_data.out_data->Metadata().var_count;
		for(int i=0; i<out_var_count; ++i)
		{
			bool has_non_zeroes = false;
			for (auto old_var:(*mapping)[i])
			{
				has_non_zeroes |= old_var.degree > 0;
			}
			EXPECT_TRUE(has_non_zeroes);
		}
	}
};
TEST_METHOD(ZeroPolys)
