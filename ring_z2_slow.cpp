#include "ring_z2_slow.h"
#include "iopolynomset.h"
#include <algorithm>
#include <cassert>
namespace
{
	template <class TMonomial>
	TMonomial Mmul(const TMonomial& m1, const TMonomial& m2)
	{
		TMonomial result = m1;
		for(auto i = m2.begin(); i!=m2.end(); ++i)
		{
			result[i->first]+=i->second;
		}
		return result;
	}

	template <class TMonomial, class TPolynomial>
	TPolynomial Pmul(const TPolynomial& p, const TMonomial& m)
	{
		TPolynomial result;
		for(auto i = p.begin(); i!=p.end(); ++i)
		{
			result.push_back(Mmul(*i, m));
		}
		return result;
	}

	template <class TMonomial>
	int MDeg(const TMonomial& m)
	{
		int result = 0;
		for(auto i = m.begin(); i!=m.end; ++i)
		{
			result+=i->second;
		}
		return result;
	}

	template <class T, class Comp = std::less<T>>
	bool unequal(const T& t1, const T& t2, bool& less, Comp compare = Comp())
	{
		if (compare(t1, t2))
		{
			less = true;
			return true;
		}
		if (compare(t2, t1))
		{
			less = false;
			return true;
		}
		return false;
	}

	template <class TMonomial>
	bool MDegRevLexless(const TMonomial& m1, const TMonomial& m2)
	{
		bool item_less;
		if (unequal(MDeg(m1), MDeg(m2), item_less))
		{
			return item_less; //x^2 < y^3
		}

		if (m1 == m2) return false; //x = x
		auto i1 = m1.rbegin();
		auto i2 = m2.rbegin();
		for(;;++i1, ++i2)
		{
			assert (i1 != m1.rend() && i2 != m2.rend()); //equal case already checked
			if (unequal(i1->first, i2->first, item_less))
			{
				return !item_less; // z < x
			}
			if (unequal(i1->second, i2->second, item_less))
			{
				return !item_less; // xz^2 < x^2z
			}
		}
	}

	template <class TMultLPoly>
	auto MultSig(const TMultLPoly& mp) -> decltype(mp.mul_by)
	{
		return Mmul(mp.mul_by, mp.poly.sig_mon);
	}
	template <class TMultLPoly>
	bool SigLess(const TMultLPoly& mp1, const TMultLPoly& mp2)
	{
		bool item_less;
		if (unequal(mp1.poly.sig_index, mp2.poly.sig_index, item_less))
		{
			return item_less; //index1 < index2
		}
		if (unequal(MultSig(mp1), MultSig(mp2), item_less))
		{
			return item_less; //index1 < index2
		}
		return false; //equal
	}

	typedef RingZ2Slow::FastAssociatedLabeledRingWithTracking FR;

	struct RingZ2SlowIoData: IOData<RingZ2Slow>
	{
		using IOData<RingZ2Slow>::in_;
		using IOData<RingZ2Slow>::in_ring_;
	};
}

FR::LPoly FR::DequeueSigSmallest(MultLPolysQueue& queue)
{
	auto min = std::min_element(queue.begin(), queue.end(), SigLess<MultLPoly>);
	assert(min != queue.end());
	LPoly result;
	result.sig_index = min->poly.sig_index;
	result.sig_mon = Mmul(min->poly.sig_mon, min->mul_by);
	result.value = Pmul(min->poly.value, min->mul_by);
	result.reconstruction_info = Pmul(min->poly.reconstruction_info, min->mul_by);
	queue.erase(min);
	return result;
}

void FR::PutInQueueExtendLabeledPolys(const PolysSet& in, MultLPolysQueue& queue)
{
	double current_sig_index = 1;
	for(auto poly = in.begin(); poly != in. end(); ++poly, current_sig_index += 1)
	{
		MultLPoly mp;		
		mp.poly.value = *poly;
		mp.poly.sig_index = current_sig_index;
		queue.push_back(mp);
	}
}

std::unique_ptr<IOData<RingZ2Slow>> RingZ2Slow::Create(const F4MPI::IOPolynomSet& in)
{
	auto* data = new RingZ2SlowIoData();
	RingZ2Slow ring;
	ring.CopyTo(data->in_ring_); 
	PolysSet in_polys;
	//TODO convert in to in_polys
	data->in_ = in_polys;
	return std::unique_ptr<IOData<RingZ2Slow>>(data);
}

F4MPI::IOPolynomSet RingZ2Slow::ConvertResult(std::unique_ptr<IOData<RingZ2Slow>>& result)
{
	result.reset();
	//TODO convert result 
	F4MPI::IOPolynomSet converted_result;
	return converted_result;
}