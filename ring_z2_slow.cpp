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
		for(auto i = m.begin(); i!=m.end(); ++i)
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

	template <class TPolynomial, class TMonomial = decltype(*TPolynomial().cbegin())>
	TMonomial HM(const TPolynomial& p) 
	{
		assert(!p.empty());
		return *std::max_element(p.begin(), p.end(), MDegRevLexless<TMonomial>);
	}

	template <class TMultLPoly, class TMonomial = decltype(TMultLPoly().mul_by)>
	TMonomial  MultSig(const TMultLPoly& mp)
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
			return item_less; //sigmon1 < sigmon2
		}
		return false; //equal
	}

	typedef RingZ2Slow::FastAssociatedLabeledRingWithTracking FR;

	struct RingZ2SlowIoData: IOData<RingZ2Slow>
	{
		using IOData<RingZ2Slow>::in_;
		using IOData<RingZ2Slow>::in_ring_;
	};
	static const char kFirstVar = 'a';
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

void FR::FillWithTrivialSyzygiesOfNonMultElements(const MultLPolysQueue& queue, LPolysResult& to_fill)
{
	for(auto i0 = queue.begin(); i0 != queue.end(); ++i0)
	{
		assert(MDeg(i0->mul_by)  == 0 && MDeg(i0->poly.sig_mon) == 0) ;
		for(auto i1 = std::next(i0); i1 != queue.end(); ++i1)
		{
			MultLPoly syz_part[2];
			syz_part[0].poly.sig_index = i0->poly.sig_index;
			syz_part[0].poly.sig_mon = HM(i1->poly.value);
			syz_part[1].poly.sig_index = i1->poly.sig_index;
			syz_part[1].poly.sig_mon = HM(i0->poly.value);;
			int greater_sig_idx = static_cast<int>(SigLess(syz_part[0], syz_part[1]));
			to_fill.push_back(syz_part[greater_sig_idx].poly);
		}
	}
}

struct RingZ2Slow::Impl
{
	int var_count;
};

RingZ2Slow::RingZ2Slow()
	:impl_(new Impl())
{}

RingZ2Slow::~RingZ2Slow()
{}

void RingZ2Slow::CopyTo(RingZ2Slow& other)const
{
	*other.impl_ = *impl_;
}

std::unique_ptr<IOData<RingZ2Slow>> RingZ2Slow::Create(const F4MPI::IOPolynomSet& in)
{
	auto* data = new RingZ2SlowIoData();
	if (
		in.field_char != 2 ||
		in.mon_order != F4MPI::CMonomial::degrevlexOrder ||
		in.type != F4MPI::FieldType::Z
		)
	{
		throw std::runtime_error("unsupported parameters for RingZ2Slow");
	}
	data->in_ring_.impl_->var_count = in.var_count;
	PolysSet in_polys;
	for(auto poly:in.polys)
	{
		Polynomial poly_in_ring;
		for(auto mon = poly.m_begin(); mon != poly.m_end(); ++mon)
		{
			Monomial mon_in_ring;
			for (int var_idx = 0; var_idx < data->in_ring_.impl_->var_count; ++var_idx)
			{
				if (int deg = mon->getDegree(var_idx))
				{
					mon_in_ring[kFirstVar] = deg;
				}
			}
			poly_in_ring.push_back(mon_in_ring);
		}
		in_polys.push_back(poly_in_ring);
	}
	data->in_ = in_polys;
	return std::unique_ptr<IOData<RingZ2Slow>>(data);
}

F4MPI::IOPolynomSet RingZ2Slow::ConvertResult(std::unique_ptr<IOData<RingZ2Slow>>& result)
{
	F4MPI::IOPolynomSet converted_result;
	converted_result.field_char = 2;
	converted_result.mon_order = F4MPI::CMonomial::degrevlexOrder;
	converted_result.type = F4MPI::FieldType::Z;
	converted_result.var_count = result->out_ring.impl_->var_count;
	for(auto poly_in_ring:result->out)
	{
		F4MPI::CPolynomial poly;
		for(auto mon_in_ring: poly_in_ring)
		{
			std::vector<F4MPI::CMonomialBase::Deg> mon(converted_result.var_count);
			for (auto var_deg_pair:mon_in_ring)
			{
				int var_idx = var_deg_pair.first - kFirstVar;
				mon[var_idx] =var_deg_pair.second;
			}
			poly.addTerm(F4MPI::CModular(1), F4MPI::CMonomial(mon));
		}
		converted_result.polys.push_back(poly);
	}
	result.reset();
	return converted_result;
}