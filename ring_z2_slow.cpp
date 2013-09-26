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
	for(auto i = m2.begin(); i!=m2.end(); ++i) {
		result[i->first]+=i->second;
	}
	return result;
}

template <class TMonomial, class TPolynomial>
TPolynomial Pmul(const TPolynomial& p, const TMonomial& m)
{
	TPolynomial result;
	for(auto i = p.begin(); i!=p.end(); ++i) {
		result.push_back(Mmul(*i, m));
	}
	return result;
}

template <class TMonomial>
int MDeg(const TMonomial& m)
{
	int result = 0;
	for(auto i = m.begin(); i!=m.end(); ++i) {
		result+=i->second;
	}
	return result;
}

template <class T, class Comp = std::less<T>>
bool unequal(const T& t1, const T& t2, bool& less, Comp compare = Comp())
{
	if (compare(t1, t2)) {
		less = true;
		return true;
	}
	if (compare(t2, t1)) {
		less = false;
		return true;
	}
	return false;
}

template <class TMonomial>
bool MDegRevLexless(const TMonomial& m1, const TMonomial& m2)
{
	bool item_less;
	if (unequal(MDeg(m1), MDeg(m2), item_less)) {
		return item_less; //x^2 < y^3
	}

	if (m1 == m2) return false; //x = x
	auto i1 = m1.rbegin();
	auto i2 = m2.rbegin();
	for(;; ++i1, ++i2) {
		assert (i1 != m1.rend() && i2 != m2.rend()); //equal case already checked
		if (unequal(i1->first, i2->first, item_less)) {
			return !item_less; // z < x
		}
		if (unequal(i1->second, i2->second, item_less)) {
			return !item_less; // xz^2 < x^2z
		}
	}
}

template <class TMonomial>
std::unique_ptr<TMonomial> DivideIfCan(const TMonomial& m1, const TMonomial& m2)
{
	std::unique_ptr<TMonomial> result(new TMonomial);
	for (auto var_deg:m1) {
		auto deg = var_deg.second;
		auto var = var_deg.first;
		assert(deg >0);
		auto var_deg_in_m2 = m2.find(var);
		if (var_deg_in_m2 == m2.end() || var_deg_in_m2 ->second < deg ) {
			return nullptr;
		}
		int deg_diff = deg - var_deg_in_m2 ->second;
		if (deg_diff == 0) {
			continue;
		}
		(*result)[var] = deg_diff;
	}
	return result;
}

template <class TPolynomial, class TConstMonomialRef = decltype(*TPolynomial().cbegin())>
TConstMonomialRef HM(const TPolynomial& p)
{
	assert(!p.empty());
	return *std::max_element(p.begin(), p.end(), MDegRevLexless<TConstMonomialRef>);
}

template <class TPolynomial, class TMonomial>
TPolynomial PSubtract(const TPolynomial& poly_to_red, const TPolynomial& by, const TMonomial& mul_by)
{
	std::set<TMonomial> presence_count;
	for (auto mon:poly_to_red) {
		presence_count.insert(mon);
	}
	for (auto mon_by_notmuled:by) {
		auto mon = Mmul(mon_by_notmuled, mul_by);
		auto pair_position_already_present_flag = presence_count.insert(mon);
		if (!pair_position_already_present_flag.second) {
			presence_count.erase(pair_position_already_present_flag.first);
		}
	}
	TPolynomial result;
	result.reserve(presence_count.size());
	result.insert(result.begin(), presence_count.begin(), presence_count.end());
	return result;
}

template <class TPolynomial, class TMonomial>
TPolynomial PReduce(const TPolynomial& poly_to_red, const TPolynomial& by, const TMonomial& mul_by)
{
	assert(HM(poly_to_red) == Mmul(HM(by), mul_by));
	return PSubtract(poly_to_red, by, mul_by);
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
	if (unequal(mp1.poly.sig_index, mp2.poly.sig_index, item_less)) {
		return item_less; //index1 < index2
	}
	if (unequal(MultSig(mp1), MultSig(mp2), item_less)) {
		return item_less; //sigmon1 < sigmon2
	}
	return false; //equal
}

template <class TMultLPoly>
bool IsSupersededBy(const TMultLPoly& maybe_supded, const TMultLPoly&  sup_by)
{
	if (sup_by.poly.value.empty()) return false;
	if (maybe_supded.poly.value.empty()) return false;
	if (maybe_supded.poly.sig_index !=  sup_by.poly.sig_index) {
		return false; //index_old != index_new
	}
	auto sig_old = MultSig(maybe_supded);
	auto sig_new = MultSig(sup_by);
	if (!DivideIfCan(sig_old, sig_new))
	{
		return false; //index_old not divisible by index_new
	}
	//auto sig_old_hm_new = 

	//TODO
	return false; //equal
}

typedef RingZ2Slow::FastAssociatedLabeledRingWithTracking FR;

struct RingZ2SlowIoData: IOData<RingZ2Slow> {
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
	result.reconstruction_info.reserve(min->poly.reconstruction_info.size());
	for (auto rec_info_poly:min->poly.reconstruction_info) {
		result.reconstruction_info.push_back(Pmul(rec_info_poly, min->mul_by));
	}
	queue.erase(min);
	return result;
}

void FR::PutInQueueExtendLabeledPolys(const PolysSet& in, MultLPolysQueue& queue)
{
	int current_poly_index = 0;
	for(auto poly: in) {
		MultLPoly mp;
		for (auto mon : poly) {
			mp.poly.value.push_back(mon);
		}
		mp.poly.reconstruction_info.resize(in.size());
		mp.poly.reconstruction_info[current_poly_index].resize(1);
		mp.poly.sig_index = double(current_poly_index) + 1;
		queue.push_back(mp);
		++current_poly_index;
	}
}

void FR::FillWithTrivialSyzygiesOfNonMultElements(const MultLPolysQueue& queue, LPolysResult& to_fill)
{
	for(auto i0 = queue.begin(); i0 != queue.end(); ++i0) {
		assert(MDeg(i0->mul_by)  == 0 && MDeg(i0->poly.sig_mon) == 0) ;
		for(auto i1 = std::next(i0); i1 != queue.end(); ++i1) {
			MultLPoly syz_part[2];
			//value and reconstruction info are zero
			syz_part[0].poly.sig_index = i0->poly.sig_index;
			syz_part[0].poly.sig_mon = HM(i1->poly.value);
			syz_part[1].poly.sig_index = i1->poly.sig_index;
			syz_part[1].poly.sig_mon = HM(i0->poly.value);;
			int greater_sig_idx = static_cast<int>(SigLess(syz_part[0], syz_part[1]));
			to_fill.push_back(syz_part[greater_sig_idx].poly);
		}
	}
}

void FR::ReduceCheckingSignatures(LPoly& poly, LPolysResult& reducers)
{
	bool failed_to_find_reducer = false;
	while(!poly.value.empty() && !failed_to_find_reducer) {
		failed_to_find_reducer = true;
		for(auto reducer:reducers) {
			if(!reducer.value.empty()) {
				auto divider = DivideIfCan(HM(poly.value), HM(reducer.value));
				if (!divider) {
					continue;
				}
				MultLPoly mult_reductor;
				mult_reductor.poly = reducer;
				mult_reductor.mul_by = *divider;
				MultLPoly to_reduce;
				to_reduce.poly = poly;
				if (!SigLess(mult_reductor, to_reduce)) {
					continue;
				}
				poly.value = PReduce(poly.value,  reducer.value, *divider);
				assert(poly.reconstruction_info.size() == reducer.reconstruction_info.size());
				for (int rec_info_idx = 0; rec_info_idx <poly.reconstruction_info.size(); ++rec_info_idx) {
					auto& rec_info_ref = poly.reconstruction_info[rec_info_idx];
					rec_info_ref = PSubtract(rec_info_ref, reducer.reconstruction_info[rec_info_idx], *divider);
				}
				failed_to_find_reducer = false;
				break;
			}
		}
	}
}

RingZ2Slow::ReconstructionInfo FR::FieldAgnosticReconstructionInfo(const LPoly& poly)
{
	ReconstructionInfo result;
	result.assign(poly.reconstruction_info.begin(),  poly.reconstruction_info.end());
	result.top = HM(poly.value);
	return result;
}

void FR::ExtendRingWithMonomialToHelpReconstruct(const LPoly& poly, LPolysResult& reducers)
{
	throw std::logic_error("ring extension requsted for Z2");
}

bool FR::IsZero(const LPoly& poly)
{
	return poly.value.empty();
}

void FR::Normalize(LPoly& poly)
{
	//Always normalized in Z2
}

void FR::InsertInResult(const LPoly& poly, LPolysResult& result)
{
	result.push_back(poly);
}


void FR::ExtendQueueBySpairPartsAndFilterUnneeded(const LPolysResult& left_parts, const LPoly& right_part, MultLPolysQueue& queue)
{
	//TODO
}

struct RingZ2Slow::Impl {
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
	) {
		throw std::runtime_error("unsupported parameters for RingZ2Slow");
	}
	data->in_ring_.impl_->var_count = in.var_count;
	PolysSet in_polys;
	for(auto poly:in.polys) {
		Polynomial poly_in_ring;
		for(auto mon = poly.m_begin(); mon != poly.m_end(); ++mon) {
			Monomial mon_in_ring;
			for (int var_idx = 0; var_idx < data->in_ring_.impl_->var_count; ++var_idx) {
				if (int deg = mon->getDegree(var_idx)) {
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
	for(auto poly_in_ring:result->out) {
		F4MPI::CPolynomial poly;
		for(auto mon_in_ring: poly_in_ring) {
			std::vector<F4MPI::CMonomialBase::Deg> mon(converted_result.var_count);
			for (auto var_deg_pair:mon_in_ring) {
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
