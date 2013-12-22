#include "ring_z2_slow.h"
#include <algorithm>
#include <cassert>
#include <unordered_set>
namespace
{

struct MonomialHash {
	template <class TMonomial>
	std::size_t operator()(TMonomial const& mon) const {
		std::size_t result = 0;
		for(auto var : mon) {
			result = result *(1 + (1 << 4)) ^ pair_hash(var);
		}
		return result;
	}

	template <class T>
	std::size_t std_hash(const T& var) const {
		return std::hash<T>()(var);
	}
	template <class Pair>
	std::size_t pair_hash(const Pair& pair) const {
		return std_hash(pair.first) + (1 + (1 << 2)) * std_hash(pair.second) ;
	}
};

template <class TMonomial>
TMonomial Mmul(const TMonomial& m1, const TMonomial& m2)
{
	TMonomial result = m1;
	for(auto i = m2.begin(); i!=m2.end(); ++i) {
		result[i->first]+=i->second;
	}
	return result;
}

template <class TMonomial>
TMonomial ToLCMMultiplier(const TMonomial& to_mul, const TMonomial& lcm_with)
{
	TMonomial result;
	for(auto var_deg_pair:lcm_with) {
		auto same_var_pos= to_mul.find(var_deg_pair.first);
		auto var_deg_result = var_deg_pair;

		if (same_var_pos != to_mul.end()) {
			var_deg_result.second -= same_var_pos->second;
		}
		if (var_deg_result.second >0) {
			result.insert(var_deg_result);
		}
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

template <class TPolynomial>
bool IsZeroImpl(const TPolynomial& p)
{
	return p.empty();
}

template <class TPolynomial, class TConstMonomialRef = decltype(*TPolynomial().cbegin())>
TConstMonomialRef HM(const TPolynomial& p)
{
	assert(!p.empty());
	return *std::max_element(p.begin(), p.end(), MDegRevLexless<TConstMonomialRef>);
}

template <class TPolynomial, class TMonomial>
TPolynomial PAdd(const TPolynomial& not_muled_item, const TPolynomial& muled_item, const TMonomial& mul_by)
{
	std::unordered_set<TMonomial, MonomialHash> presence_count;
	for (auto mon:not_muled_item) {
		presence_count.insert(mon);
	}
	for (auto mon_by_notmuled:muled_item) {
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
	return PAdd(poly_to_red, by, mul_by);
}

template <class TMultLPoly, class TMonomial = decltype(TMultLPoly().mul_by)>
TMonomial  MultSig(const TMultLPoly& mp)
{
	return Mmul(mp.mul_by, mp.poly.sig_mon);
}

template <class TMultLPoly, class TMonomial = decltype(TMultLPoly().mul_by)>
TMonomial  MultHM(const TMultLPoly& mp)
{
	return Mmul(mp.mul_by, HM(mp.poly.value));
}

template <class TMultLPoly>
bool SigLess(const TMultLPoly& mp1, const TMultLPoly& mp2)
{
	bool item_less;
	if (unequal(mp1.poly.sig_index, mp2.poly.sig_index, item_less)) {
		return item_less; //index1 < index2
	}
	auto sig1 = MultSig(mp1);
	auto sig2 = MultSig(mp2);
	if (unequal(sig1, sig2, item_less, MDegRevLexless<decltype(sig1)>)) {
		return item_less; //sigmon1 < sigmon2
	}
	return false; //equal
}

template <class TMultLPoly>
bool IsSupersededBy(const TMultLPoly& maybe_supded, const TMultLPoly& sup_by)
{
	if (maybe_supded.poly.sig_index !=  sup_by.poly.sig_index) {
		return false; //can't supersede wth different indexes
	}
	auto sig_supded = MultSig(maybe_supded);
	auto sig_by = MultSig(sup_by);
	if (!DivideIfCan(sig_supded, sig_by)) {
		return false; //index_old not divisible by index_new
	}
	if (IsZeroImpl(sup_by.poly.value)) return true; //zero polynomial with dividing sig
	if (IsZeroImpl(maybe_supded.poly.value)) return false; //zero polynomial with sig divisible by non-zero

	auto sig_supded_hm_by = Mmul(sig_supded, MultHM(sup_by));
	auto sig_by_hm_supded = Mmul(sig_by, MultHM(maybe_supded));
	if (MDegRevLexless(sig_by_hm_supded, sig_supded_hm_by)) {
		return false; //HM/S for maybe_supded is smaller
	}
	return true; // HM/S for maybe_supded is not smaller
}

static const char kFirstVar = 'a';
}

struct RingZ2SlowBase::Impl
{
};

struct RingZ2SlowBase::InPolysSetWithOrigMetadata
{
};

struct RingZ2SlowBase::OutPolysSetForVariyingMetadata
{
};

RingZ2SlowBase::RingZ2SlowBase(int var_count):
	impl_(create_deleter_ptr(new Impl))
{}

bool RingZ2SlowBase::ConstructAndInsertNormalizedImpl(const unique_deleter_ptr<const InPolysSetWithOrigMetadata>& prepared_input, 
		Enumerator<CrossRingInfo::PerVariableData> top_info,  
		Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input_polys_mons, 
		const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& result)
{
	//TODO: solve via matrix and fast field
	for (auto top_var: top_info)
	{
		
	}
	for (auto poly: input_polys_mons)
	{
		for (auto mon: poly)
		{
			for (auto var: mon)
			{
				
			}
		}				
	}
	return true;
}

void RingZ2SlowBase::ConvertResultToFixedMetadataImpl(const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& constructed_result, CrossRingInfo::MonomialListListWithCoef<ImplementedOrder, ImplementedField>& basic_result)
{
	//TODO
}

void RingZ2SlowBase::ExtendWithMonomialImpl(Enumerator<CrossRingInfo::PerVariableData> info)
{
	//TODO
}

unique_deleter_ptr<RingZ2SlowBase::OutPolysSetForVariyingMetadata> RingZ2SlowBase::PrepareEmptyResult()
{
	return create_deleter_ptr(new OutPolysSetForVariyingMetadata());
}

unique_deleter_ptr<const RingZ2SlowBase::InPolysSetWithOrigMetadata> RingZ2SlowBase::PrepareForReconstructionImpl(Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input)
{
	//TODO
	auto result = create_deleter_ptr(new InPolysSetWithOrigMetadata());
	return result;
}
/*
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
				for (int rec_info_idx = 0; rec_info_idx < int(poly.reconstruction_info.size()); ++rec_info_idx) {
					auto& rec_info_ref = poly.reconstruction_info[rec_info_idx];
					rec_info_ref = PAdd(rec_info_ref, reducer.reconstruction_info[rec_info_idx], *divider);
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
	return IsZeroImpl(poly.value);
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
	assert(!IsZeroImpl(right_part.value));
	for (auto left_part:left_parts) {
		if (IsZeroImpl(left_part.value)) {
			continue;
		}
		MultLPoly left;
		left.poly = left_part;
		left.mul_by = ToLCMMultiplier(HM(left_part.value), HM(right_part.value));
		MultLPoly right;
		right.poly = right_part;
		right.mul_by = ToLCMMultiplier(HM(right_part.value), HM(left_part.value));
		MultLPoly &new_lpoly = SigLess(left, right) ? right : left;
		bool was_supeseded;
		for(auto existing_lpoly : queue) {
			if (IsSupersededBy(new_lpoly, existing_lpoly)) {
				was_supeseded = true;
				break;
			}
		}
		if (was_supeseded) {
			continue;
		}
		queue.erase(
		    std::remove_if(
		        queue.begin(), queue.end(), std::bind(IsSupersededBy<MultLPoly>, std::placeholders::_1, new_lpoly)
		    ),
		    queue.end()
		);
		queue.push_back(new_lpoly);
	}
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

bool RingZ2Slow::ConstructAndInsertNormalized(const PolysSet& in, const ReconstructionInfo& info, PolysSet& out)
{
	Polynomial result;
	assert(info.size() == in.size());
	auto in_poly = in.cbegin();
	for(auto rec_info : info)
	{
		assert(in_poly != in.cend());
		for(auto mon : rec_info)
		{
			result = PAdd(result, *in_poly, mon);
		}
		++in_poly;
	}
	assert(!IsZeroImpl(result));
	assert(HM(result) == info.top);
	out.push_back(result);
	return true;
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
*/
