#include "ringfast_z2_simpledegrevlex.h"
#include "simple_mon.h"
#include "simple_sig.h"
#include <cassert>
#include <list>
#include <deque>

typedef std::uint64_t SigIdx;

struct Signature {
	SigIdx index;
	SimpleMon mon;
};

struct FastPolynomial : std::vector<SimpleMon> {};

struct RingFastZ2SimpleDegrevlexBase::LPolyImpl {
	FastPolynomial value;
	std::vector<FastPolynomial> reconstruction_info;
	Signature sig;
};


struct MultLPoly
{
	MultLPoly(const RingFastZ2SimpleDegrevlexBase::LPoly& a_poly):
		poly(a_poly)
	{}
	const RingFastZ2SimpleDegrevlexBase::LPoly& poly;
	SimpleMon mul_by;
};

struct RingFastZ2SimpleDegrevlexBase::MultLPolysQueueImpl : std::list<MultLPoly>
{};

struct RingFastZ2SimpleDegrevlexBase::LPolysResultImpl : std::deque<RingFastZ2SimpleDegrevlexBase::LPoly>
{};

struct RingFastZ2SimpleDegrevlexBase::Impl
{
	std::set<SigIdx> all_indices;
	std::list<LPoly> initial_polys;
};

RingFastZ2SimpleDegrevlexBase::RingFastZ2SimpleDegrevlexBase(){}

struct LessComparer
{
	typedef bool (&CompareMethod)(const SimpleMon& m1, const SimpleMon& m2);
	const CompareMethod compare_method;

	bool operator()(const SimpleMon& m1, const SimpleMon& m2)const
	{
		return compare_method(m1, m2);
	}

	bool operator()(const MultLPoly& p1, const MultLPoly& p2)const
	{
		return SigLess(p1, p2, *this);
	}

};
void RingFastZ2SimpleDegrevlexBase::AddLabeledPolyBeforeImpl(int new_var_index, int new_poly_index_in_rec_basis, Enumerator<CrossRingInfo::PerVariableData> monomial, LPolysResult& reducers, const LPoly& poly_before)
{
	//a new polynomial that would allow reducing monomial would be added
	LPoly& new_poly = emplaced_back(*reducers);
	SimpleMon old_mons;
	SimpleMon new_mon;
	for (auto var: monomial) {
		old_mons[var.index] = var.degree;
	}
	new_mon[ new_var_index] = 1;
	//set new_poly->value to (old_mon - new_mon)
	new_poly->value.push_back(old_mons);
	new_poly->value.push_back(new_mon);
	//reconstruction info - exactly corresponds to added polynomial. So use single "1" monomial with "1" coefficient
	new_poly->reconstruction_info.resize(new_poly_index_in_rec_basis + 1);
	new_poly->reconstruction_info[new_poly_index_in_rec_basis].emplace_back();
	//set signature to new value, just before poly_before
	//monomial is the same, since poly_before is already multiplied polynomial
	new_poly->sig.mon = poly_before->sig.mon;
	const auto greater_idx_pos = impl_->all_indices.find(new_poly->sig.index);
	assert(greater_idx_pos != impl_->all_indices.end());
	SigIdx smaller_idx =
		(greater_idx_pos == impl_->all_indices.begin()) ?
			std::numeric_limits<SigIdx>::min() :
			*std::prev(greater_idx_pos);
	new_poly->sig.index  = smaller_idx + (*greater_idx_pos - smaller_idx)/2;
	//those assertions can fail during correct operation of an algorithm if too much additions were performed in same place
	assert(new_poly->sig.index > smaller_idx);
	assert(new_poly->sig.index < *greater_idx_pos);
	impl_->all_indices.insert(new_poly->sig.index);
}

RingFastZ2SimpleDegrevlexBase::LPoly RingFastZ2SimpleDegrevlexBase::DequeueSigSmallest(MultLPolysQueue& queue)
{
	LessComparer ms_less{MDegRevLexless<SimpleMon>};
	auto& queue_impl = *(queue);
	assert(!queue_impl.empty());
	auto min_it = std::min_element(queue_impl.begin(), queue_impl.end(), ms_less);
	assert(min_it != queue_impl.end());
	LPoly result;
	result->sig.index = min_it->poly->sig.index;
	result->sig.mon = Mmul(min_it->poly->sig.mon, min_it->mul_by);
	result->value = Pmul(min_it->poly->value, min_it->mul_by);
	result->reconstruction_info.reserve(min_it->poly->reconstruction_info.size());
	for (auto rec_info_poly:min_it->poly->reconstruction_info) {
		result->reconstruction_info.push_back(Pmul(rec_info_poly, min_it->mul_by));
	}
	queue_impl.erase(min_it);
	return result;
}

void RingFastZ2SimpleDegrevlexBase::ExtendQueueBySpairPartsAndFilterUnneeded(const LPolysResult& left_parts, const LPoly& right_part, MultLPolysQueue& queue)
{
	LessComparer ms_less{MDegRevLexless<SimpleMon>};
	assert(!IsZero(right_part));
	for (const auto& left_part:*left_parts) {
		if (IsZero(left_part)) {
			continue;
		}
		MultLPoly left(left_part);
		left.mul_by = ToLCMMultiplier(HM(left_part->value), HM(right_part->value));
		MultLPoly right(right_part);
		right.mul_by = ToLCMMultiplier(HM(right_part->value), HM(left_part->value));
		MultLPoly &new_lpoly = ms_less(left, right) ? right : left;
		bool was_supeseded;
		auto& queue_impl = *queue;
		for(auto existing_lpoly : queue_impl) {
			if (IsSupersededBy(new_lpoly, existing_lpoly, ms_less)) {
				was_supeseded = true;
				break;
			}
		}
		if (was_supeseded) {
			continue;
		}
		queue_impl.remove_if(std::bind(IsSupersededBy<MultLPoly, LessComparer>, std::placeholders::_1, new_lpoly, ms_less));
		queue_impl.emplace_back(new_lpoly);
	}
}

CrossRingInfo::PerVariableData VarDataFromMapItem(const std::pair<const int, int>& index_degree_pair)
{
	return CrossRingInfo::PerVariableData::FromDI(index_degree_pair.second, index_degree_pair.first);
}


Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> RingFastZ2SimpleDegrevlexBase::FieldAgnosticReconstructionInfoPolysImpl(const LPoly& poly)
{
	auto poly_enumerator = FullRangeEnumerator(poly->reconstruction_info);
	typedef decltype(poly_enumerator.GetAndMove()) PolynomialAsRange;
	auto poly_with_mon_enumerator = ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullRangeEnumerator<PolynomialAsRange>)>(poly_enumerator);
	typedef decltype(poly_with_mon_enumerator.GetAndMove().GetAndMove()) MonomialAsRange;

	auto poly_with_mon_with_varpairs_enumerator =
		ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM((
					ConverterEnumeratorCFunc<
					STATIC_WITHTYPE_AS_TEMPLATE_PARAM(FullNonSizedRangeEnumerator<MonomialAsRange>),
					MonomialAsRange>
				))>(poly_with_mon_enumerator);


	typedef decltype(poly_with_mon_with_varpairs_enumerator.GetAndMove().GetAndMove()) MonomialAsEnum;
	typedef decltype(poly_with_mon_with_varpairs_enumerator.GetAndMove().GetAndMove().GetAndMove()) VarAsMapItem;

	return ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM((
				ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM((
						ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM(VarDataFromMapItem), VarAsMapItem>
						)),
				MonomialAsEnum>
			))>(poly_with_mon_with_varpairs_enumerator);
}

Enumerator<CrossRingInfo::PerVariableData> RingFastZ2SimpleDegrevlexBase::FieldAgnosticReconstructionInfoTopImpl(const LPoly& poly)
{
	return ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM(VarDataFromMapItem)>(FullRangeEnumerator(HM(poly->value)));
}

RingFastZ2SimpleDegrevlexBase::LPolysResult RingFastZ2SimpleDegrevlexBase::FillWithTrivialSyzygiesOfNonMultElements(const MultLPolysQueue& queue)
{
	LPolysResult result;
	LessComparer ms_less{MDegRevLexless<SimpleMon>};
	const auto& queue_impl = *queue;
	for(auto i0 = queue_impl.begin(); i0 != queue_impl.end(); ++i0) {
		assert(MDeg(i0->mul_by)  == 0 && MDeg(i0->poly->sig.mon) == 0) ;
		for(auto i1 = std::next(i0); i1 != queue_impl.end(); ++i1) {
			LPoly syz_part[2];
			//value and reconstruction info are zero
			syz_part[0]->sig.index = i0->poly->sig.index;
			syz_part[0]->sig.mon = HM(i1->poly->value);
			syz_part[1]->sig.index = i1->poly->sig.index;
			syz_part[1]->sig.mon = HM(i0->poly->value);
			//mul_by is identity nomomial
			MultLPoly mult[2] = {{syz_part[0]}, {syz_part[1]}};
			int greater_sig_idx = ms_less(mult[0], mult[1]) ? 1 : 0;
			InsertInResult(std::move(syz_part[greater_sig_idx]), result);
		}
	}
	return result;
}

void RingFastZ2SimpleDegrevlexBase::InsertInResult(LPoly&& poly, LPolysResult& result)
{
	result->emplace_back(std::move(poly));
}

bool RingFastZ2SimpleDegrevlexBase::IsZero(const LPoly& poly)
{
	return poly->value.empty();
}

void RingFastZ2SimpleDegrevlexBase::Normalize(LPoly& poly)
{
	IgnoreIfUnused(poly);
}

RingFastZ2SimpleDegrevlexBase::MultLPolysQueue RingFastZ2SimpleDegrevlexBase::PutInQueueExtendLabeledPolysImpl(Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input)
{
	const SigIdx polys_count = input.size();
	static const SigIdx kMaxIdx = std::numeric_limits<SigIdx>::max();
	const SigIdx kStepPerIdx = kMaxIdx/polys_count;
	int cur_poly_index = 0;
	MultLPolysQueue result;
	for (auto poly:input)
	{
		++cur_poly_index;
		LPoly& labeled_poly = emplaced_back(impl_->initial_polys);
		for (auto vars_colection : poly){
			SimpleMon mon;
			for(auto var: vars_colection)
			{
				mon[var.index] = var.degree;
			}
			labeled_poly->value.push_back(mon);
		}
		labeled_poly->reconstruction_info.resize(input.size());
		labeled_poly->reconstruction_info[cur_poly_index].resize(1);
		labeled_poly->sig.index = static_cast<SigIdx>(cur_poly_index)*kStepPerIdx;
		result->emplace_back(MultLPoly{labeled_poly});
		impl_->all_indices.insert(labeled_poly->sig.index);
	}
	return result;
}

bool RingFastZ2SimpleDegrevlexBase::QueueEmpty(const MultLPolysQueue& queue)
{
	return queue->empty();
}

void RingFastZ2SimpleDegrevlexBase::ReduceCheckingSignatures(LPoly& poly, LPolysResult& reducers)
{
	LessComparer ms_less{MDegRevLexless<SimpleMon>};
	bool failed_to_find_reducer = false;
	while(!poly->value.empty() && !failed_to_find_reducer) {
		failed_to_find_reducer = true;
		for(const auto& reducer:*reducers) {
			if(!reducer->value.empty()) {
				auto divider = DivideIfCan(HM(poly->value), HM(reducer->value));
				if (!divider) {
					continue;
				}
				MultLPoly mult_reductor{reducer};
				mult_reductor.mul_by = *divider;
				MultLPoly to_reduce{poly};

				if (!ms_less(mult_reductor, to_reduce)) {
					continue;
				}
				failed_to_find_reducer = false;
				poly->value = PReduceZ2(poly->value,  reducer->value, *divider);
				assert(poly->reconstruction_info.size() == reducer->reconstruction_info.size());
				for (int rec_info_idx = 0; rec_info_idx < int(poly->reconstruction_info.size()); ++rec_info_idx) {
					auto& rec_info_ref = poly->reconstruction_info[rec_info_idx];
					rec_info_ref = PAddZ2(rec_info_ref, reducer->reconstruction_info[rec_info_idx], *divider);
				}
				break;
			}
		}
	}
}
