#include "ring_z2_slow.h"
#include "sparse_matrix.h"
#include <algorithm>
#include <deque>
#include <cassert>
#include <unordered_set>
#include <list>
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

template <class TMultLPoly, class TMonomial = decltype(std::declval<TMultLPoly>().mul_by)>
TMonomial  MultSig(const TMultLPoly& mp)
{
	return Mmul(mp.mul_by, mp.poly.impl_->sig.mon);
}

template <class TMultLPoly, class TMonomial = decltype(std::declval<TMultLPoly>().mul_by)>
TMonomial  MultHM(const TMultLPoly& mp)
{
	return Mmul(mp.mul_by, HM(mp.poly.impl_->value));
}

template <class TMultLPoly, class MonCompare>
bool SigLess(const TMultLPoly& mp1, const TMultLPoly& mp2, MonCompare monLess)
{
	bool item_less;
	if (unequal(mp1.poly.impl_->sig.index, mp2.poly.impl_->sig.index, item_less)) {
		return item_less; //index1 < index2
	}
	auto sig1 = MultSig(mp1);
	auto sig2 = MultSig(mp2);
	if (unequal(sig1, sig2, item_less, monLess)) {
		return item_less; //sigmon1 < sigmon2
	}
	return false; //equal
}

template <class TMultLPoly, class MonCompare>
bool IsSupersededBy(const TMultLPoly& maybe_supded, const TMultLPoly& sup_by, MonCompare monLess)
{
	if (maybe_supded.poly.impl_->sig.index !=  sup_by.poly.impl_->sig.index) {
		return false; //can't supersede wth different indexes
	}
	auto sig_supded = MultSig(maybe_supded);
	auto sig_by = MultSig(sup_by);
	if (!DivideIfCan(sig_supded, sig_by)) {
		return false; //index_old not divisible by index_new
	}
	if (IsZeroImpl(sup_by.poly.impl_->value)) return true; //zero polynomial with dividing sig
	if (IsZeroImpl(maybe_supded.poly.impl_->value)) return false; //zero polynomial with sig divisible by non-zero

	auto sig_supded_hm_by = Mmul(sig_supded, MultHM(sup_by));
	auto sig_by_hm_supded = Mmul(sig_by, MultHM(maybe_supded));
	if (monLess(sig_by_hm_supded, sig_supded_hm_by)) {
		return false; //HM/S for maybe_supded is smaller
	}
	return true; // HM/S for maybe_supded is not smaller
}

static const char kFirstVarOnDebugOutput = 'a';
}

bool MonomialLessDegRevLex(const BaseMon& m1, const BaseMon& m2)
{
	return MDegRevLexless(m1, m2);
}

struct SlowMon : BaseMon {};
struct SlowPolynomial : std::vector<SlowMon> {};

typedef FastZ2SlowBasedRingBase::FastMonomial FastMonomial;
typedef std::uint64_t SigIdx;

struct Signature {
	SigIdx index;
	FastMonomial mon;
};

//struct SLowPol : std::vector<SlowMon> {};

struct RingZ2SlowBase::Impl {
	//this struct shoulld contain data corresponding to polynomial ring, but not to some polynomials of this ring
	Impl(int keeped_vars_count)
		:keeped_vars_count_(keeped_vars_count)
	{}
	const int keeped_vars_count_;
	std::vector<SlowMon> new_variables;
};

//input polynomials in ring with original monomial count
struct RingZ2SlowBase::InPolysSetWithOrigMetadata : std::vector<SlowPolynomial> {};

//input polynomials in ring with new monomial count
struct RingZ2SlowBase::OutPolysSetForVariyingMetadata : std::vector<SlowPolynomial> {};

RingZ2SlowBase::RingZ2SlowBase(int var_count):
	impl_(new Impl(var_count))
{}

bool RingZ2SlowBase::ConstructAndInsertNormalizedImpl(const InPolysSetWithOrigMetadata& reconstruction_basis,
        Enumerator<CrossRingInfo::PerVariableData> top_info,
        Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input_polys_mons,
        const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& result)
{
	SlowMon top;
	for (auto top_var: top_info) {
		top[top_var.index] = top_var.degree;
	}
	auto mon_less = [this](const SlowMon& m1, const SlowMon& m2){return this->MonomialLess(m1, m2);};

	struct MultSlowPolyIterator
	{
		MultSlowPolyIterator(const SlowPolynomial& a_poly)
			:poly(a_poly)
		{}
		SlowMon mul_by;
		SlowPolynomial::const_iterator it;
		SlowPolynomial::const_iterator first_smaller_top_it;//iterator to first element smaller than top (or end)
		const SlowPolynomial& poly;
		SlowMon GetMultMon()
		{
			return Mmul(mul_by, *it);
		}
	};
	//prepare:
	//collection of multiplied input polynomials
	//this associates each polynomial in question (corresponding to monomialss of input_polys_mons) with unique number, corresponding to column number
	std::deque<MultSlowPolyIterator> mult_inputs_with_iters;
	//a sorted collection of their monomials in question (only gretear-or-equal than top_info) - each corresponds to matrix row
	std::set<SlowMon, decltype(mon_less)> high_mons(mon_less);
	auto mul_by_it = input_polys_mons.begin();
	auto orig_poly_it = reconstruction_basis.begin();
	for (;orig_poly_it != reconstruction_basis.end(); ++mul_by_it, ++orig_poly_it)
	{
		assert(mul_by_it != input_polys_mons.end());
		for (auto mon: *mul_by_it) {
			mult_inputs_with_iters.emplace_back(*orig_poly_it);
			auto &mult = mult_inputs_with_iters.back();
			SlowMon mul_by;
			for (auto var: mon) {
				mult.mul_by[var.index] = var.degree;
			}
			mult.it = mult.poly.begin();
			for(; mult.it != mult.poly.end(); ++mult.it )
			{
				auto orig_mon = mult.GetMultMon();
				if (mon_less(orig_mon, top))
				{
					break;
				}
				high_mons.insert(orig_mon);
			}
			mult.first_smaller_top_it = mult.it;
		}
	}
	//assert here that top correspond to smallest monomial in high_mons
	assert(!high_mons.empty() && high_mons.find(top) == high_mons.begin());
	ImplementedField field;
	std::vector<std::vector<SparseMatrix::Element<ImplementedField>>> matrix;
	matrix.reserve(high_mons.size());
	//populate matrix rows with (coef from ImplementedField; int column number)
	//move iterator to element smaller than top
	for(auto mult:mult_inputs_with_iters)
	{
		mult.it = mult.first_smaller_top_it;
	}
	for(auto mon:high_mons)
	{
		ImplementedField::Value one;
		field.SetOne(one);
		matrix.emplace_back();
		auto row = matrix.back();
		for(int col = 0; col < (int)mult_inputs_with_iters.size(); ++col)
		{
			auto mult = mult_inputs_with_iters[col];
			for(;mult.it != mult.poly.begin(); --mult.it)
			{
				auto mult_mon = mult.GetMultMon();
				if (!mon_less(mult_mon, mon))
				{
					if (mon_less(mon, mult_mon))
					{
						break;//too big
					}
					decltype(row)::value_type elem;
					elem.value = one;
					elem.column = col;
					row.push_back(elem);
					//eqaul - add to matrix
				}
				//too small find greater
			}
		}
	}
	std::vector<SparseMatrix::Element<ImplementedField>> solution;
	int max_diferent_numbers_in_coefficients = std::accumulate(reconstruction_basis.begin(), reconstruction_basis.end(), 0, [](int sum, const SlowPolynomial& poly){return sum + poly.size();});
	//send rows collection to solver that shoud assume that right-side column has 1 in first cell and last rows.size()-1 zeroes
	//solver returns only non-zeros - list of pairs (coef from ImplementedField; int column number)
	SparseMatrix::SolveWithRightSideContainigSingleOne(field, matrix, solution, max_diferent_numbers_in_coefficients);
	if(solution.empty())
	{
		//if solver fails - return false
		return false;
	}
	result->emplace_back();
	auto& new_poly = result->back();
	//add top_info with coef 1.
	new_poly.push_back(top);
	for (auto result_item:solution)
	{
		//calculate sum with coefs given from solver for monomials smaller than top_info, add it to result
		auto& to_mul = mult_inputs_with_iters[result_item.column];
		SlowPolynomial smaller_than_top_part;
		smaller_than_top_part.insert(smaller_than_top_part.begin(), to_mul.first_smaller_top_it, to_mul.poly.end());
		new_poly = PAdd(new_poly, smaller_than_top_part, to_mul.mul_by);
	}	
	return true;
}

void RingZ2SlowBase::ConvertResultToFixedMetadataImpl(const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& constructed_result, CrossRingInfo::MonomialListListWithCoef<ImplementedOrder, ImplementedField>& basic_result)
{
	ImplementedField::Value one;
	basic_result.Field().SetOne(one);
	for (auto poly: *constructed_result) {
		basic_result.BeginPolynomialConstruction(poly.size());
		for (auto mon: poly) {
			for (auto var: mon) {
				basic_result.AddVariable(CrossRingInfo::PerVariableData::FromDI(var.second, var.first));
			}
			basic_result.MonomialAdditionDone(one);
		}
	}
}

int RingZ2SlowBase::VarMappingImplReturningOldVarCount(std::vector<int>& new_monomial_vars) const
{
	
	new_monomial_vars.resize(impl_->keeped_vars_count_ * impl_->new_variables.size());
	int cur_new_var_start = 0;
	for(auto new_var:impl_->new_variables)
	{
		for (auto old_var:new_var)
		{
			new_monomial_vars[cur_new_var_start + old_var.first] = old_var.second;
		}
		cur_new_var_start += impl_->keeped_vars_count_;
	}
	assert(cur_new_var_start == int(new_monomial_vars.size()));
	return impl_->keeped_vars_count_;
}


RingZ2SlowBase::NewIndices RingZ2SlowBase::ExtendRingWithMonomialToHelpReconstructImpl(const unique_deleter_ptr<InPolysSetWithOrigMetadata>& reconstruction_basis, Enumerator<CrossRingInfo::PerVariableData> info)
{
	impl_->new_variables.emplace_back();
	SlowMon& old_mons = impl_->new_variables.back();
	for (auto var:info)
	{
		int& new_degree =old_mons[var.index];
		assert(0 == new_degree);
		new_degree = var.degree;
	}
	const int new_var_index = impl_->keeped_vars_count_ + impl_->new_variables.size() - 1;
	SlowMon new_mon;
	new_mon[new_var_index] = 1;
	reconstruction_basis->emplace_back();
	//new polynomial = old_mons - new_mon
	auto new_poly = std::prev(reconstruction_basis->end());	
	new_poly->push_back(old_mons);
	new_poly->push_back(new_mon);
	return Initialized<NewIndices>(&NewIndices::new_var_index, new_var_index, &NewIndices::new_poly_index, std::distance(reconstruction_basis->begin(), new_poly));
}

unique_deleter_ptr<RingZ2SlowBase::OutPolysSetForVariyingMetadata> RingZ2SlowBase::PrepareEmptyResult()
{
	return MoveToResultType(new OutPolysSetForVariyingMetadata());
}

unique_deleter_ptr<RingZ2SlowBase::InPolysSetWithOrigMetadata> RingZ2SlowBase::PrepareForReconstructionImpl(Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input)
{
	auto result = as_deleter_ptr(new InPolysSetWithOrigMetadata());
	for (auto poly: input) {
		result->emplace_back();
		auto& poly_to_populate = result->back();
		for (auto mon: poly) {
			poly_to_populate.emplace_back();
			auto& mon_to_populate = poly_to_populate.back();
			for (auto var: mon) {
				mon_to_populate[var.index] = var.degree;
			}
		}
	}
	return MoveToResultType(result);
}

struct FastPolynomial : std::vector<FastMonomial> {};

struct FastZ2SlowBasedRingBase::LPoly::Impl {	
	FastPolynomial value;
	std::vector<FastPolynomial> reconstruction_info;
	Signature sig;
};


struct MultLPoly
{
	MultLPoly(const FastZ2SlowBasedRingBase::LPoly& a_poly):
		poly(a_poly)
	{}
	const FastZ2SlowBasedRingBase::LPoly& poly;
	FastMonomial mul_by;
};


struct FastZ2SlowBasedRingBase::MultLPolysQueue::Impl : std::list<MultLPoly>
{};

struct FastZ2SlowBasedRingBase::LPolysResult::Impl : std::deque<FastZ2SlowBasedRingBase::LPoly>
{};

struct FastZ2SlowBasedRingBase::Impl
{
	std::set<SigIdx> all_indices;
};

FastZ2SlowBasedRingBase::FastZ2SlowBasedRingBase(){}

struct LessComparer
{
	typedef bool (FastZ2SlowBasedRingBase::*CompareMethod)(const FastMonomial& m1, const FastMonomial& m2) const;
	const CompareMethod compare_method;
	const FastZ2SlowBasedRingBase& ring_object;
	
	bool operator()(const FastMonomial& m1, const FastMonomial& m2)const
	{
		return (ring_object.*compare_method)(m1, m2);
	}

	bool operator()(const MultLPoly& p1, const MultLPoly& p2)const
	{
		return SigLess(p1, p2, *this);
	}

};
void FastZ2SlowBasedRingBase::AddLabeledPolyBeforeImpl(int new_var_index, int new_poly_index_in_rec_basis, Enumerator<CrossRingInfo::PerVariableData> monomial, LPolysResult& reducers, const LPoly& poly_before)
{
	//a new polynomial that would allow reducing monomial would be added
	reducers.impl_->emplace_back();
	LPoly& new_poly = reducers.impl_->back();
	new_poly.impl_.reset(new LPoly::Impl());
	FastMonomial old_mons;
	FastMonomial new_mon;
	for (auto var: monomial) {
		old_mons[var.index] = var.degree;
	}
	new_mon[ new_var_index] = 1;
	//set new_poly.impl_->value to (old_mon - new_mon)
	new_poly.impl_->value.push_back(old_mons);
	new_poly.impl_->value.push_back(new_mon);
	//reconstruction info - exactly corresponds to added polynomial. So use single "1" monomial with "1" coefficient
	new_poly.impl_->reconstruction_info.resize(new_poly_index_in_rec_basis + 1);
	new_poly.impl_->reconstruction_info[new_poly_index_in_rec_basis].emplace_back();
	//set signature to new value, just before poly_before
	//monomial is the same, since poly_before is already multiplied polynomial
	new_poly.impl_->sig.mon = poly_before.impl_->sig.mon;
	const auto greater_idx_pos = impl_->all_indices.find(new_poly.impl_->sig.index);
	assert(greater_idx_pos != impl_->all_indices.end());
	SigIdx smaller_idx = 
		(greater_idx_pos == impl_->all_indices.begin()) ? 
			std::numeric_limits<SigIdx>::min() : 
			*std::prev(greater_idx_pos);
	new_poly.impl_->sig.index  = smaller_idx + (*greater_idx_pos - smaller_idx)/2;
	//those assertions can fail during correct operation of an algorithm if too much additions were performed in same place
	assert(new_poly.impl_->sig.index > smaller_idx);
	assert(new_poly.impl_->sig.index < *greater_idx_pos);
	impl_->all_indices.insert(new_poly.impl_->sig.index);
}

FastZ2SlowBasedRingBase::LPoly FastZ2SlowBasedRingBase::DequeueSigSmallest(MultLPolysQueue& queue)
{
	LessComparer ms_less{&FastZ2SlowBasedRingBase::MonomialLess, *this};
	auto& queue_impl = *(queue.impl_);
	assert(!queue_impl.empty());
	auto min_it = std::min_element(queue_impl.begin(), queue_impl.end(), ms_less);
	assert(min_it != queue_impl.end());
	LPoly result;
	result.impl_->sig.index = min_it->poly.impl_->sig.index;
	result.impl_->sig.mon = Mmul(min_it->poly.impl_->sig.mon, min_it->mul_by);
	result.impl_->value = Pmul(min_it->poly.impl_->value, min_it->mul_by);
	result.impl_->reconstruction_info.reserve(min_it->poly.impl_->reconstruction_info.size());
	for (auto rec_info_poly:min_it->poly.impl_->reconstruction_info) {
		result.impl_->reconstruction_info.push_back(Pmul(rec_info_poly, min_it->mul_by));
	}
	 queue_impl.erase(min_it);
	return result;
}

void FastZ2SlowBasedRingBase::ExtendQueueBySpairPartsAndFilterUnneeded(const LPolysResult& left_parts, const LPoly& right_part, MultLPolysQueue& queue)
{
	LessComparer ms_less{&FastZ2SlowBasedRingBase::MonomialLess, *this};
	assert(!IsZero(right_part));
	for (const auto& left_part:*left_parts.impl_) {
		if (IsZero(left_part)) {
			continue;
		}
		MultLPoly left(left_part);
		left.mul_by = ToLCMMultiplier(HM(left_part.impl_->value), HM(right_part.impl_->value));
		MultLPoly right(right_part);
		right.mul_by = ToLCMMultiplier(HM(right_part.impl_->value), HM(left_part.impl_->value));
		MultLPoly &new_lpoly = ms_less(left, right) ? right : left;
		bool was_supeseded;
		auto& queue_impl = *queue.impl_;
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


Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> FastZ2SlowBasedRingBase::FieldAgnosticReconstructionInfoPolysImpl(const LPoly& poly)
{
	auto poly_enumerator = FullRangeEnumerator(poly.impl_->reconstruction_info);
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

Enumerator<CrossRingInfo::PerVariableData> FastZ2SlowBasedRingBase::FieldAgnosticReconstructionInfoTopImpl(const LPoly& poly)
{
	return ConverterEnumeratorCFunc<STATIC_WITHTYPE_AS_TEMPLATE_PARAM(VarDataFromMapItem)>(FullRangeEnumerator(HM(poly.impl_->value)));
}

FastZ2SlowBasedRingBase::LPolysResult FastZ2SlowBasedRingBase::FillWithTrivialSyzygiesOfNonMultElements(const MultLPolysQueue& queue)
{
	LPolysResult result;
	LessComparer ms_less{&FastZ2SlowBasedRingBase::MonomialLess, *this};
	const auto& queue_impl = *queue.impl_;
	for(auto i0 = queue_impl.begin(); i0 != queue_impl.end(); ++i0) {
		assert(MDeg(i0->mul_by)  == 0 && MDeg(i0->poly.impl_->sig.mon) == 0) ;
		for(auto i1 = std::next(i0); i1 != queue_impl.end(); ++i1) {		
			LPoly syz_part[2];
			//value and reconstruction info are zero
			syz_part[0].impl_->sig.index = i0->poly.impl_->sig.index;
			syz_part[0].impl_->sig.mon = HM(i1->poly.impl_->value);
			syz_part[1].impl_->sig.index = i1->poly.impl_->sig.index;
			syz_part[1].impl_->sig.mon = HM(i0->poly.impl_->value);
			//mul_by is identity nomomial
			MultLPoly mult[2] = {{syz_part[0]}, {syz_part[1]}};
			int greater_sig_idx = ms_less(mult[0], mult[1]) ? 1 : 0;
			InsertInResult(std::move(syz_part[greater_sig_idx]), result);
		}
	}
	return result;
}

void FastZ2SlowBasedRingBase::InsertInResult(LPoly&& poly, LPolysResult& result)
{
	result.impl_->emplace_back(std::move(poly));
}

bool FastZ2SlowBasedRingBase::IsZero(const LPoly& poly)
{
	//TODO
	return true;
}

void FastZ2SlowBasedRingBase::Normalize(LPoly& poly)
{
}

FastZ2SlowBasedRingBase::MultLPolysQueue FastZ2SlowBasedRingBase::PutInQueueExtendLabeledPolysImpl(Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input)
{
	
	//TODO
	//TODO: fill all_indices with indices
	return MultLPolysQueue();
}

bool FastZ2SlowBasedRingBase::QueueEmpty(const MultLPolysQueue& queue)
{
	//TODO
	return true;
}

void FastZ2SlowBasedRingBase::ReduceCheckingSignatures(LPoly& poly, LPolysResult& reducers)
{
	//TODO
}

/*
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
		mp.poly.sig_index = SigIdx(current_poly_index) + 1;
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
