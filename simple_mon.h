#pragma once
#include <map>
#include <unordered_set>
#include <memory>
#include <algorithm>

#include "utils.h"


struct SimpleMon : std::map<int,int> {
	friend bool operator<(const SimpleMon&, const SimpleMon&); //undefined
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
TPolynomial PAddZ2(const TPolynomial& not_muled_item, const TPolynomial& muled_item, const TMonomial& mul_by)
{
	std::unordered_set<TMonomial, decltype(SmallCollectionHash)> presence_count;
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
TPolynomial PReduceZ2(const TPolynomial& poly_to_red, const TPolynomial& by, const TMonomial& mul_by)
{
	assert(HM(poly_to_red) == Mmul(HM(by), mul_by));
	return PAddZ2(poly_to_red, by, mul_by);
}
