#pragma once
#include "utils.h"
template <class TMultLPoly, class TMonomial = decltype(std::declval<TMultLPoly>().mul_by)>
TMonomial  MultSig(const TMultLPoly& mp)
{
	return Mmul(mp.mul_by, mp.poly->sig.mon);
}

template <class TMultLPoly, class TMonomial = decltype(std::declval<TMultLPoly>().mul_by)>
TMonomial  MultHM(const TMultLPoly& mp)
{
	return Mmul(mp.mul_by, HM(mp.poly->value));
}

template <class TMultLPoly, class MonCompare>
bool SigLess(const TMultLPoly& mp1, const TMultLPoly& mp2, MonCompare monLess)
{
	bool item_less;
	if (unequal(mp1.poly->sig.index, mp2.poly->sig.index, item_less)) {
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
	if (maybe_supded.poly->sig.index !=  sup_by.poly->sig.index) {
		return false; //can't supersede wth different indexes
	}
	auto sig_supded = MultSig(maybe_supded);
	auto sig_by = MultSig(sup_by);
	if (!DivideIfCan(sig_supded, sig_by)) {
		return false; //index_old not divisible by index_new
	}
	if (IsZeroImpl(sup_by.poly->value)) return true; //zero polynomial with dividing sig
	if (IsZeroImpl(maybe_supded.poly->value)) return false; //zero polynomial with sig divisible by non-zero

	auto sig_supded_hm_by = Mmul(sig_supded, MultHM(sup_by));
	auto sig_by_hm_supded = Mmul(sig_by, MultHM(maybe_supded));
	if (monLess(sig_by_hm_supded, sig_supded_hm_by)) {
		return false; //HM/S for maybe_supded is smaller
	}
	return true; // HM/S for maybe_supded is not smaller
}

