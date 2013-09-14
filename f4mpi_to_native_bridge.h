#pragma once
template <class TRing, template <class> class TAlgoT>
PolynomSet GBwithRingAlgo(PolynomSet F, const F4AlgData* /*f4options*/){
	auto io_data = TRing::Create(F);
	TAlgoT<TRing>::Do(*io_data);
	return TRing::ConvertResult(io_data);
}