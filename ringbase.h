#pragma once
template <class TRing>
struct IOData
{
	const TRing& in_ring;
	const typename TRing::PolysSet& in;
	TRing out_ring;
	typename TRing::PolysSet out;
	virtual ~IOData(){};
protected:
	TRing in_ring_;
	typename TRing::PolysSet in_;
	IOData():
		in_ring(in_ring_), in(in_)
	{}
};
