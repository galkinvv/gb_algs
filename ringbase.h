#pragma once
#include "cross_ring_info.h"
template <class TMonomialMetadata, class TField>
struct Ring
{
	typedef TField Field;
	typedef TMonomialMetadata MonomialMetadata;

	struct IOData
	{
protected:
	Field field_;
	MonomialMetadata in_order_;
	CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field> in_data_;
	IOData():
		field(field_), in_order(in_order_), in_data(in_data_)
	{}
	
	public:
		const Field& field;
		const MonomialMetadata& in_order;
		const CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field>& in_data;
		MonomialMetadata out_order;
		CrossRingInfo::MonomialListListWithCoef<MonomialMetadata, Field> out_data;
		virtual ~IOData(){};
	};
};
