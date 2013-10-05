#pragma once
#include <vector>
#include <cassert>

namespace CrossRingInfo
{
	enum class MonomialOrder
	{
		DegRevLex
	};
	struct MonimailMetaDataWithoutOrder
	{
		int var_count;
	};
	template <MonomialOrder t_order>
	struct MonimailMetaData:MonimailMetaDataWithoutOrder
	{
		static const MonomialOrder order = t_order;
	};
	
	typedef std::vector<int> DegreesContainer;
	typedef int ArrayPos;
	struct ArrayInterval
	{
		ArrayPos begin, end;
		std::size_t size(){
			return end - begin;
		}
	};
	
	
	struct PerVariableData
	{
		PerVariableData(int a_degree, int a_index):
			degree(a_degree),
			index(a_index)
		{}
		int degree;
		int index;
	};
	
	template <class Container>
	struct MonomialData{
		MonomialData(Container& container, ArrayInterval pos)
			:container_(container)
			,pos_(pos)
		{
		}
		
		void Add(const PerVariableData& var)
		{
			assert(var.index >= 0);
			assert(var.index < pos.size());
			assert(var.degree > 0);
			auto& container_value = container[pos.begin + var.index];
			assert(container_value  == 0);
			container_value = var.degree;
		}
		
		typedef ZeroSkippingIndexedIterator<Container, PerVariableData> const_iterator;
		const_iterator begin() const
		{
			return const_iterator::BeginOf(container, pos_);
		}

		const_iterator end() const
		{
			return const_iterator::EndOf(container, pos_);
		}

	protected:
		Container& container_;
		ArrayInterval pos_;
		MonomialData(const MonomialData&) = delete;
		MonomialData operator=(const MonomialData&) = delete;
	};
	template <class Container>
	class MonomialCollection: private std::vector<MonomialData<Container>>{
		
	};
	typedef   
	class BasisElementReconstructionInfo
	{
		
	};
	class InputElementConstructionInfo
	{
		
	};
}

