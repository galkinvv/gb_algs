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
	
	template<class Container, class Item>
	struct ZeroSkippingIndexedIterator
	{	
		ZeroSkippingIndexedIterator():
			index_(0), 
			container_(nullptr)
		{
		}
			
		explicit ZeroSkippingIndexedIterator(const Container& container):
			index_(0), 
			container_(&container)
		{
			GoToNextNonzero();
		}
			
		Item operator*()
		{
			assert(container_);
			assert(index_ < container_->size());
			return Item((*container_)[index_], index_);
		}
		
		//iterate over zero values
		void operator++()
		{
			GoNext();
			GoToNextNonzero();
		}
		
		bool operator !=(const ZeroSkippingIndexedIterator& other)
		{
			//assume compare with end;
			assert(!other.container_);
			assert(container_);
			return index_ != container_->size();
		}
	private:
		void GoNext()
		{
			++index_;
		}
		void GoToNextNonzero()
		{	
			assert(container_);
			while(index_ < container_->size())
			{
				if ((*container_)[index_])
				{
					break;
				}
				GoNext();
			}
		}
		int index_;
		const Container* container_;
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
	
	struct MonomialData{
		MonomialData(const MonimailMetaDataWithoutOrder& metadata)
		{
			data.resize(metadata.var_count);
		}
		
		void Add(const PerVariableData& var)
		{
			assert(var.index < data.size());
			assert(data[var.index] == 0);
			data[var.index] = var.degree;
		}
		typedef ZeroSkippingIndexedIterator<std::vector<int>, PerVariableData> const_iterator;
		const_iterator begin() const
		{
			return const_iterator(data);
		}

		const_iterator end() const
		{
			return const_iterator();
		}

	protected:
		std::vector<int> data;
		template <class TMonimailMetaData>
		MonomialData(const MonomialData&) = delete;
		MonomialData operator=(const MonomialData&) = delete;
	};
	class MonomialCollection: private std::vector<MonomialData>{
		
	};
	class BasisElementReconstructionInfo
	{
		
	};
	class InputElementConstructionInfo
	{
		
	};
}

