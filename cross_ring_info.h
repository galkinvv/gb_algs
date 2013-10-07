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
	
	struct PerVariableData
	{
		PerVariableData(int a_degree, int a_index):
			degree(a_degree),
			index(a_index)
		{}
		int degree;
		int index;
	};
	
	template <class ContainerIterator>
	struct MonomialData{
		MonomialData(const ContainerIterator& begin, const ContainerIterator& end)
			:begin_(begin)
			,end_(end)
		{
		}
		
		typedef ZeroSkippingIndexedIterator<ContainerIterator, PerVariableData> const_iterator;
		const_iterator begin() const
		{
			return const_iterator::BeginOf(begin_, end_);
		}

		const_iterator end() const
		{
			return const_iterator::EndOf(begin_, end_);
		}

	protected:
		const ContainerIterator &begin_;
		const ContainerIterator &end_;
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

	struct MonomialCollection{
		MonomialCollection(int expected_monomial_count, DegreesContainer& container, const MonimailMetaDataWithoutOrder& monomial_metadata)
			:interval{container.size(), container.size()}
			,monomial_metadata_(monomial_metadata)
			,container_(container)
		{
			container_.resize(interval.end + expected_monomial_count*monomial_metadata.var_count);
		}
		void MonomialAdditionDone()
		{
			interval.end += monomial_metadata.var_count;
			assert(interval.end <= container_.size());
		}
		void AddVariable(const PerVariableData& data)
		{
			int  index = interval.end + data.index;
			assert(index <= container_.size());
			auto& value_to_change = container_[index];
			assert(0 == value_to_change );
			value_to_change  = data.degree;
		}
		template <class BaseIterator = DegreesContainer::const_iterator>
		struct const_iterator
		{
			BaseIterator position;
			const MonimailMetaDataWithoutOrder& monomial_metadata;
			
			Item operator*()
			{
				return MonomialData<BaseIterator>(position, NextPosition());
			}
			
			void operator++()
			{
				position = NextPosition();
			}
			
			bool operator !=(const const_iterator& other)
			{
				return position != other.position;
			}
		private:
			BaseIterator NextPosition()
			{
				return position + monomial_metadata.var_count;
			}
		};
		const_iterator begin()const
		{
			return const_iterator{container_.begin() + interval.begin * monomial_metadata_.var_count , monomial_metadata_};
		}
		const_iterator end()const
		{
			return const_iterator{container_.begin() + interval.end* monomial_metadata_.var_count , monomial_metadata_};
		}
	private:
		ArrayInterval interval;
		const MonimailMetaDataWithoutOrder& monomial_metadata_;
		DegreesContainer& container_;
	};
	
	
	class BasisElementReconstructionInfo
	{
		
	};
	class InputElementConstructionInfo
	{
		
	};
}

