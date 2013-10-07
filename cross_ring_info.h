#pragma once
#include <vector>
#include <cassert>
#include "zero_skipping_indexed_iterator.h"

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
		const int degree;
		const int index;
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

	private:
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

	class MonomialCollection{
	  protected:
		MonomialCollection(int expected_monomial_count, DegreesContainer& container, const MonimailMetaDataWithoutOrder& monomial_metadata)
			:interval{int(container.size()), int(container.size())}
			,monomial_metadata_(monomial_metadata)
			,container_(container)
		{
			container_.resize(interval.end + expected_monomial_count*monomial_metadata.var_count);
		}
		void MonomialAdditionDone()
		{
			interval.end += monomial_metadata_.var_count;
			assert(interval.end <= int(container_.size()));
		}
		void AddVariable(const PerVariableData& data)
		{
			int  index = interval.end + data.index;
			assert(index <= int(container_.size()));
			auto& value_to_change = container_[index];
			assert(0 == value_to_change );
			value_to_change  = data.degree;
		}
		int size()const
		{
			return (interval.end - interval.begin)/monomial_metadata_.var_count;
		}
		template <class BaseIterator>
		struct Iterator
		{
			BaseIterator position;
			const MonimailMetaDataWithoutOrder& monomial_metadata;
			typedef MonomialData<BaseIterator> Item;
			Item operator*()
			{
				return Item(position, NextPosition());
			}
			
			void operator++()
			{
				position = NextPosition();
			}
			
			bool operator !=(const Iterator& other)
			{
				return position != other.position;
			}
		private:
			BaseIterator NextPosition()
			{
				return position + monomial_metadata.var_count;
			}
		};
	  public:
		typedef Iterator<DegreesContainer::const_iterator> const_iterator;
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

	struct MonomialCollectionImpl:MonomialCollection
	{
		template <class... TA>
		MonomialCollectionImpl(TA&&... args):
			MonomialCollection(std::forward<TA>(args)...)
		{}
		using MonomialCollection::MonomialAdditionDone;
		using MonomialCollection::AddVariable;
		using MonomialCollection::size;
	};

	template <class MonomialMetadata>
	struct BasisElementReconstructionInfo
	{
		BasisElementReconstructionInfo(const MonomialMetadata& metadata, int expected_monomial_count)
		:metadata_(metadata)
		,poly_info_(expected_monomial_count, degrees_, metadata_)
		{}

		void AddVariable(const PerVariableData& data)
		{
			poly_info_.AddVariable(data);
		}
		void MonomialAdditionDone()
		{
			poly_info_.MonomialAdditionDone();
		}
		const MonomialMetadata& MetaData()const
		{
			return metadata_;
		}
		MonomialCollectionImpl::const_iterator begin()const
		{
			return poly_info_.begin();
		}
		MonomialCollectionImpl::const_iterator end()const
		{
			assert(poly_info_.size() * metadata_.var_count == degrees_.size());
			return poly_info_.end();
		}
	private:
		const MonomialMetadata& metadata_;
		DegreesContainer degrees_;
		MonomialCollectionImpl poly_info_;
	};

	template <class MonomialMetadata>
	struct InputElementConstructionInfo
	{
		InputElementConstructionInfo(const MonomialMetadata& metadata)
			:metadata_(metadata)
		{}
		
	private:
		const MonomialMetadata& metadata_;
		DegreesContainer degrees_;
		std::vector<MonomialCollectionImpl> input_poly_infos_;
	};
}

