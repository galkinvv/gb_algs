#pragma once
#include <vector>
#include <algorithm>
#include <cassert>
#include <ostream>
#include "zero_skipping_indexed_iterator.h"
#include "utils.h"

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
	std::ostream& operator << (std::ostream& s, const PerVariableData& data)
	{
		s << "x_" << data.index << "^" << data.degree;
		return s;
	}
	
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
		const ContainerIterator begin_;
		const ContainerIterator end_;
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
			assert(data.index >= 0);
			assert(data.index < monomial_metadata_.var_count);
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
			Item operator*()const
			{
				return Item(position, NextPosition());
			}
			
			PseudoPointer<Item> operator->()const
			{
				return PseudoPointer<Item>(position, NextPosition());
			}

			void operator++()
			{
				position = NextPosition();
			}
			
			bool operator !=(const Iterator& other)const
			{
				return position != other.position;
			}
		private:
			BaseIterator NextPosition()const
			{
				return position + monomial_metadata.var_count;
			}
		};
	  public:
		typedef Iterator<DegreesContainer::const_iterator> const_iterator;
		const_iterator begin()const
		{
			return const_iterator{container_.begin() + interval.begin, monomial_metadata_};
		}
		const_iterator end()const
		{
			return const_iterator{container_.begin() + interval.end, monomial_metadata_};
		}
	private:
		ArrayInterval interval;
		const MonimailMetaDataWithoutOrder& monomial_metadata_;
		DegreesContainer& container_;
	};

	std::ostream& operator << (std::ostream& s, const MonomialCollection& data)
	{
		s << '[';
		OutputContainer(s, data, ", ");
		s << ']';
		return s;
	}

	struct MonomialCollectionImpl:MonomialCollection
	{
		DECLARE_FORWARDING_CONSTRUCTOR(MonomialCollectionImpl, MonomialCollection)

		using MonomialCollection::MonomialAdditionDone;
		using MonomialCollection::AddVariable;
		using MonomialCollection::size;
	};

	template <class MonomialMetadata>
	struct BasisElementReconstructionInfo
	{
		BasisElementReconstructionInfo(const MonomialMetadata& metadata, int monomial_count)
		:metadata_(metadata)
		,poly_info_(monomial_count, degrees_, metadata_)
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
	std::ostream& operator << (std::ostream& s, const BasisElementReconstructionInfo<MonomialMetadata>& data)
	{
		s << '[';
		OutputContainer(s, data, ", ");
		s << ']';
		return s;
	}

	template <class MonomialMetadata>
	struct InputElementsConstructionInfo
	{
		InputElementsConstructionInfo(const MonomialMetadata& metadata)
			:metadata_(metadata)
		{}
		void BeginPolynomialConstruction(int monomial_count)
		{
			assert(TotalMonomials() * metadata_.var_count == degrees_.size());
			input_poly_infos_.emplace_back(monomial_count, degrees_, metadata_);
		}
		void AddVariable(const PerVariableData& data)
		{
			BuiltPolynomial().AddVariable(data);
		}
		void MonomialAdditionDone()
		{
			BuiltPolynomial().MonomialAdditionDone();
		}
		const MonomialMetadata& MetaData()const
		{
			return metadata_;
		}
		
		const MonomialCollection* begin()const
		{
			return input_poly_infos_.data();
		}
		const MonomialCollection* end()const
		{
			return input_poly_infos_.data() + input_poly_infos_.size();
		}

	private:
		int TotalMonomials()const
		{
			return std::accumulate(input_poly_infos_.begin(), input_poly_infos_.end(), 0, CalcSizeSum);
		}
		static int CalcSizeSum(int  res, const MonomialCollectionImpl& poly)
		{
			return res + poly.size();
		}
		const MonomialMetadata& metadata_;
		DegreesContainer degrees_;
		std::vector<MonomialCollectionImpl> input_poly_infos_;
		MonomialCollectionImpl& BuiltPolynomial()
		{
			assert(!input_poly_infos_.empty());
			return input_poly_infos_.back();
		}
	};
	
	template <class MonomialMetadata>
	std::ostream& operator << (std::ostream& s, const InputElementsConstructionInfo<MonomialMetadata>& data)
	{
		s << '[';
		OutputContainer(s, data, ", ");
		s << ']';
		return s;
	}


}

