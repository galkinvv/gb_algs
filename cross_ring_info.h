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
	
	inline std::ostream& operator << (std::ostream& s, const PerVariableData& data)
	{
		return s << "x_" << data.index << "^" << data.degree;
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

	template <class ContainerIterator>
	std::ostream& operator << (std::ostream& s, const MonomialData<ContainerIterator>& data)
	{
		return s << "(" << OutputContainer(data, "*") << ")";
	}

	template <class ContainerIterator, class Coef>
	struct MonomialDataWithCoef:MonomialData<ContainerIterator>{
		typedef MonomialData<ContainerIterator> Base;
		MonomialDataWithCoef(const ContainerIterator& begin, const ContainerIterator& end, const Coef& coef)
			:Base(begin,end)
			,coef_(coef)
		{
		}
		const Coef& coef()
		{
			return coef_;
		}
	private:
	  const Coef& coef_;
	};

	template <class ContainerIterator, class Coef>
	std::ostream& operator << (std::ostream& s, const MonomialDataWithCoef<ContainerIterator, Coef>& data)
	{
		return s << "(" << data.coef() << " * " << OutputContainer(data, "*") << ")";
	}

	typedef std::vector<int> DegreesContainer;
		
	template <class MonomialMetadata>
	class SingleMonomial
	{
		typedef MonomialData<DegreesContainer::const_iterator> MonData;
		typedef decltype(std::declval<MonData>().begin()) VarIterator;
	public:
		SingleMonomial(const MonomialMetadata& metadata):
			metadata_(metadata),
			data_(metadata.var_count)
		{}
		
		VarIterator begin()const
		{
			return AsMonomialData().begin();
		}

		VarIterator end()const
		{
			return AsMonomialData().end();
		}
		const MonomialMetadata& MetaData()const
		{
			return metadata_;
		}		
		void AddVariable(const PerVariableData& data)
		{
			assert(data.index >= 0);
			assert(data.index < metadata_.var_count);
			auto& value_to_change = data_[data.index];
			assert(0 == value_to_change );
			value_to_change  = data.degree;
		}
	private:
		template <class MonomialMetadata2>
		friend std::ostream& operator << (std::ostream& s, const SingleMonomial<MonomialMetadata2>& data);
		MonData AsMonomialData()const
		{
			return MonData(data_.begin(), data_.end());
		}
		const MonomialMetadata metadata_;
		DegreesContainer data_;
	};

	template <class MonomialMetadata>
	std::ostream& operator << (std::ostream& s, const SingleMonomial<MonomialMetadata>& data)
	{
		return s << data.AsMonomialData();
	}

	typedef int ArrayPos;
	struct ArrayInterval
	{
		ArrayPos begin, end;
		std::size_t size(){
			return end - begin;
		}
	};

	class MonomialCollectionBase{
	  protected:
		MonomialCollectionBase(int expected_monomial_count, DegreesContainer& container, const MonimailMetaDataWithoutOrder& monomial_metadata)
			:interval{int(container.size()), int(container.size())}
			,monomial_metadata_(monomial_metadata)
			,container_(container)
		{
			container_.resize(interval.end + expected_monomial_count*monomial_metadata.var_count);
		}
		void AddVariable(const PerVariableData& data)
		{
			assert(data.index >= 0);
			assert(data.index < monomial_metadata_.var_count);
			int  index = interval.end + data.index;
			assert(index < int(container_.size()));
			auto& value_to_change = container_[index];
			assert(0 == value_to_change );
			value_to_change  = data.degree;
		}
		int size()const
		{
			return (interval.end - interval.begin)/monomial_metadata_.var_count;
		}
	protected:
		template <class DegIterator>
		struct IteratorBase
		{
			IteratorBase(DegIterator a_position, const MonimailMetaDataWithoutOrder& a_monomial_metadata):
				position(a_position), monomial_metadata(a_monomial_metadata)
			{}
			
			DegIterator position;
			const MonimailMetaDataWithoutOrder& monomial_metadata;

			void operator++()
			{
				position = NextPosition();
			}
			
			bool operator !=(const IteratorBase& other)const
			{
				return position != other.position;
			}
		protected:
			DegIterator NextPosition()const
			{
				return position + monomial_metadata.var_count;
			}
		};

		ArrayInterval interval;
		const MonimailMetaDataWithoutOrder& monomial_metadata_;
		DegreesContainer& container_;
	};

	class MonomialCollection:protected MonomialCollectionBase
	{
		public:
		using MonomialCollectionBase::MonomialCollectionBase;
		void MonomialAdditionDone()
		{
			interval.end += monomial_metadata_.var_count;
			assert(interval.end <= int(container_.size()));
		}
		
		template <class DegIterator>
		struct Iterator:IteratorBase<DegIterator>
		{
			using IteratorBase<DegIterator>::IteratorBase;
			typedef MonomialData<DegIterator> Item;
			
			Item operator*()const
			{
				return Item(this->position, this->NextPosition());
			}
			
			PseudoPointer<Item> operator->()const
			{
				return PseudoPointer<Item>(this->position, this->NextPosition());
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
	};
	
	template <class Coef>
	class MonomialCollectionWithCoef:protected MonomialCollectionBase
	{
		public:
		typedef std::vector<Coef> CoefContainer;
		MonomialCollectionWithCoef(CoefContainer& coef_container, int expected_monomial_count, DegreesContainer& container, const MonimailMetaDataWithoutOrder& monomial_metadata):
			MonomialCollectionBase(expected_monomial_count, container, monomial_metadata), coef_data_(coef_container), first_coef_pos_(coef_data_.size())
			{
				coef_data_.reserve(first_coef_pos_ + expected_monomial_count);
			}
		void MonomialAdditionDone(const Coef& coef)
		{
			interval.end += monomial_metadata_.var_count;
			assert(interval.end <= int(container_.size()));
			assert(coef_data_.size() < coef_data_.capacity());
			coef_data_.push_back(coef);
		}
		
		template <class DegIterator, class CoefIterator>
		struct Iterator:IteratorBase<DegIterator>
		{
			typedef IteratorBase<DegIterator> Base;
			Iterator(CoefIterator a_coef_postition, DegIterator a_position, const MonimailMetaDataWithoutOrder& a_monomial_metadata):
				Base(a_position, a_monomial_metadata), coef_postition(a_coef_postition)
			{}
			typedef MonomialDataWithCoef<DegIterator, Coef> Item;
			
			Item operator*()const
			{
				return Item(this->position, this->NextPosition(), *coef_postition);
			}
			
			PseudoPointer<Item> operator->()const
			{
				return PseudoPointer<Item>(this->position, this->NextPosition(), *coef_postition);
			}
			
			void operator++()
			{
				Base::operator ++();
				++coef_postition;
			}

			CoefIterator coef_postition;
		};
	public:
		typedef Iterator<DegreesContainer::const_iterator, typename CoefContainer::const_iterator> const_iterator;
		const_iterator begin()const
		{
			return const_iterator{coef_data_.begin() + first_coef_pos_, container_.begin() + interval.begin, monomial_metadata_};
		}
		const_iterator end()const
		{
			return const_iterator{coef_data_.begin() + first_coef_pos_, container_.begin() + interval.end, monomial_metadata_};
		}
	  private:
	   CoefContainer& coef_data_;
	   int first_coef_pos_;
	};
	
	
	inline std::ostream& operator << (std::ostream& s, const MonomialCollection& data)
	{
		return s << '(' << OutputContainer(data, "+") << ')';
	}

	struct MonomialCollectionImpl:MonomialCollection
	{
		DECLARE_FORWARDING_CONSTRUCTOR(MonomialCollectionImpl, MonomialCollection)

		using MonomialCollection::MonomialAdditionDone;
		using MonomialCollection::AddVariable;
		using MonomialCollection::size;
	};

	template <class MonomialMetadata>
	struct MonomialListListBase
	{
		protected:
		MonomialListListBase(const MonomialMetadata& metadata)
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

		DECLARE_FORWARDING_METHOD(void, MonomialAdditionDone, BuiltPolynomial())

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
		const MonomialMetadata metadata_;
		DegreesContainer degrees_;
		std::vector<MonomialCollectionImpl> input_poly_infos_;
		MonomialCollectionImpl& BuiltPolynomial()
		{
			assert(!input_poly_infos_.empty());
			return input_poly_infos_.back();
		}
	};
	
	template <class MonomialMetadata>
	struct MonomialListList:private MonomialListListBase<MonomialMetadata>
	{	
		typedef MonomialListListBase<MonomialMetadata> Base;
		DECLARE_FORWARDING_CONSTRUCTOR(MonomialListList, Base);
		using Base::begin;
		using Base::end;
		using Base::AddVariable;
		using Base::MonomialAdditionDone;
		using Base::BeginPolynomialConstruction;
		using Base::MetaData;
	};

	
	template <class MonomialMetadata>
	struct MonomialListListWithTopInfo: private MonomialListList<MonomialMetadata>
	{
		typedef  MonomialListList<MonomialMetadata> Base;
		MonomialListListWithTopInfo(const MonomialMetadata& metadata)
			:Base(metadata)
		{
			BeginPolynomialConstruction(1);
		}
		void TopInfoAdditionDone()
		{
			Base::MonomialAdditionDone();
			assert((Base::begin() + 1) == Base::end());
		}
		const MonomialCollection* begin()const
		{
			auto result = Base::begin();
			assert(result != Base::end());
			return result + 1;
		}
		using Base::end;
		using Base::MetaData;
		using Base::BeginPolynomialConstruction;
		using Base::AddVariable;
		void MonomialAdditionDone()
		{
			assert((Base::begin() + 1) != Base::end());
			Base::MonomialAdditionDone();
		}
		decltype(*(std::declval<Base>().begin()->begin())) TopInfo()const
		{
			auto result_poly = Base::begin();
			assert(result_poly != Base::end());
			return *result_poly->begin();
		}
	};
	
	template <class MonomialMetadata>
	std::ostream& operator << (std::ostream& s, const MonomialListList<MonomialMetadata>& data)
	{
		return s << '{' << OutputContainer(data, ", ") << '}';
	}

	template <class MonomialMetadata>
	std::ostream& operator << (std::ostream& s, const MonomialListListWithTopInfo<MonomialMetadata>& data)
	{
		return s << "{data:{" << OutputContainer(data, ", ") << "}, top:" << OutputContainer(data.TopInfo(), ", ") << "}";
	}

}

