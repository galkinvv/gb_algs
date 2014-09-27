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
	struct MonomialMetadataWithoutOrder
	{
		int var_count;
	};
	
	template <MonomialOrder t_order>
	struct MonomialMetadata:MonomialMetadataWithoutOrder
	{
		static const MonomialOrder order = t_order;
	};
	
	struct PerVariableData
	{
		PerVariableData(int a_degree, int a_index):
			degree(a_degree),
			index(a_index)
		{}
		static PerVariableData FromDI(int a_degree, int a_index)
		{
			return PerVariableData(a_degree, a_index);
		}
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
		const Coef& coef()const
		{
			return coef_;
		}
	private:
	  const Coef& coef_;
	};

	template <class ContainerIterator, class Coef>
	std::ostream& operator << (std::ostream& s, const MonomialDataWithCoef<ContainerIterator, Coef>& data)
	{
		s << "(" << data.coef();
		if (data.end() != data.begin())
		{
			s << " " << OutputContainer(data, "*");
		}
		return s << ")";
	}

	typedef std::vector<int> DegreesContainer;
		
	template <class MonomialMetadata>
	class AddedVarInfo
	{
		typedef MonomialData<DegreesContainer::const_iterator> MonData;
		typedef decltype(std::declval<MonData>().begin()) VarIterator;
	public:
		AddedVarInfo(const MonomialMetadata& metadata, int addedMonomialIndex,  int addedPolynomialIndex):
			metadata_(metadata),
			data_(metadata.var_count),
			addedMonomialIndex_(addedMonomialIndex),
			addedPolynomialIndex_(addedPolynomialIndex)
		{
			assert(addedMonomialIndex  <= metadata.var_count); //check that all previous indices are known
		}
		
		VarIterator begin()const
		{
			return AsMonomialData().begin();
		}

		VarIterator end()const
		{
			return AsMonomialData().end();
		}
		const MonomialMetadata& Metadata()const
		{
			return metadata_;
		}		
		int AddedMonomialIndex()const
		{
			return addedMonomialIndex_;
		}
		int AddedPolynomialIndex()const
		{
			return addedPolynomialIndex_;
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
		friend std::ostream& operator << (std::ostream& s, const AddedVarInfo<MonomialMetadata2>& data);
		MonData AsMonomialData()const
		{
			return MonData(data_.begin(), data_.end());
		}
		const MonomialMetadata metadata_;
		DegreesContainer data_;
		const int addedMonomialIndex_;
		const int addedPolynomialIndex_;
	};

	template <class MonomialMetadata>
	std::ostream& operator << (std::ostream& s, const AddedVarInfo<MonomialMetadata>& data)
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
		MonomialCollectionBase(int expected_monomial_count, DegreesContainer& container, const MonomialMetadataWithoutOrder& monomial_metadata)
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
		template <class DegIterator>
		struct IteratorBase
		{
			IteratorBase(DegIterator a_position, const MonomialMetadataWithoutOrder& a_monomial_metadata):
				position(a_position), monomial_metadata(a_monomial_metadata)
			{}
			
			DegIterator position;
			std::reference_wrapper<const MonomialMetadataWithoutOrder> monomial_metadata;

			void operator++()
			{
				position = NextPosition(1);
			}
			
			bool operator !=(const IteratorBase& other)const
			{
				return position != other.position;
			}
		protected:
			DegIterator NextPosition(int steps)const
			{
				return position + monomial_metadata.get().var_count * steps;
			}
		};

		ArrayInterval interval;
		const MonomialMetadataWithoutOrder& monomial_metadata_;
		DegreesContainer& container_;
	};

	class MonomialCollectionPlain:protected MonomialCollectionBase
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
				return Item(this->position, this->NextPosition(1));
			}
			
			PseudoPointer<Item> operator->()const
			{
				return PseudoPointer<Item>(this->position, this->NextPosition(1));
			}
		};
		typedef Iterator<DegreesContainer::const_iterator> const_iterator;
		const_iterator begin()const
		{
			return const_iterator{container_.begin() + interval.begin, monomial_metadata_};
		}
		const_iterator end()const
		{
			return const_iterator{container_.begin() + interval.end, monomial_metadata_};
		}
		using MonomialCollectionBase::size;
	};
	
	inline std::ostream& operator << (std::ostream& s, const MonomialCollectionPlain& data)
	{
		return s << '(' << OutputContainer(data, "+") << ')';
	}

	template <class Coef>
	class MonomialCollectionWithCoef:protected MonomialCollectionBase
	{
		public:
		typedef std::vector<Coef> CoefContainer;
		MonomialCollectionWithCoef(CoefContainer& coef_container, int expected_monomial_count, DegreesContainer& container, const MonomialMetadataWithoutOrder& monomial_metadata):
			MonomialCollectionBase(expected_monomial_count, container, monomial_metadata), coef_data_(coef_container), first_coef_pos_(coef_data_.size())
		{
			coef_data_.reserve(first_coef_pos_ + expected_monomial_count);
		}
	
		template <class... TA>
		void MonomialAdditionDone(TA&&... args)
		{
			interval.end += monomial_metadata_.var_count;
			assert(interval.end <= int(container_.size()));
			assert(coef_data_.size() < coef_data_.capacity());
			coef_data_.emplace_back(std::forward<TA>(args)...);
		}
		
		template <class DegIterator, class CoefIterator>
		struct Iterator:IteratorBase<DegIterator>
		{
			typedef IteratorBase<DegIterator> Base;
			Iterator(CoefIterator a_coef_postition, DegIterator a_position, const MonomialMetadataWithoutOrder& a_monomial_metadata):
				Base(a_position, a_monomial_metadata), coef_postition_(a_coef_postition)
			{}
			typedef MonomialDataWithCoef<DegIterator, Coef> Item;
			
			Item operator*()const
			{
				return Item(this->position, this->NextPosition(1), *coef_postition_);
			}
			
			PseudoPointer<Item> operator->()const
			{
				return PseudoPointer<Item>(this->position, this->NextPosition(1), *coef_postition_);
			}
			
			void operator++()
			{
				Base::operator ++();
				++coef_postition_;
			}
		  private:
			CoefIterator coef_postition_;
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
		using MonomialCollectionBase::size;
	  private:
	   CoefContainer& coef_data_;
	   int first_coef_pos_;
	};
		
	template <class Coef>
	std::ostream& operator << (std::ostream& s, const MonomialCollectionWithCoef<Coef>& data)
	{
		return s << '(' << OutputContainer(data, "+") << ')';
	}

	template<class SpecificMonomialCollection>
	struct MonomialCollectionImpl:SpecificMonomialCollection
	{
		DECLARE_FORWARDING_CONSTRUCTOR(MonomialCollectionImpl, SpecificMonomialCollection)

		using SpecificMonomialCollection::MonomialAdditionDone;
		using SpecificMonomialCollection::AddVariable;
		using SpecificMonomialCollection::size;
	};

	template <class MonomialMetadata>
	struct VariableMapping: private DegreesContainer, private MonomialCollectionPlain
	{
		VariableMapping(const MonomialMetadata& oldMetaData, int newMonomialsCount)
			:MonomialCollectionPlain(newMonomialsCount, *this, oldMetaData)
			,oldMetadata_(oldMetaData)
		{}

		using MonomialCollectionPlain::MonomialAdditionDone;
		using MonomialCollectionPlain::AddVariable;
		using MonomialCollectionPlain::size;
		
		const MonomialMetadata& OldMetadata()const
		{
			return oldMetadata_;
		}

	  private:
		struct SingleMonIteratorPosition:MonomialCollectionPlain::const_iterator
		{
			SingleMonIteratorPosition(const DegreesContainer::const_iterator& begin, const MonomialMetadataWithoutOrder& a_monomial_metadata, int monomial_idx)
				:MonomialCollectionPlain::const_iterator(begin, a_monomial_metadata)
			{
				NextPosition(monomial_idx);
			}
		};
	  public:
		decltype(*std::declval<SingleMonIteratorPosition>()) operator[](int idx)const
		{
			assert(idx >= 0);
			assert(idx < size());
			return *SingleMonIteratorPosition(static_cast<const DegreesContainer*>(this)->begin(), oldMetadata_, idx);
		}

	  
	  private:
		const MonomialMetadata oldMetadata_;
	};
	
	template <class MonomialMetadata, class SpecificMonomialCollection>
	struct MonomialListListBase
	{
		MonomialListListBase(const MonomialMetadata& metadata)
			:metadata_(metadata)
		{}
		void AddVariable(const PerVariableData& data)
		{
			BuiltPolynomial().AddVariable(data);
		}
		DECLARE_FORWARDING_METHOD(void, MonomialAdditionDone, BuiltPolynomial())

		const MonomialMetadata& Metadata()const
		{
			return metadata_;
		}		
		const SpecificMonomialCollection* begin()const
		{
			return input_poly_infos_.data();
		}
		const SpecificMonomialCollection* end()const
		{
			return input_poly_infos_.data() + input_poly_infos_.size();
		}
		size_t size()const
		{
			return input_poly_infos_.size();
		}
		int TotalMonomials()const
		{
			return std::accumulate(input_poly_infos_.begin(), input_poly_infos_.end(), 0, CalcSizeSum);
		}
		static int CalcSizeSum(int res, const MonomialCollectionImpl<SpecificMonomialCollection>& poly)
		{
			return res + poly.size();
		}
	  protected:
		const MonomialMetadata metadata_;
		DegreesContainer degrees_;
		std::vector<MonomialCollectionImpl<SpecificMonomialCollection>> input_poly_infos_;
		MonomialCollectionImpl<SpecificMonomialCollection>& BuiltPolynomial()
		{
			assert(!input_poly_infos_.empty());
			return input_poly_infos_.back();
		}
	};
	
	template <class MonomialMetadata>
	struct MonomialListListWithTopInfo: private MonomialListListBase<MonomialMetadata, MonomialCollectionPlain>
	{
		typedef MonomialListListBase<MonomialMetadata, MonomialCollectionPlain> Base;
		MonomialListListWithTopInfo(const MonomialMetadata& metadata)
			:Base(metadata)
		{
			BeginPolynomialConstruction(kMonomialInPolynomialForTopInfo);
		}
		void TopInfoAdditionDone()
		{
			Base::MonomialAdditionDone();
			assert(Base::begin() + kPolynomialsForTopInfo == Base::end());
		}
		const MonomialCollectionPlain* begin()const
		{
			auto result = Base::begin();
			assert(result != Base::end());
			return result + kPolynomialsForTopInfo;
		}
		size_t size()const
		{
			assert(Base::size() >= kPolynomialsForTopInfo);
			return Base::size() - kPolynomialsForTopInfo;
		}

		void BeginPolynomialConstruction(int monomial_count)
		{
			assert(this->TotalMonomials() * this->metadata_.var_count == this->degrees_.size());
			this->input_poly_infos_.emplace_back(monomial_count, this->degrees_, this->metadata_);
		}

		using Base::end;
		using Base::Metadata;
		using Base::AddVariable;
		void MonomialAdditionDone()
		{
			assert((Base::begin() + kPolynomialsForTopInfo) != Base::end());
			Base::MonomialAdditionDone();
		}
		decltype(*(std::declval<Base>().begin()->begin())) TopInfo()const
		{
			auto result_poly = Base::begin();
			assert(result_poly != Base::end());
			return *(result_poly->begin());
		}
	  private:
		static const int kPolynomialsForTopInfo = 1;
		static const int kMonomialInPolynomialForTopInfo = 1;
	};
	
	template <class MonomialMetadata, class TField>
	struct MonomialListListWithCoef:private MonomialListListBase<MonomialMetadata, MonomialCollectionWithCoef<typename TField::Value>>
	{	
		typedef typename TField::Value Coef;
		typedef MonomialListListBase<MonomialMetadata, MonomialCollectionWithCoef<Coef>> Base;
		MonomialListListWithCoef(const MonomialMetadata& metadata, const TField& field)
			:Base(metadata)
			,field_(field)
		{}
		using Base::size;
		using Base::begin;
		using Base::end;
		using Base::AddVariable;
		using Base::MonomialAdditionDone;
		void BeginPolynomialConstruction(int monomial_count)
		{
			assert(this->TotalMonomials() * this->metadata_.var_count == this->degrees_.size());
			this->input_poly_infos_.emplace_back(coefs_, monomial_count, this->degrees_, this->metadata_);
		}
		using Base::Metadata;
		const TField& Field()const
		{
			return field_;
		};
	  private:
		TField field_;
		std::vector<Coef> coefs_;
	};


	template <class MonomialMetadata, class Coef>
	std::ostream& operator << (std::ostream& s, const MonomialListListWithCoef<MonomialMetadata, Coef>& data)
	{
		return s << '{' << OutputContainer(data, ", ") << '}';
	}

	template <class MonomialMetadata>
	std::ostream& operator << (std::ostream& s, const MonomialListListWithTopInfo<MonomialMetadata>& data)
	{
		return s << "{data:{" << OutputContainer(data, ", ") << "}, top:" << OutputContainer(data.TopInfo(), ", ") << "}";
	}

	template <class MonomialMetadata, class Coef>
	std::istream& operator >> (std::istream& s, MonomialListListWithCoef<MonomialMetadata, Coef>& data)
	{
		struct
		{
			MonomialListListWithCoef<MonomialMetadata, Coef>& data_;
			bool TryReadPoly(std::istream& s)
			{
				//data_.BeginPolynomialConstruction(.. read from stteam first..);
				//TODO
				return false;
			}
		} reader = {data};
		return s >> FrontCharReader<MonomialListListWithCoef<MonomialMetadata, Coef>, '{'>() >> InputContainer(MEM_BIND(reader, TryReadPoly, _1), ',') >> CharReader<'}'>();
	}

}
