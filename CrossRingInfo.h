namespace CrossRingInfo
{
	class enum MonomialOrder
	{
		DegRevLex
	};
	struct MonimailMetaDataWithoutOrder
	{
		int var_count;
	}
	template <MonomialOrder t_order>
	struct MonimailMetaData:MonimailMetaDataWithoutOrder
	{
		static const MonomialOrder order = t_order;
	}
	struct MonomialData{
		struct VariableData
		{
			int var_index, var_deg;
		}
		template<class Container>
		struct zero_skipping_const_iterator
		{
			explicit zero_skipping_const_iterator(const Container& cont):
				cont_iterator*(cont.cbegin())
				
			//iterate overzerovalues
			void operator++()
			{
				GoNext();
				GoToNextNonzero();
			}
			bool operator !=(const zero_skipping_const_iterator& other)
			{
				return items_left_to_end != other.items_left_to_end;
			}
		private:
		    void GoNext()
			{
					++cont_iterator;
					--items_left_to_end;
			}
			void GoToNextNonzero()
			{
				do{
					if (items_left_to_end == left_at_end)
					{
						break;
					}
					GoNext();
				}while(!*cont_iterator);
				assert(items_left_to_end > = kLeftAtEnd);
				assert(items_left_to_end == kLeftAtEnd || bool(*cont_iterator));
			}
			const int kLeftAtEnd= 0;
			typename Container::const_reverse_iterator cont_iterator;
			int items_left_to_end;
		}
		MonomialData(const MonimailMetaDataWithoutOrder& metadata)
		{
			data.resize(metadata.var_count);
		}
		
	protected:
		std::vector<int> data;
		template <class TMonimailMetaData>
		MonomialData(const MonomialData&) = delete;
		MonomialData operator=(const MonomialData&) = delete;
	};
	class MonomialCollection: private std::vector<Monomial>{
		
	};
	class BasisElementReconstructionInfo
	{
		
	};
	class InputElementConstructionInfo
	{
		
	};
}