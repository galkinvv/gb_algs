#pragma once
template<class BaseIterator, class Item>
	struct ZeroSkippingIndexedIterator
	{
		static ZeroSkippingIndexedIterator BeginOf(const BaseIterator& begin, const BaseIterator& end)
		{
			return ZeroSkippingIndexedIterator(begin, end, begin);
		}
			
		static ZeroSkippingIndexedIterator EndOf(const BaseIterator& begin, const BaseIterator& end)
		{
			return ZeroSkippingIndexedIterator(begin, end, end);
		}

		Item operator*()
		{
			assert(!AtEnd());
			return Item(CurrentValue(), current_ - begin_);
		}
		
		//iterate over zero values
		void operator++()
		{
			GoNext();
			GoToNextNonzero();
		}
		
		bool operator !=(const ZeroSkippingIndexedIterator& other)
		{
			assert(other.begin_ == begin_);
			assert(other.end_ == end_);
			return other.current_ != current_;
		}
	private:
		ZeroSkippingIndexedIterator(const BaseIterator begin, const BaseIterator end, BaseIterator current)
			:begin_(begin), end_(end), current_(current)
		{
			GoToNextNonzero();
		}
		
		decltype(*(BaseIterator())) CurrentValue()
		{
			return *current_;
		}
		
		bool AtEnd()
		{
			return current_ == end_;
		}
		void GoNext()
		{
			++current_;
		}
		void GoToNextNonzero()
		{	
			while(!AtEnd())
			{
				if (CurrentValue())
				{
					break;
				}
				GoNext();
			}
		}
		BaseIterator  begin_;
		BaseIterator  end_;
		BaseIterator  current_;
	};
	