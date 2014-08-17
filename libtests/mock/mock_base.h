#pragma once
#include <string>
#include <stdexcept>
#include <sstream>

template<class InData = const char *, class OutData = const char *>
struct AnwserData
{
	InData in;
	OutData out;
};

//selects answer from data based on in parameter. data is pointer to array which ends in element with out evaluating to false in boolean context
template<class TAnwserData>
auto SelectAnwser(const TAnwserData data[], decltype(data->in) in) -> decltype(data->out)
{
	for(; data->out; ++data)
	{
		if (data->in == in)
		{
			return data->out;
		}
	}
	std::ostringstream s;
	s << "Failed to find mock anwser for input " << in;
	throw std::logic_error(s.str());
}