#include "reducebyset.h"
#include "commonpolyops.h"
#include <algorithm>
using namespace std;
namespace F4MPI
{

ReduceBySet::ReduceBySet(const PolynomSet& aReducers)
{
	//MEASURE_TIME_IN_BLOCK("ReduceBySet::ReduceBySet");
	initialReducers.assign(aReducers.begin(), aReducers.end());
	for(auto i = initialReducers.begin(); i != initialReducers.end(); ++i)
	{
		i->normalize();
	}
	//sort(initialReducers.rbegin(),initialReducers.rend(), cmpForReduceBySize);
	//sort(initialReducers.rbegin(),initialReducers.rend(), cmpForReduceByOrder);
	//sort(initialReducers.begin(),initialReducers.end(), cmpForReduceBySize);
	sort(initialReducers.begin(),initialReducers.end(), cmpForReduceByOrder);
}

void ReduceBySet::reduceTopWithCache(CPolynomial& polyToReduce)
{
	do {} while(reduceTopOnceWithCache(polyToReduce));
}

const CPolynomial* ReduceBySet::findTopReducer(const CPolynomial& poly, CMonomial& mulby)
{
	return GetReducerAndMulFromSearchCache(poly.HM(), mulby);
}

bool ReduceBySet::reduceTopOnceWithCache(CPolynomial& polyToReduce)
{
	if (polyToReduce.empty()) return false;
	const CMonomial cMon = polyToReduce.HM();
	/*
	auto itPos = fullReducedCache.find(MonomialPtr(cMon));
	if (itPos != fullReducedCache.end())
	{
		const auto reducedMonPoly = itPos->second;
		if (!reducedMonPoly) return false;
		auto reducer = *reducedMonPoly;
		reducer *= polyToReduce.HC();
		polyToReduce = polyToReduce + reducer;
		return true;
	}
	*/
	CMonomial mulby;
	const auto* reductorPtr = GetReducerAndMulFromSearchCache(polyToReduce.HM(), mulby);
	if (!reductorPtr) return false;
	ReduceSinglePoly(polyToReduce, *reductorPtr, mulby);
	return true;
}

void ReduceBySet::reduceWithSearchCache(CPolynomial& polyToReduce)
{
	int i = 0;
	while(i < int(polyToReduce.size()))
	{
		const auto& mon = polyToReduce.getMon(i);
		CMonomial mulby;
		const auto reducedMonPtr = GetReducerAndMulFromSearchCache(mon, mulby);
		if (!reducedMonPtr)
		{
			++i;
		}
		else
		{
			polyToReduce = polyToReduce.GetLinearComb(*reducedMonPtr, -polyToReduce.getCoeff(i), mulby);
		}
	}
}

const CPolynomial* ReduceBySet::GetReducerAndMulFromSearchCache(const CInternalMonomial& m, CMonomial& mul_by)
{
	const CPolynomial* result = 0;
	CMonomial cMon = m;
	auto itPos = reducerSearchCache.find(MonomialPtr(cMon));
	if (itPos != reducerSearchCache.end())
	{
		if (itPos->second)
		{
			bool divisible = m.tryDivide(itPos->second->HM(), mul_by);
			assert(divisible);
		}
		return itPos->second;
	}
	for(auto redIt = initialReducers.cbegin(); redIt != initialReducers.cend(); ++redIt)
	{
		if (m.tryDivide(redIt->HM(), mul_by))
		{
			result = &*redIt;
			break;
		}
	}

	//запомнить надо в любом случае
	knownMonomialContainer.push_back(cMon);
	reducerSearchCache.insert(make_pair(MonomialPtr(knownMonomialContainer.back()), result));
	return result;
}

void ReduceBySet::freeCache()
{
	knownMonomialContainer.clear();
	//knownReducedContainer.clear();
	reducerSearchCache.clear();
	//fullReducedCache.clear();
}

void ReduceBySet::appendReducerGreater(const CPolynomial& reducer)
{
	assert(initialReducers.empty() || initialReducers.back().HM().compareTo(reducer.HM()) < 0);
	initialReducers.push_back(reducer);
	//freeCache();
	const auto reducerPtr = &initialReducers.back();
	reducerPtr->normalize();
	const auto& hm = reducer.HM();
	//исправление кешей редукторов с учётом нового многочлена
	for(auto i = reducerSearchCache.begin(); i!=reducerSearchCache.end();++i)
	{
		if (i->second == 0 && i->first.ptr->divisibleBy(hm))
		{
			reducerSearchCache[i->first] = reducerPtr;
		}
	}
	//knownReducedContainer.clear();
	//fullReducedCache.clear();
}

const SortedReducersSet& ReduceBySet::allReducers()
{
	return initialReducers;
}

void ReduceSinglePoly(CPolynomial& polyToReduce, const CPolynomial& reductor, const CMonomial& mulby)
{
	polyToReduce = polyToReduce.GetLinearComb(reductor, -polyToReduce.HC() *  CModular::inverseMod(reductor.HC()), mulby);
}

}
