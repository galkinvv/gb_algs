#pragma once
#include <unordered_map>
#include "types.h"
#include "monomialmap.h"
#include <deque>
namespace F4MPI{

struct ReduceBySet
{
	ReduceBySet(){}
	ReduceBySet(const PolynomSet& aReducers);
	void appendReducerGreater(const CPolynomial& reducer);
	const SortedReducersSet& allReducers();
	void reduceTopWithCache(CPolynomial& polyToReduce);
	const CPolynomial* findToPReduceZ2r(const CPolynomial& poly, CMonomial& mulby);
	bool reduceTopOnceWithCache(CPolynomial& polyToReduce);
	//void reduceWithFullCache(CPolynomial& polyToReduce);
	void reduceWithSearchCache(CPolynomial& polyToReduce);
	void freeCache();
private:
	const CPolynomial* GetReducerAndMulFromSearchCache(const CInternalMonomial& m, CMonomial& mulby);

	std::deque<CMonomial> knownMonomialContainer;
	//SortedReducersSet knownReducedContainer;
	SortedReducersSet initialReducers;
	std::unordered_map<MonomialPtr, const CPolynomial*, monomialHasher> reducerSearchCache;
};

void ReduceSinglePoly(CPolynomial& polyToReduce, const CPolynomial& reductor, const CMonomial& mulby);

}
