#pragma once
#include "outputroutines.h"
#include "conversions.h"

namespace F4MPI{
void Normalize(PolynomSet& polys);
SPair MakeSPair(const CPolynomial& p1, const CPolynomial& p2);

/**
Удаляет из вектора неотмеченныые элементы.
Все элементы вектора \a A, для которых в соответсвующей позиции в векторе \a Mark стоит \c false удаляются.
Размер вектора уменьшается. Время работы не более O(A.size()*сложность swap для элементов A).
*/
template<class X>
void EraseAll(std::vector<X> &A, std::vector<bool> Mark)
{	
	int L = 0;
	int R = A.size()-1;
	int count = 0;
	for(int i = 0; i<Mark.size(); i++)
		if(Mark[i])
			count++;
	while(L < R)
	{		
		while(L < R && Mark[L])L++;
		while(L < R && !Mark[R])R--;
		if(L < R)
		{
			swap(A[L], A[R]);
			L++;
			R--;
		}
	}
	A.resize(count);
}

template <class Dividers>
bool CheckMonIsDivisibleBySome(const CMonomial& mon, const Dividers& dividers)
{
	for(auto i = dividers.begin();i!=dividers.end();++i)
	{
		if (mon.divisibleBy(*i))
		{
			//cout << "division found for " << mon.toString() << " by " << i->toString() << endl;
			return true;
		}
		else
		{
			//cout << "NOT division for " << mon.toString() << " by " << i->toString() << endl;
		}
	}
	return false;
}

void GetBasisTops(const PolynomSet &basis, std::vector<CMonomial>& basisTops);
void Preprocess (PolynomSet& polys, PolynomSet& reducers);
bool cmpForReduceBySize(const CPolynomial& a, const CPolynomial &b);
bool cmpForReduceByOrder(const CPolynomial& a, const CPolynomial &b);
void AutoReduceBasis(PolynomSet& basis, const F4AlgData* f4options);
void AutoReduceSetWithBasisFull(const PolynomSet& set, const F4AlgData* f4options, ReduceBySet& result);
bool cmpForUpdaters(const CPolynomial& a, const CPolynomial &b);
}
