/**
\file
Реализация внешних функций работы с MonomialMap
*/
#include "monomialmap.h"
#include "types.h"

namespace F4MPI{

void storeMonomialsFromPoly(MonomialMap &M, const CPolynomial &p){
	for(CPolynomial::m_const_iterator i = p.m_begin(); i!=p.m_end(); ++i)		
			M.storeMonomial(*i);
}



void storeNotProcessedMonomialsFromPoly(MonomialMap &M, MonomialMap& notThere, const CPolynomial& p){
	for(CPolynomial::m_const_iterator i = p.m_begin(); i!=p.m_end(); ++i)
		if(!notThere.containsMonomial(*i))
			M.storeMonomial(*i);
}
} //namespace F4MPI
