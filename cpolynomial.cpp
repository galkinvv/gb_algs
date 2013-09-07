/**
\file
Реализация не-inline методов CPlainPolynomial
*/

#include "cpolynomial.h"
#include "cmodular.h"

#include <string>
#include <cmath>
using namespace std;
namespace F4MPI{

/**прибавление многочлена с коэффициентом.
В \a Result записывается многочлен равный \a p1 + \a p2 * \a addCoeff.
Исользуется для реализации операций + и -, поскольку их производительность не играет серьёзной роли.
*/
void CPlainPolynomial::addCPolynomialMultiplied(const CPlainPolynomial& p1, const CPlainPolynomial& p2, CModular addCoeff, CPlainPolynomial& Result){
	//MEASURE_TIME_IN_BLOCK("addCPolynomialMultiplied");
	int i1 = 0;
	int i2 = 0;
	int i1f = p1.size();
	int i2f = p2.size();

	CMonomial M;
	Result.reserve(i1f + i2f);

	while(i1!=i1f && i2!=i2f){			
		int dif = p1.getMon(i1).compareTo(p2.getMon(i2));
		CModular c;
		if(dif==0){				
			c = p1.getCoeff(i1) +addCoeff*p2.getCoeff(i2);
			M = p1.getMon(i1);
			++i1;
			++i2;
		}		
		else if(dif>0){
			c = p1.getCoeff(i1);
			M = p1.getMon(i1);				
			++i1;
		}		
		else{
			c = addCoeff*p2.getCoeff(i2);
			M = p2.getMon(i2);				
			++i2;
		}
		if(c!=0){
			Result.pushTermBack(c, M);
		}
	}
	while(i1!=i1f){			
		CModular c;
		c = p1.getCoeff(i1);
		M = p1.getMon(i1);				
		if(c!=0){
			Result.pushTermBack(c, M);
		}
		++i1;
	}

	while(i2!=i2f){			
		CModular c;
		c = addCoeff*p2.getCoeff(i2);
		M = p2.getMon(i2);
		if(c!=0){
			Result.pushTermBack(c, M);
		}
		++i2;
	}	
	if (Result.empty()) Result.clear();
}

void CPlainPolynomial::addCPolynomialMultipliedWithMonom(const CPlainPolynomial& p1, const CPlainPolynomial& p2, CModular addCoeff, const CMonomial& addMon, CPlainPolynomial& Result){
	//MEASURE_TIME_IN_BLOCK("addCPolynomialMultipliedWithMonom");
	if (addMon.isOne())
	{
		addCPolynomialMultiplied(p1, p2, addCoeff, Result);
		return;
	}
	int i1 = 0;
	int i2 = 0;
	int i1f = p1.size();
	int i2f = p2.size();

	CMonomial M, mon2;
	Result.reserve(i1f + i2f);

	if (i2!=i2f) mon2 = p2.getMon(i2) * addMon;
	while(i1!=i1f && i2!=i2f){
		int dif = p1.getMon(i1).compareTo(mon2);
		CModular c;
		if(dif==0){				
			c = p1.getCoeff(i1) +addCoeff*p2.getCoeff(i2);
			M = mon2;
			++i1;
			++i2;
			if (i2!=i2f) mon2 = p2.getMon(i2) * addMon;
		}		
		else if(dif>0){
			c = p1.getCoeff(i1);
			M = p1.getMon(i1);				
			++i1;
		}		
		else{
			c = addCoeff*p2.getCoeff(i2);
			M = mon2;				
			++i2;
			if (i2!=i2f) mon2 = p2.getMon(i2) * addMon;
		}
		if(c!=0){
			Result.pushTermBack(c, M);
		}
	}
	while(i1!=i1f){			
		CModular c;
		c = p1.getCoeff(i1);
		M = p1.getMon(i1);				
		Result.pushTermBack(c, M);
		++i1;
	}

	while(i2!=i2f){			
		CModular c;
		c = addCoeff * p2.getCoeff(i2);
		M = addMon * p2.getMon(i2);
		Result.pushTermBack(c, M);
		++i2;
	}
	if (Result.empty()) Result.clear();
}

/*
CPlainPolynomial CPlainPolynomial::operator+(const CPlainPolynomial& p)const{
	CPlainPolynomial Result;
	addCPolynomialMultiplied(*this,p,CModular(1),Result);
	return Result;	
}


CPlainPolynomial CPlainPolynomial::operator-(const CPlainPolynomial& p)const{
	CPlainPolynomial Result;
	addCPolynomialMultiplied(*this,p,CModular(-1),Result);
	return Result;	
}
*/

//Реализация не очень эффективна, но используется только в вводе данных
CPlainPolynomial CPlainPolynomial::operator*(const CPlainPolynomial& p)const{
	CPlainPolynomial result;
	for (int i=0;i<p.size();++i){
		CPlainPolynomial tcopy=*this;
		tcopy*=p.getMon(i);
		tcopy*=p.getCoeff(i);
		CPlainPolynomial added;
		addCPolynomialMultiplied(result, tcopy, CModular(1), added);
		result=added;
	}
	return result;	
}

void CPlainPolynomial::addTerm(CModular c, const CMonomial& m){
	m_iterator i=mons.begin();
	while (i!=mons.end() && (m.compareTo(*i) < 0)) ++i;
	coeffs.insert(coeffs.begin()+(i-mons.begin()),c);
	mons.insert(i,m);
}

void CPlainPolynomial::pushTermBack(CModular c, const CMonomial& m){
	mons.push_back(m);
	coeffs.push_back(c);
}

void CPlainPolynomial::printPolynomial(ostream& output, ParserVarNames* names) const {
	if (empty())
	{
		output<<"empty_poly";
		return;
	}
	int i=0;
	while(i!=size()){
		string term= "";
		string mon = getMon(i).toString(names);
		CModular coeff=getCoeff(i);
		//if(coeff == -1) term = "-";
		//bool nonOne=std::abs(coeff.toint())!=1;
		bool nonOne=coeff.toint()!=1;
		if(nonOne || mon=="") term = coeff.toString();
		if(nonOne && mon!=""){
			term += '*';
		}
		term+=mon;
		output<<term;
		++i;
		if(i!=size())
			output<<"+";
	}
}
} //namespace F4MPI
