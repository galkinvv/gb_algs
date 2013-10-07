/** \file
процедуры для преобразование в текст.
*/

#include "outputroutines.h"
#include "conversions.h"

using namespace std;
namespace F4MPI{


void PrintPolynomSet(ostream& output,PolynomSet& polys, ParserVarNames* names){
	if(polys.empty()){
		output<<"empty\n";
		return;
	}
	//сортировка многочленов для удобства проверки ответа
	typedef set<CPolynomial> PolyPrintSet;
	PolyPrintSet sortedset(polys.begin(),polys.end());//этот set использует другую сортировку (operator<)
	bool first = true;
	for(PolyPrintSet::iterator i = sortedset.begin(); i!=sortedset.end(); ++i){		
		if (!first){
			output<<",\n";
		}else first=false;
		i->printPolynomial(output,names);
	}
	output<<"\n";
}

void PrintPolynomSetRaw(ostream& output,PolynomSet& polys, ParserVarNames* names){
	if(polys.empty()){
		output<<"empty\n";
		return;
	}

	bool first = true;
	for(PolynomSet::iterator i = polys.begin(); i!=polys.end(); ++i){		
		if (!first){
			output<<",\n";
		}else first=false;
		i->printPolynomial(output,names);
	}
	output<<"\n";
}


void PrintMatrixAsPolynomials(ostream& output,CMatrix& matrix, ParserVarNames* names){
	if (matrix.empty()){
		output<<"empty\n";
		return;
	}
	for(int i = 0; i<int(matrix.size()); ++i){
		CPolynomial P;
		MatrixRow R = matrix.getRow(i);
		rowToPolynomial(R,matrix.getMonomialMap(),P);
		P.printPolynomial(output,names);
		output<<"\n";
	}
}

//печать множества S-пар
void PrintSPairSet(ostream& output,SPairSet& pairs){
	if (pairs.empty()){
		output<<"empty\n";
		return;
	}
	for(SPairSet::iterator i = pairs.begin(); i!=pairs.end(); ++i){		
		output<<"{";
		i->first.printPolynomial(output);
		output<<", ";
		i->second.printPolynomial(output);
		output<<"}\n";
	}
}

// a вот это недокументированная функция
CMonomial getMonomial(int i, int j, int k, int l, int m, int n, int o, int p){
	vector<CMonomialBase::Deg> ints1;
	ints1.push_back(i);
	if(CMonomial::theNumberOfVariables>1)
		ints1.push_back(j);
	if(CMonomial::theNumberOfVariables>2)
		ints1.push_back(k);
	if(CMonomial::theNumberOfVariables>3)
		ints1.push_back(l);
	if(CMonomial::theNumberOfVariables>4)
		ints1.push_back(m);
	if(CMonomial::theNumberOfVariables>5)
		ints1.push_back(n);
	if(CMonomial::theNumberOfVariables>6)
		ints1.push_back(o);
	if(CMonomial::theNumberOfVariables>7)
		ints1.push_back(p);
	CMonomial M1(ints1);
	return M1;
}
} //namespace F4MPI
