/**
\file
Реализация алгоритма F4
*/
#include "f4main.h"
#include "commonpolyops.h"
#include "outputroutines.h"
#include "conversions.h"
#include "matrixinfoimpl.h"

using namespace std;
namespace F4MPI{

/**реализация алгоритма Update.
\param G промежуточный базис. В него добавляется новый многочлен и выкидываются некоторые старые
\param P набор рассматриваемых S-пар, который модифицируется в соответсвии с критериями Бухбергера
\param h новый многочлен, добавляемый в промежуточный базис
*/
void Update(PolynomSet& G, SPairSet& P, const CPolynomial& h)
{	
	//MEASURE_TIME_IN_BLOCK("Update");
	vector<int> D;
	D.reserve(G.size());	
	vector<int> Pnew;
	Pnew.reserve(G.size());
	CMonomial HMh = h.HM();
	CMonomial HMg1;
	CMonomial LCM_HMh_HMg1;
	CMonomial HMg2;
	CMonomial LCM_HMh_HMg2;
	for(int i = 0; i<G.size(); i++)
	{
		HMg1 = G[i].HM();
		bool Condition = false;
		bool gcdOne = false;
		if(CMonomial::gcd(HMh, HMg1).getDegree()==0)
		{
			Condition = true;
			gcdOne = true;
		}
		else{		
			Condition = true;
			LCM_HMh_HMg1 = CMonomial::lcm(HMh, HMg1);			
			for(int j = i+1; Condition && j<G.size(); ++j)
			{
				HMg2 = G[j].HM();
				LCM_HMh_HMg2 = CMonomial::lcm(HMh, HMg2);
				if(LCM_HMh_HMg1.divisibleBy(LCM_HMh_HMg2))
					Condition = false;
			}			
			for(int j = 0; Condition && j<D.size(); j++)
			{
				HMg2 = G[D[j]].HM();
				LCM_HMh_HMg2 = CMonomial::lcm(HMh, HMg2);
				if(LCM_HMh_HMg1.divisibleBy(LCM_HMh_HMg2))
					Condition = false;
			}
		}
		if(Condition)
		{
			D.push_back(i);			
			if(!gcdOne)
				Pnew.push_back(i);				
		}
	}	
	vector<bool> Mark(P.size()+Pnew.size());
	fill(Mark.begin(), Mark.end(), true);

	{
		//MEASURE_TIME_IN_BLOCK("SecondCriteria");
		CMonomial LCM_HMg1_HMg2;
		for(int i = 0; i<P.size(); i++)
		{
			LCM_HMg1_HMg2 = CMonomial::lcm(P[i].first.HM(), P[i].second.HM());
			if( LCM_HMg1_HMg2.divisibleBy(HMh) &&        
				CMonomial::lcm(HMh, P[i].first.HM())!=LCM_HMg1_HMg2 &&
				CMonomial::lcm(HMh, P[i].second.HM())!=LCM_HMg1_HMg2)
			{
				Mark[i] = false;
				continue;
			}		
		}
	}

	{
		//MEASURE_TIME_IN_BLOCK("MakeNewSPairs");
		for(int i = 0; i<Pnew.size(); i++)
			P.push_back(MakeSPair(h, G[Pnew[i]]));
	}
	{
		//MEASURE_TIME_IN_BLOCK("EraseMarkedFromP");
		EraseAll<SPair>(P, Mark);
	}
	{
		//MEASURE_TIME_IN_BLOCK("CheckDivisibility");
		Mark.resize(G.size());
		for(int i = 0; i<G.size(); i++)
		{
			if(!G[i].HM().divisibleBy(HMh))
				Mark[i] = true;
			else
				Mark[i] = false;
		}
	}
	{
		//MEASURE_TIME_IN_BLOCK("EraseMarkedFromG");
		G.push_back(h);
		Mark.push_back(true);
		EraseAll<CPolynomial>(G, Mark);
	}
}

	
///Приводит матрицу \a m к ступенчатому/сильно ступенчатому виду в соответствии с \a f4options
void doReduceMatrix(CMatrix& m, const F4AlgData* f4options){
	if(f4options->diagonalEachStep){
		m.toDiagonalNormalForm(f4options);
	}else{
		m.toRowEchelonForm(f4options);
	}
}

///Редуцирует матрицу, собирая статистику по ней при необходимости
void ReduceMatrix(CMatrix& m, const F4AlgData* f4options){
	f4options->stats->totalNumberOfReducedMatr++;
	m.doMatrixStatsPre(f4options);
	doReduceMatrix(m,f4options);
	m.doMatrixStatsPost(f4options);
	if (f4options->mpi_start_info.isMainProcess() && f4options->showInfoToStdout){
		printf("%d matrices complete\n", f4options->stats->totalNumberOfReducedMatr);
		fflush(stdout);
	}
}

/**
Подготавливает S-пару к обработке.
Домножает многочлены S-пары на (минимально возможные) мономы таким образом, чтоб старшие их мономы стали равны друг другу
*/
void SPolynomial2(SPair& sp)
{
	CMonomial lcm = CMonomial::lcm(sp.first.HM(), sp.second.HM());
	CMonomial M1;
	CMonomial M2;
	
	lcm.tryDivide(sp.first.HM(), M1);
	lcm.tryDivide(sp.second.HM(), M2);
	sp.first*= M1;
	sp.second*= M2;	
}

/**
Выбирает S-пары для рассмотрения на следующем шаге.
На основе S-пар выбранных из множества \a sPairs формируются многочлены, записываемые в \a ret.
Из исходного множеcтва S-пар выбранные выкидываются.
*/
void SelectSPairs(SPairSet &sPairs, PolynomSet& ret)
{		
	//MEASURE_TIME_IN_BLOCK("SelectSPairs");
	vector<bool> Mark(sPairs.size());
	int minDeg = 2000000000;
	
	for(int i = 0; i<sPairs.size(); i++)
	{
		int d = CMonomial::lcm(sPairs[i].first.HM(), sPairs[i].second.HM()).getDegree();	
		if(d<minDeg)minDeg = d;
	}
	int count = 0;

	for(int i = 0; i<sPairs.size(); i++)
	{		
		if(CMonomial::lcm(sPairs[i].first.HM(), sPairs[i].second.HM()).getDegree() == minDeg)
		{		
			count++;
			Mark[i] = false;
		}
		else
			Mark[i] = true;
	}	
	for(int i = 0; i<sPairs.size(); i++)if(!Mark[i])
	{
		SPolynomial2(sPairs[i]);
		ret.push_back(sPairs[i].first);
		ret.push_back(sPairs[i].second);
	}
	EraseAll<SPair>(sPairs, Mark);		
}

/**реализация алгоритма Reduce (см теоретическую документацию).
\param polysToReduce множество многочленов, представляющее S-пары, которое нужно редуцировать.
Значение аргумента после возврата из процедуры неопределено (портится)
\param reducers множество многочленов-редукторов
\param result место для записи результата
\param f4options параметры F4: порядок сортировки многочленов перед помещением в матрицу и параметры матричных операций.
*/
void ReduceF4(PolynomSet& polysToReduce, PolynomSet& reducers, PolynomSet& result, const F4AlgData* f4options)
{	
	//MEASURE_TIME_IN_BLOCK("Reduce");
	Preprocess(polysToReduce, reducers);
	
	MonomialMap preprocessedHM;

	{
		//MEASURE_TIME_IN_BLOCK("storeMonomial");
		for(PolynomSet::iterator i = polysToReduce.begin(); i!=polysToReduce.end(); ++i)
		{
			preprocessedHM.storeMonomial(i->HM());
		}
	}

	CMatrix mainMatrix;

	{
		//MEASURE_TIME_IN_BLOCK("sort");
		if(f4options->useSizesForSelectingRow){
			sort(polysToReduce.begin(),polysToReduce.end(),cmpForReduceBySize);
		}else{
			sort(polysToReduce.begin(),polysToReduce.end(),cmpForReduceByOrder);
		}
	}
	{
		//MEASURE_TIME_IN_BLOCK("polyToMatrix");
		polyToMatrix(polysToReduce, mainMatrix);
	}

	{
		//MEASURE_TIME_IN_BLOCK("reduceMatrix");
		ReduceMatrix(mainMatrix, f4options);
	}

	{
		//MEASURE_TIME_IN_BLOCK("matrixToPoly");
		matrixToPoly(mainMatrix, result, preprocessedHM);
	}
}


/**реализация алгоритма F4
\param F множество многочленов, для которого нужно найти базис
\param f4options параметры различных этапов алгоритма и сбора статистики.
\retval базис Грёбнера для \a F
*/
PolynomSet F4(PolynomSet F, const F4AlgData* f4options){
	//MEASURE_TIME_IN_BLOCK("F4");
	PolynomSet basis;
	basis.reserve(100000);
	SPairSet sPairs;
	sPairs.reserve(100000);
	for(PolynomSet::iterator i = F.begin(); i!=F.end(); ++i)
	{	
		Update(basis, sPairs, *i);
	}
	PolynomSet newBasisElements;
	newBasisElements.reserve(100000);
	PolynomSet sPolynomials;
	sPolynomials.reserve(100000);
		
	while(!sPairs.empty())
	{		
		// Selecting SPairs		
		sPolynomials.clear();
		SelectSPairs(sPairs, sPolynomials);
		
		newBasisElements.clear();		
		ReduceF4(sPolynomials, basis, newBasisElements, f4options);		
		sort(newBasisElements.begin(),newBasisElements.end(),cmpForUpdaters);
		// Updating basis and sPairs
		for(int i = 0; i<newBasisElements.size(); i++)		
			Update(basis, sPairs, newBasisElements[i]);				
		
	}
	if(f4options->autoReduceBasis){
		AutoReduceBasis(basis, f4options);
	}
	Normalize(basis);
	return basis;
}

} //namespace F4MPI
