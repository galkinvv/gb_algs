#include "commonpolyops.h"
#include "reducebyset.h"
using namespace std;
namespace F4MPI
{
/**
Нормализует все многочлены множества.
Домножает каждый многочлен множества на множитель, обратный старшему коэффициенту этого многочлена
*/
void Normalize(PolynomSet& polys){
	for(PolynomSet::iterator i=polys.begin();i!=polys.end();++i){
		i->normalize();
	}
}

///Возвращает S-пару, составленную из \a p1 и \a p2
SPair MakeSPair(const CPolynomial& p1, const CPolynomial& p2)
{
	if(p1.compareTo(p2)<0)
		return SPair(p1, p2);
	return SPair(p2, p1);
}

void GetBasisTops(const PolynomSet &basis, vector<CMonomial>& basisTops)
{
	for(PolynomSet::const_iterator i = basis.begin();i!=basis.end();++i)
	{
		if (!i->empty())
		{
			CMonomial top = i->HM();
			basisTops.push_back(top);
		}
	}
}


///добавляет в множество \a Result все мономы всех многочленов из \a polys
void MonomialsOf(MonomialMap& Result, PolynomSet& polys)
{	
	for(PolynomSet::const_iterator i = polys.begin(); i!=polys.end(); ++i)
		storeMonomialsFromPoly(Result, *i);
}

/**функция сравнения для препроцессинга.
Перед препроцессингом многочлены упорядочиваются в соотвествии с этим сравнением,
чтоб при возможности использования нескольких многочленов, как препроцессирующих
был выбран оптимальный с точки зрения данного критерия.
Сейчас в этой функции используется критерий минимальности числа ненулевых мономов.
*/
bool cmpForPreprocess(const CPolynomial& a, const CPolynomial &b)
{
	int x=a.size()-b.size();
	if (x) return x<0;
	return a.compareTo(b) < 0;
}

/**\details
Подготавливает многочлены к использованию в качестве препроцессоров - сортирует и убирает одинаковые.
*/
void Unique(PolynomSet& polys)
{
	sort(polys.begin(), polys.end(), &cmpForPreprocess);
	polys.resize(unique(polys.begin(), polys.end())-polys.begin());
}

/**реализация препроцессинга (см теоретическую документацию)
\param polys множество многочленов, которое будет редуцироваться.
В результате препроцессинга к нему добавляются новые многочлены.
\param reducers множество возможных многочленов-препроцессоров.
Должно быть предварительно осортировано в соответствии с критерием оптимальности для использования в препроцессинге.
Более оптимаотные препроцессорв должнв стоять в начале.
*/
void Preprocess (PolynomSet& polys, PolynomSet& reducers)
{
	//MEASURE_TIME_IN_BLOCK("Preprocess");
	MonomialMap processed;
	MonomialMap monsToProcess;
	{
		//MEASURE_TIME_IN_BLOCK("MonomialsOf");
		MonomialsOf (monsToProcess, polys);	
	}
	
	CMonomial mon, mulby, HMR;	

	{
		//MEASURE_TIME_IN_BLOCK("Unique");
		Unique(reducers);
	}

	while(!monsToProcess.empty())
	{			
		mon = monsToProcess.selectMonoimial();
		monsToProcess.erase(mon);
				
		processed.storeMonomial(mon);

		for(int i = 0; i<reducers.size(); i++)
		{			
			HMR = reducers[i].HM();
			if(mon.tryDivide(HMR, mulby))
			{
				polys.push_back(reducers[i]);												
				polys.back()*=mulby;				
				{
					//MEASURE_TIME_IN_BLOCK("storeNotProcessedCMonomials");
					storeNotProcessedMonomialsFromPoly(monsToProcess, processed, polys.back());
				}
				break;
			}
			
		}		
	}		
}

///Оператор сравнения по числу ненулевых мономов
bool cmpForReduceBySize(const CPolynomial& a, const CPolynomial &b)
{
	int x=a.size()-b.size();
	if (x) return x<0;
	return a.compareTo(b) < 0;
}

///Оператор сравнения по упорядочиванию на многочленах
bool cmpForReduceByOrder(const CPolynomial& a, const CPolynomial &b)
{
	int x=a.compareTo(b);
	if (x) return x<0;
	return a.size()-b.size()>0;
}


/**авторедукция множества \a basis.
Проводит полную авторедукцию множества многочленов, используя те же матричные операции, что и Reduce().
(И, как следствие, эта функцтя тоже является распааллеленнной). Требует, чтоб аргумент \a basis представлял собой множество многочленов старшие мономы которых не делятся друг на друга (иначе резкльтат не определён). Используется для авторедукции полученного в итоге базиса Грёбнера для обеспечения его единственности.
*/
void AutoReduceBasis(PolynomSet& basis, const F4AlgData* f4options){	
	//MEASURE_TIME_IN_BLOCK("AutoReduce");
	PolynomSet minimalBasis;
	vector<bool> Mark;
	int Size = 0;
	CMonomial HMf;
	for(int i = 0; i<basis.size(); i++)
	{
		HMf = basis[i].HM();		
		for(int j = 0; j<minimalBasis.size(); j++)if(Mark[j]){
			if(minimalBasis[j].HM().divisibleBy(HMf)){
				Mark[j] = false, Size--;
			}
		}
		bool Condition = true;
		for(int j = 0; Condition && j<minimalBasis.size(); j++)if(Mark[j]){
			if(HMf.divisibleBy(minimalBasis[j].HM())){
				Condition = false;
			}
		}
		if(Condition)
		{
			minimalBasis.push_back(basis[i]);		
			Mark.push_back(true);
			Size++;
		}
	}
	PolynomSet reducers;
	reducers.reserve(Size);	
	MonomialMap minimalBasisMons;	
	for(int j = 0; j<minimalBasis.size(); j++)if(Mark[j])
	{
		reducers.push_back(minimalBasis[j]);
		minimalBasisMons.storeMonomial(minimalBasis[j].HM());
	}	
	minimalBasis.resize(reducers.size());
	for(int i = 0; i<reducers.size(); i++)
	{
		minimalBasis[i] = reducers[i];
	}
	
	Preprocess(minimalBasis, reducers);
	CMatrix matrix;
	polyToMatrix(minimalBasis,matrix);
	matrix.toDiagonalNormalForm(f4options);
	PolynomSet reducedPolys;
	MonomialMap empty;
	matrixToPoly(matrix,reducedPolys,empty);
	basis.clear();
	for(PolynomSet::iterator i = reducedPolys.begin(); i!=reducedPolys.end(); ++i)
	{
		if(minimalBasisMons.containsMonomial(i->HM()))
		{
			basis.push_back(*i);
		}
	}
}

/**функция сравнения для подготовки к Update.
Результат вызова Update для множества многочленов может существенно зависеть от того,
в каком порядке к этим многочленам Update применяется.
Сейчас в качестве критерия взято число ненулевых элементов.
*/
bool cmpForUpdaters(const CPolynomial& a, const CPolynomial &b){
	//Полное упорядочивание для гарантии одинакового результата вне зависимости от порядка строк в матрице
	int x=a.size()-b.size();
	if (x) return x<0;
	return a.compareTo(b)>0;
}

void SortByOrder(PolynomSet& set)
{
	//MEASURE_TIME_IN_BLOCK("SortByOrder");
	sort(set.begin(), set.end(), cmpForReduceByOrder);
}
void AutoReduceSetWithBasisFull(const PolynomSet& set, const F4AlgData* f4options, ReduceBySet& result)
{
	//MEASURE_TIME_IN_BLOCK("AutoReduceSetWithBasisFull");
	auto mySet = set;
	result = ReduceBySet();
	SortByOrder(mySet);
	sort(mySet.begin(), mySet.end(), cmpForReduceByOrder);
	for(auto greaterIt = mySet.cbegin(); greaterIt != mySet.cend(); ++greaterIt)
	{
		bool needAdd = true;
		const auto& resultSet = result.allReducers();
		for(auto smallerIt = resultSet.cbegin(); smallerIt != resultSet.cend(); ++smallerIt)
		{
			if (greaterIt->HM().divisibleBy(smallerIt->HM()))
			{
				needAdd = false;
				break;
			}
		}
		if (needAdd)
		{
			auto polyReduced = *greaterIt;
			result.reduceWithSearchCache(polyReduced);
			result.appendReducerGreater(polyReduced);
		}
	}
}
}
