#include "f4main.h"
#include "commonpolyops.h"
#include "outputroutines.h"
#include "conversions.h"
#include "reducebyset.h"
#include <iostream>
#include <list>
#include <cassert>
#include <map>
using namespace std;
namespace F4MPI
{

namespace PerryOrig
{

typedef CMonomialForLabel Label;

struct LabeledPolynom
{
	CPolynomial poly;
	Label label;
	void printTo(ostream& s)const
	{
		//poly.printPolynomial(s);
		s << poly.HM().toString();
		s << "@[" << label.mon().toString() << "]";
	}
};

typedef LabeledPolynom* LPolyPtr;
struct SPairPart
{
	bool isNew;
	union
	{
		struct
		{
			int polyIdx;
		} old;
		struct
		{
			LPolyPtr polyPtr;
		} n;
	} p;
	CMonomial multiplier;
};

struct SP;
struct F5Data;
struct BySigComparator
{
	BySigComparator(const F5Data& aData):data(aData){}
	const F5Data& data;
	bool operator()(LPolyPtr polyIdx1, LPolyPtr polyIdx2)const;
};



//S-пары должны быть отсортированы таким оьразом, чтоб можно было брать набор минимальных по степени
typedef multimap<int, SP> SortedSPairs;

//S-полиномы должны быть отсортированы таким оьразом, чтоб можно было брать минимальный по сигнатуре
typedef set<LPolyPtr, BySigComparator> SortedSPols;

//typedef map<int, int> SortedSPols;


struct F5Data
{
	F5Data(ReduceBySet& basis, const CPolynomial& extraCurrent):reduceHelper(basis)
	{
		const auto& allPolys = reduceHelper.allReducers();
		prevBasis.assign(allPolys.cbegin(), allPolys.cend());
		rulesOld.resize(prevBasis.size());

		//заполнение правил для существующих многочленов базиса
		for(int i = prevBasis.size() - 1; i>=0; --i)
		{
			//cout << "rules old for ";
			//basis[i].printPolynomial(cout);
			//cout << endl;
			CMonomial iTop = prevBasis[i].HM();
			for(int k = 0; k < i; ++k)
			{
				CMonomial kTop = prevBasis[k].HM();
				CMonomial kiLCM;
				CMonomial::lcm(iTop, kTop).tryDivide(iTop, kiLCM);
				rulesOld[i].push_back(kiLCM);
				//cout << kiLCM.toString() << " ";

			}
			//cout << " --" << endl;
		}
		GetBasisTops(prevBasis, prevBasisTops);
		//метка для первоначального многочлена - текущий индекс с единичным мономом
		gCurr.push_back(addPolyToEnd(extraCurrent));
	}
	typedef LPolyPtr RuleItem;
	
	bool checkCurReducer(const CMonomial& monToMulLabel, LPolyPtr poly)const
	{
		auto mon = monToMulLabel * poly->label.mon(); 
		if (CheckMonIsDivisibleBySome(mon, prevBasisTops))
		{
			//cout << "divisible by prev tops\n";
			//сигнатура домноженного делится на одну из сигнатур предыдущего базиса
			return false;
		}
		//правило отбрасывает редукцию, если нашлось правило, более позднее, чем соответствующее текущему полиному
		for(auto i = rulesCurrent.rbegin(); i != rulesCurrent.rend(); ++i)
		{
			if (poly == *i)
			{
				return true;
			}
			if (mon.divisibleBy((*i)->label.mon()))
			{
				//cout << "applying for idx "<<idxInPolys<<" rule with idx "<<i->idxInPolys <<" and mon "<<i->monomial.mon().toString()<<endl;
				return false;
			}
		}
		return true;
	}
	bool checkSPairPart(const SPairPart& part)const
	{
		if (part.isNew) return checkCurReducer(part.multiplier, part.p.n.polyPtr);
		else return checkOldReducer(part.p.old.polyIdx, part.multiplier);
	}
	const CPolynomial& polyForPart(const SPairPart& part)const
	{
		if (part.isNew) return part.p.n.polyPtr->poly;
		else return prevBasis[part.p.old.polyIdx];
	}
	bool reducePolynomiallyTopOnceByPrevBasis(CPolynomial& poly)
	{
		//MEASURE_TIME_IN_BLOCK("reducePolynomiallyTopOnceByPrevBasis");
		return reduceHelper.reduceTopOnceWithCache(poly);
		//return ReducePolynomiallyTopOnce(poly, prevBasis);
	}
	void reducePolynomiallyFullByPrevBasis(CPolynomial& poly)
	{
		//MEASURE_TIME_IN_BLOCK("reducePolynomiallyFullByPrevBasis");
		reduceHelper.reduceWithSearchCache(poly);
		//ReducePolynomiallyFull(poly, prevBasis);
	}
	void addNewPolyToEndAndRule(const Label& label, const CPolynomial& poly, SortedSPols& sPols)
	{
		//MEASURE_TIME_IN_BLOCK("addNewPolyToEndAndRule");
		polys.push_back(LabeledPolynom());
		auto& back = polys.back();
		back.label = label;
		back.poly = poly;
		//reducePolynomiallyFullByPrevBasis(back.poly);
		addRuleCurrent(&back);
		sPols.insert(&back);
	}
	LPolyPtr addPolyToEnd(const CPolynomial& poly)
	{
		//MEASURE_TIME_IN_BLOCK("addPolyToEnd");
		polys.push_back(LabeledPolynom());
		auto& back = polys.back();
		back.poly = poly;
		reducePolynomiallyFullByPrevBasis(back.poly);
		return &back;
	}
	const Label labelForPart(const SPairPart& part)const
	{
		Label result;
		if (part.isNew)
		{
			result = part.p.n.polyPtr->label * part.multiplier;
		}
		else
		{
			assert(!"labelForPart for isNew = false");
		}
		return result;
	}
	void addRuleCurrent(RuleItem item)
	{
		rulesCurrent.push_back(item);
	}
	PolynomSet prevBasis;
	list<LabeledPolynom> polys;
	vector<LPolyPtr> gCurr;
	ReduceBySet& reduceHelper;
private:
	vector<RuleItem> rulesCurrent;
	vector<vector<CMonomial>> rulesOld;
	vector<CMonomial> prevBasisTops;
	bool checkOldReducer(int polyInBasis, const CMonomial& mon)const
	{
		//cout << "checkOldReducer running for  " << mon.toString() << " polyIdx = " << polyInBasis << " poly = ";
		//prevBasis[polyInBasis].printPolynomial(cout);
		//cout<<endl;
		bool result = !CheckMonIsDivisibleBySome(mon, rulesOld[polyInBasis]);
		//cout<< "checkOldReducer result " << result <<endl;

		return result;
	}
};

bool BySigComparator::operator()(LPolyPtr polyIdx1, LPolyPtr polyIdx2)const
{
	return polyIdx1->label.compareToByLabelOrder(polyIdx2->label)<0;
}


struct SP
{
	CMonomial m;
	SPairPart s[2];
	static SP fromNewOld(LPolyPtr idx1, int idx2, const F5Data& data)
	{
		SP result;
		result.s[0].isNew = true;
		result.s[0].p.n.polyPtr = idx1;
		result.s[1].isNew = false;
		result.s[1].p.old.polyIdx = idx2;
		result.fillPair(data);
		return result;
	}
	static SP from2New(LPolyPtr idx1, LPolyPtr idx2, const F5Data& data)
	{
		SP result;
		result.s[0].isNew = true;
		result.s[0].p.n.polyPtr = idx1;
		result.s[1].isNew = true;
		result.s[1].p.n.polyPtr = idx2;
		result.fillPair(data);
		return result;
	}
	bool isUseFul(const F5Data& data)const
	{
		return data.checkSPairPart(s[0]) && data.checkSPairPart(s[1]);
	}
	void addSpolyIfUseful(F5Data& data, SortedSPols& sPols)const
	{
		//cout << "Analizing S-pair s0={" << s[0].isNew << ' ' << s[0].polyIdx << ' ' << s[0].multiplier.toString() << "}   s1={" << s[1].isNew << ' ' << s[1].polyIdx << ' ' << s[1].multiplier.toString() << "}" << endl;
		//cout << "thinking about spoly..." <<endl;
		bool skip = !isUseFul(data);
		if (skip) return;
		//добавление одного элемента к data.polys
		auto label = data.labelForPart(s[0]);
		if (s[1].isNew)
		{
			auto monlabel2 = data.labelForPart(s[1]);
			if (monlabel2.compareToByLabelOrder(label) > 0) label = monlabel2;
		}
		auto newPoly = data.polyForPart(s[0]);
		auto p1 = data.polyForPart(s[1]);
		p1 *= s[1].multiplier;
		p1 *= newPoly.HC()*CModular::inverseMod(p1.HC());
		newPoly *= s[0].multiplier;
		newPoly = newPoly - p1;
		data.addNewPolyToEndAndRule(label, newPoly, sPols);
		
		/*
		if (skip)
		{
			cout << "skipping s-poly";
			newPoly.poly.printPolynomial(cout);
			cout << endl;
			data.polys.resize(newPolyPos);
			return;
		}
		*/
		
		//cout << "added rule in addSpolyIfUseful polypos = "<<newPolyPos<<" mon = " <<newPoly.label.monomial.mon().toString() <<endl;
		//data.addRuleCurrent(label.monomial, newPolyPos);

		//для корректности сортировки в spos необходимо сначала заполнить метку в data.polys
		//sPols.insert(newPolyPos);
		//cout << "added S-Poly from it: ";
		//newPoly.poly.printPolynomial(cout);
		//cout << endl;
	}
	int deg()const
	{
		return m.getDegree();
	}
private:
	void fillPair(const F5Data& data)
	{
		const auto& p0 = data.polyForPart(s[0]);
		const auto& p1 = data.polyForPart(s[1]);
		const auto& h0 = p0.HM();
		const auto& h1 = p1.HM();
		m = CMonomial::lcm(h0, h1);
		m.tryDivide(h0, s[0].multiplier);
		m.tryDivide(h1, s[1].multiplier);
	}
	SP(){}
};

void AddSPtoSorted(const SP& sPair, SortedSPairs& pairs)
{
	//cout<<"before insert "<<pairs.size()<<endl;
	pairs.insert(make_pair(sPair.deg(), sPair));
	//cout<<"after insert "<<pairs.size()<<endl;
}

const LabeledPolynom* F5PerryFindReductor(const CInternalMonomial& monomial, const Label& sig, const F5Data& data)
{
	//MEASURE_TIME_IN_BLOCK("F5PerryFindReductor");
	SPairPart spp;
	spp.isNew = true;
	for(auto j = data.gCurr.crbegin(); j != data.gCurr.crend(); ++j)
	{
		auto& reductor = **j;
		const auto& hm = reductor.poly.HM();
		CMonomial u;
		if (monomial.tryDivide(hm, u))
		{
			spp.p.n.polyPtr = &reductor;
			spp.multiplier = u;
			if (0 != sig.compareToByLabelOrder(reductor.label * u) && data.checkSPairPart(spp)) return &reductor;
			//if (0 < sig.compareToByLabelOrder(reductor.label * u) && data.checkSPairPart(spp)) return &reductor;
		}
	}
	return 0;
}

void F5PerryTopReduction(LPolyPtr k, F5Data& data, SortedSPols& sPols)
{
	//MEASURE_TIME_IN_BLOCK("F5PerryTopReduction");
	auto& labeledPoly = *k;
	//cout << "analizing...";
	//labeledPoly.poly.printPolynomial(cout);
	//cout << endl;

	const LabeledPolynom* reductor;
	do
	{
		if (labeledPoly.poly.empty()) return;
		reductor = F5PerryFindReductor(labeledPoly.poly.HM(), labeledPoly.label, data);		
	}while(reductor == 0 && data.reducePolynomiallyTopOnceByPrevBasis(labeledPoly.poly));
	if (reductor == 0)
	{
		data.reducePolynomiallyFullByPrevBasis(labeledPoly.poly);
		//cout << "cant be reduced\n";
		labeledPoly.poly.normalize();
		//k нельзя отредуцировать, он добавляется к готовым
		data.gCurr.push_back(k);
		return;
	}
	auto reductorPoly = reductor->poly;
	CMonomial u;
	labeledPoly.poly.HM().tryDivide(reductorPoly.HM(), u);
	reductorPoly *= u;
	reductorPoly *= labeledPoly.poly.HC() * CModular::inverseMod(reductorPoly.HC());
	data.reducePolynomiallyFullByPrevBasis(reductorPoly);
	const auto reductorMonomialLabel = reductor->label * u;
	auto reducedPoly = labeledPoly.poly - reductorPoly;
	int compareResult = labeledPoly.label.compareToByLabelOrder(reductorMonomialLabel);
	//cout << "Comparing sigs to reduce " << labeledPoly.label.monomial.mon().toString() << " and reductor's " << reductorMonomialLabel.mon().toString() << "  result = " << compareResult << endl;
	if (compareResult > 0)
	{
		//разрешается обычная редукция, поскольку сигнатура рдуцируемого больше сигнатуры редуктора
		labeledPoly.poly = reducedPoly;
		//cout << "cnanged S-Poly in F5PerryTopReduction ";
		//reducedPoly.printPolynomial(cout);
		//cout << endl;

		sPols.insert(k);
		return;
	}
	data.addNewPolyToEndAndRule(reductorMonomialLabel, reducedPoly, sPols);
	//cout << "added rule in F5PerryTopReduction polypos = "<<newPolyPos<<" mon = " <<reductorMonomialLabel.mon().toString() <<endl;
	//data.addRuleCurrent(reductorMonomialLabel, newPolyPos);
	sPols.insert(k);	
	//cout << "added S-Poly in F5PerryTopReduction ";
	//newPoly.poly.printPolynomial(cout);
	//cout << endl;

	//sPols.insert(newPolyPos);
}

void F5PerryReduce(F5Data& data, SortedSPols& sPols)//может менять второй аргумент, добавляя S-полиномы в процессе работы
{
	//MEASURE_TIME_IN_BLOCK("F5PerryReduce");
	while(!sPols.empty())
	{
		auto curPoly = *sPols.begin();
		sPols.erase(curPoly);
		//полином редуцируется по старому базису по ссылке, непосредственно там где он лежит
		//cout << "before reducing by prev basis before F5PerryTopReduction ";
		//data.polys[curIdx].poly.printPolynomial(cout);
		//cout << endl;
		//data.reducePolynomiallyFullByPrevBasis(curPoly->poly);
		//все многочлены в Spols уже редуцированы относительно старого базиса
		if (curPoly->poly.empty())
		{
			cout << "reduction to zero before F5PerryTopReduction" <<endl;
			continue;
		}
		//cout << "after reducing by prev basis before F5PerryTopReduction ";
		//data.polys[curIdx].poly.printPolynomial(cout);
		F5PerryTopReduction(curPoly, data, sPols);
	}
}

void AddFromNewOldForNewPoly(const F5Data& data, LPolyPtr newPoly, SortedSPairs& sPairs)
{
	for(int i=0;i<int(data.prevBasis.size());++i)
	{
		AddSPtoSorted(SP::fromNewOld(newPoly, i, data), sPairs);
	}
}

void F5ByPerryStep(ReduceBySet& oldBasis, PolynomSet &basis, const CPolynomial& extraPolynom, const F4AlgData* /*f4options*/)
{
	//MEASURE_TIME_IN_BLOCK("F5ByPerryStep");
	F5Data data(oldBasis, extraPolynom);
	SortedSPairs sPairsToProcess;
	AddFromNewOldForNewPoly(data, &data.polys.front(), sPairsToProcess);
	auto sPols = SortedSPols(BySigComparator(data));
	while (!sPairsToProcess.empty())
	{
		const auto smallestKey = sPairsToProcess.cbegin()->first;
		for (auto sPairIt = sPairsToProcess.cbegin(); sPairIt != sPairsToProcess.cend() && sPairIt->first == smallestKey; ++sPairIt)
		{
			const auto& sPair = sPairIt->second;
			sPair.addSpolyIfUseful(data, sPols);
		}
		sPairsToProcess.erase(smallestKey);
		int prevGCurrSize = data.gCurr.size();
		F5PerryReduce(data, sPols);
		for (int k = prevGCurrSize; k < int(data.gCurr.size()); ++k)
		{
			AddFromNewOldForNewPoly(data, data.gCurr[k], sPairsToProcess);
			for (int j = 0; j < k; ++j)
			{
				//cout<<"new "<<j<<" "<<k<<endl;
				AddSPtoSorted(SP::from2New(data.gCurr[k], data.gCurr[j], data), sPairsToProcess);
			}
		}
		cout << "Computed basis till degree " << smallestKey << " gcurr size = " << data.gCurr.size() << endl;

		/*
		for(auto j = data.gCurr.cbegin(); j != data.gCurr.cend(); ++j)
		{
			(*j)->printTo(cout);
			auto k = j;
			for(++k; k != data.gCurr.cend(); ++k)
			{
				auto& pj = **j;
				auto& pk = **k;
				if (true || pj.poly.HM().divisibleBy(pk.poly.HM()))
				{
					pj.printTo(cout);
					cout << " is  divisible by ";
					pk.printTo(cout);
					cout << endl;
				}
			}
		}
		*/
		//for (int i = lastGCurrSize; i < data.gCurr.size(); ++i)
		//{
			//const auto& newp = data.polys[data.gCurr[i]];
			//cout<<"HM="<<newp.poly.HM().toString()<<" label="<<newp.label.monomial.mon().toString()<< " idx in polys= "<<data.gCurr[i]<<endl;
		//}
		//cout << "Spairs: " << sPairsToProcess.size() << endl;
		//data.reduceHelper.freeCache();
	}
	basis = data.prevBasis;
	for(auto newInBasis = data.gCurr.cbegin(); newInBasis != data.gCurr.cend(); ++newInBasis)
	{
		basis.push_back((*newInBasis)->poly);
	}
}

PolynomSet F5ByPerry(PolynomSet F, const F4AlgData* f4options){
	//MEASURE_TIME_IN_BLOCK("F5ByPerry");
	Normalize(F);
	ReduceBySet reducedBasis;
	for(auto i = F.begin(); i != F.end(); ++i)
	{
		PolynomSet newBasis;
		F5ByPerryStep(reducedBasis, newBasis, *i, f4options);
		//cout << "size before reduction " << basis.size() << endl;
		//PrintPolynomSetRaw(cout, basis);
		//PrintPolynomSetRaw(cout, basis);
		AutoReduceSetWithBasisFull(newBasis, f4options, reducedBasis);
		cout << "size after reduction " << reducedBasis.allReducers().size() << endl;
		//PrintPolynomSetRaw(cout, basis);
	}
	const auto& result = reducedBasis.allReducers();
	return PolynomSet(result.begin(), result.end());
}
}
}
