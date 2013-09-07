#include "f4main.h"
#include "commonpolyops.h"
#include "outputroutines.h"
#include "conversions.h"
#include "reducebyset.h"
#include <vector>
#include <iostream>
#include <deque>
#include <cassert>
#include <array>
#include <map>


static const bool more_output = false;
using namespace std;
namespace F4MPI
{

namespace F5Orig
{

struct Label : CMonomialForLabel
{
	int index;
};

struct LabeledPolynom
{
	CPolynomial poly;
	Label label;
	void printTo(ostream& s)const
	{
		poly.printPolynomial(s);
		//s << poly.HM().toString();
		s << "@["<< label.index << "." << label.mon().toString() << "]";
	}
};

typedef LabeledPolynom* LPolyPtr;

std::string cutbraces(const char *multilinestr)
{
	return std::string(multilinestr +1, strlen(multilinestr) - 2);
}
#define MULTILINESTR(x) (cutbraces(#x))

static const char* const LatexNewLine = "\\\\*\n";
std::string LatexVarName(int idx)
{
	const std::array<const char*, 5> first_vars = {"x","y","z","t","w"};
	if (idx < first_vars.size()) return first_vars[idx];
	ostringstream s;
	s << "v_{" << idx << "}";
	return s.str();
}
std::string LatexOrderDescr()
{
	ostringstream s;
	s << "\\mbox{" << getMonomialOrderName(CMonomial::getOrder()) << "}(";
	for (int i=0;i<CMonomial::theNumberOfVariables;++i)
	{
		if (i>0) s << ">";
		s << LatexVarName(i);
	}
	s << ")";
	return s.str();
}
std::string Latex(const CMonomial& mon)
{
	if (mon.isOne()) return "\\monone";
	ostringstream s;
	for (int i = 0; i < CMonomial::theNumberOfVariables; ++i)
	{
		if (int d = mon.getDegree(i))
		{
			s << LatexVarName(i);
			if (d > 1)
			{
				s << "^{" << d << "}";
			}
		}
	}
	return s.str();
}

std::string Latex(const CPolynomial& p)
{
	if (p.empty())
	{
		return "0";
	}
	ostringstream s;
	int i = 0;
	for(;;){
		int coeff=p.getCoeff(i).toint();
		assert(coeff != 0);
		if (coeff != 1)
		{
			s << coeff;
		}
		s << Latex(p.getMon(i));
		++i;
		if(i==p.size())
		{
			break;
		}
		s<<"+";
	}
	return s.str();
}

std::string Latex(LPolyPtr lpoly)
{
	ostringstream s;
	s << "\\lpoly{" << Latex(lpoly->label.mon()) <<"}{"<< lpoly->label.index << "}{" << Latex(lpoly->poly) <<"}";
	return s.str();
}

std::string Latex(Label label)
{
	ostringstream s;
	s << "\\Sss{" << Latex(label.mon()) <<"}{"<< label.index << "}";
	return s.str();
}

typedef std::vector<LPolyPtr> LabeledPolyVector;
struct SPairPart
{
	LPolyPtr polyPtr;
	CMonomial multiplier;
};

struct SP;
struct BySigComparator
{
	bool operator()(LPolyPtr polyIdx1, LPolyPtr polyIdx2)const
	{
		assert(polyIdx1->label.index == polyIdx2->label.index);
		return polyIdx1->label.compareToByLabelOrder(polyIdx2->label)<0;
	}
};


//S-пары должны быть отсортированы таким оьразом, чтоб можно было брать набор минимальных по степени
typedef multimap<int, SP> SortedSPairs;

//S-полиномы должны быть отсортированы таким оьразом, чтоб можно было брать минимальный по сигнатуре
typedef set<LPolyPtr, BySigComparator> SortedSPols;

struct F5Globals
{
	deque<LabeledPolynom> R;
	vector<deque<LPolyPtr>> Rule;
	int m() //number of polynomials in input
	{
		return Rule.size();
	}
	ostream* latex_log;
};

int PosForLogging(LPolyPtr lpoly, const F5Globals& f5_globals)
{
	auto i = f5_globals.R.begin();
	while (&*i != lpoly) ++i;
	return 1 + (i - f5_globals.R.begin());
}


std::string LatexRI(LPolyPtr lpoly, const F5Globals& f5_globals)
{
	ostringstream s;
	s << "R["<< PosForLogging(lpoly, f5_globals)<<"]";
	return s.str();
}


bool RejectedByF5Criteria(const CMonomial& monToMulLabel, LPolyPtr poly, ReduceBySet& phi, const F5Globals& f5_globals)
{
	//polynomial is rejected if its signature is reducible by phi
	CPolynomial poly_to_test;
	poly_to_test.pushTermBack(1, monToMulLabel * poly->label.mon());
	CMonomial mulby_unused;
	auto result = phi.findTopReducer(poly_to_test, mulby_unused);
	if (result && f5_globals.latex_log)
	{
		*f5_globals.latex_log << " отбрасывается $\\varphi$-проверкой, ибо моном сигнатуры $"<<Latex(monToMulLabel)<<"\\cdot"<<Latex(poly->label)
		<<"$ делится на старший моном $"<<Latex(result->HM())<<"$ базиса прошлого шага.";
	}
	
	return bool(result);
}

bool IsRewritten(const CMonomial& monToMulLabel, LPolyPtr poly, const F5Globals& f5_globals)
{
	auto mon = monToMulLabel * poly->label.mon();
	//правило отбрасывает редукцию, если нашлось правило, более позднее, чем соответствующее текущему полиному
	const auto& rulesCurrent = f5_globals.Rule[poly->label.index - 1];
	for(auto i = rulesCurrent.rbegin(); i != rulesCurrent.rend(); ++i)
	{
		if (poly == *i)
		{
			return false;
		}
		if (mon.divisibleBy((*i)->label.mon()))
		{
			if (f5_globals.latex_log)
			{
				*f5_globals.latex_log << "отбрасывается критерием перезаписи, поскольку сигнатура $"<<Latex(monToMulLabel)<<"\\cdot"<<Latex(poly->label)
				<<"$ делится на сигнатуру элемента $"<<LatexRI(*i, f5_globals)<<"$.";
			}
			return true;
		}
	}
	return false;
}

struct SP
{
	CMonomial m;
	SPairPart s[2];

	int deg()const
	{
		return m.getDegree();
	}
};

void AddSPairs(const std::vector<SP>& sPairs, SortedSPairs& P)
{
	for(const SP& sPair:sPairs){
		//cout<<"before insert "<<pairs.size()<<endl;
		P.insert(make_pair(sPair.deg(), sPair));
		//cout<<"after insert "<<pairs.size()<<endl;
	}
}

LPolyPtr AddNewPolyAndRule(const Label& label, const CPolynomial& poly, F5Globals& f5_globals)
{
	//MEASURE_TIME_IN_BLOCK("addNewPolyInToDoAndRule");
	f5_globals.R.push_back(LabeledPolynom());
	auto& back = f5_globals.R.back();
	back.label = label;
	back.poly = poly;
	//PhiNormalForm_ByPrevBasis(back.poly);
	auto& cur_rule = f5_globals.Rule[label.index - 1];
	cur_rule.push_back(&back);
	if (f5_globals.latex_log)
	{
		*f5_globals.latex_log << "добавляется как $"<<LatexRI(&back, f5_globals)<<"="<<Latex(&back)<<"$ и попадает в правила, делая Rule["<<label.index<<"] равным $[";
		bool first = true;
		for(auto i = cur_rule.rbegin(); i != cur_rule.rend(); ++i)
		{
			if (i != cur_rule.rbegin())
			{
				*f5_globals.latex_log << ", ";
			}
			*f5_globals.latex_log << PosForLogging(*i, f5_globals);
		}
		
		*f5_globals.latex_log << "]$.";
	}
	return &back;
}

LPolyPtr F5IsReducible(LPolyPtr poly, const LabeledPolyVector& GWithTodo, const F5Globals& f5_globals, ReduceBySet& phi)
{
	//MEASURE_TIME_IN_BLOCK("F5IsReducible");
	if (more_output) {
		cout << "Searching reductor for ";
		poly->printTo(cout);
		cout << endl;
	}
	const CInternalMonomial& monomial = poly->poly.HM();
	const Label& sig = poly->label;
	for(auto& reductor : GWithTodo)
	{
		const auto& hm = reductor->poly.HM();
		CMonomial mulby;
		if (!monomial.tryDivide(hm, mulby)) continue;
		if (f5_globals.latex_log)
		{
			*f5_globals.latex_log << LatexNewLine << "TopReduction: потенциальный редуктор $" << LatexRI(reductor, f5_globals) << "$ ";
		}
		assert(sig.index == reductor->label.index);
		const CMonomial u = mulby;
		if (more_output) {
			cout << " testing:";
			reductor->printTo(cout);
		}
		if (RejectedByF5Criteria(u, reductor, phi, f5_globals))
		{
			if (more_output) {
				cout << "rejected by f5 criteria. " << endl;
			}
			continue;
		}
		if (IsRewritten(u, reductor, f5_globals))
		{
			if (more_output) {
				cout << "rejected by rewritten criteria. " << endl;
			}
			continue;
		}
		if (0 == sig.compareToByLabelOrder(reductor->label * u))
		{
			if (more_output) {
				cout << "rejected by equal signatures. " << endl;
			}
			continue;
		}
		if (more_output) {
			cout << " found." << endl;
		}
		return reductor;
	}
	if (more_output) {
		cout << " not found\n";
	}
	return 0;
}


struct TopReductionResult
{
	LabeledPolyVector ToDo, Done;
};

TopReductionResult F5TopReduction(LPolyPtr r0, F5Globals& f5_globals, const LabeledPolyVector& GWithTodo, ReduceBySet& phi)
{
	TopReductionResult result;
	//MEASURE_TIME_IN_BLOCK("F5TopReduction");
	auto& labeledPoly = *r0;
	if (labeledPoly.poly.empty())
	{
		cout << "Zero polynomial got to F5TopReduction" <<endl;
		return result;
	}
	if (f5_globals.latex_log)
	{
		*f5_globals.latex_log << LatexNewLine << "TopReduction: поиск редукторов $" << LatexRI(&labeledPoly, f5_globals) <<"$.";
	}
	
	//cout << "analizing...";
	//labeledPoly.poly.printPolynomial(cout);
	//cout << endl;


	LPolyPtr reductor = F5IsReducible(&labeledPoly, GWithTodo, f5_globals, phi);

	if (reductor == 0)
	{
		if (f5_globals.latex_log)
		{
			*f5_globals.latex_log << LatexNewLine << "TopReduction: $" << LatexRI(&labeledPoly, f5_globals) << "$ не может быть далее редуцирован и добавляется в Done.";
		}
		
		//F5 does not full reduce data.PhiNormalForm_ByPrevBasis(labeledPoly.poly);
		//cout << "cant be reduced\n";
		//r0 нельзя отредуцировать, он добавляется к готовым
		result.Done.push_back(r0);
		return result;
	}
	auto reductorPoly = reductor->poly;
	CMonomial u;
	labeledPoly.poly.HM().tryDivide(reductorPoly.HM(), u);
	reductorPoly *= u;
	reductorPoly *= labeledPoly.poly.HC() * CModular::inverseMod(reductorPoly.HC());
	Label reductorMonomialLabel = reductor->label;
	reductorMonomialLabel *= u;
	auto reducedPoly = labeledPoly.poly - reductorPoly;
	if (!reducedPoly.empty()) reducedPoly.normalize();
	int compareResult = labeledPoly.label.compareToByLabelOrder(reductorMonomialLabel);
	//cout << "Comparing sigs to reduce " << labeledPoly.label.monomial.mon().toString() << " and reductor's " << reductorMonomialLabel.mon().toString() << "  result = " << compareResult << endl;
	if (compareResult > 0)
	{
		//разрешается обычная редукция, поскольку сигнатура рдуцируемого больше сигнатуры редуктора
		labeledPoly.poly = reducedPoly;
		if (f5_globals.latex_log)
		{
			*f5_globals.latex_log << "редуцирует $" << LatexRI(&labeledPoly, f5_globals) << "$ до $"<<Latex(&labeledPoly)<<"$ и результат остаётся в ToDo.";
		}
		//cout << "cnanged S-Poly in F5TopReduction ";
		//reducedPoly.printPolynomial(cout);
		//cout << endl;
	}
	else
	{
		if (f5_globals.latex_log)
		{
			*f5_globals.latex_log << "не может произвести редукцию с сохранением сигнатуры. S-многочлен $" << LatexRI(&labeledPoly, f5_globals) <<"$ и $"<< LatexRI(reductor, f5_globals) <<"$ ";
		}
		result.ToDo.push_back(AddNewPolyAndRule(reductorMonomialLabel, reducedPoly, f5_globals));
	}
	//cout << "added rule in F5TopReduction polypos = "<<newPolyPos<<" mon = " <<reductorMonomialLabel.mon().toString() <<endl;
	//data.addRuleCurrent(reductorMonomialLabel, newPolyPos);
	
	result.ToDo.push_back(r0);
	return result;
	//cout << "added S-Poly in F5TopReduction ";
	//newPoly.poly.printPolynomial(cout);
	//cout << endl;

	//sPolsToDo.insert(newPolyPos);
}

LabeledPolyVector F5Reduction(LabeledPolyVector &G, F5Globals& f5_globals, ReduceBySet& phi, SortedSPols& sPolsToDo)//может менять второй аргумент, добавляя S-полиномы в процессе работы
{
	//MEASURE_TIME_IN_BLOCK("F5Reduction");
	LabeledPolyVector GWithDone = G;
	LabeledPolyVector Done;
	while(!sPolsToDo.empty())
	{
		auto h = *sPolsToDo.begin();
		sPolsToDo.erase(h);
		//cout << "before reducing by prev basis before F5TopReduction ";
		//data.R[curIdx].poly.printPolynomial(cout);
		//cout << endl;
		//полином редуцируется по старому базису по ссылке, непосредственно там где он лежит
		phi.reduceTopWithCache(h->poly);
		//cout << "after reducing by prev basis before F5TopReduction ";
		//data.R[curIdx].poly.printPolynomial(cout);
		auto top_red_resut = F5TopReduction(h, f5_globals, GWithDone, phi);
		if (more_output) {
			cout << "F5TopReduction return size todo:" << top_red_resut.ToDo.size() << "  done:" << top_red_resut.Done.size() << endl;
		}
		sPolsToDo.insert(top_red_resut.ToDo.begin(), top_red_resut.ToDo.end());
		Done.insert(Done.end(), top_red_resut.Done.begin(), top_red_resut.Done.end());
		GWithDone.insert(GWithDone.end(), top_red_resut.Done.begin(), top_red_resut.Done.end());
	}
	return Done;
}

std::vector<SP> F5CritPairs(LPolyPtr r1, const LabeledPolyVector &r2_parts, int index_k, ReduceBySet& phi, F5Globals& f5_globals)
{
	std::vector<SP> result;
	for(auto r2 : r2_parts)
	{
		SP sp;
		sp.s[0].polyPtr = r1;
		sp.s[1].polyPtr = r2;
		const auto& h0 = sp.s[0].polyPtr->poly.HM();
		const auto& h1 = sp.s[1].polyPtr->poly.HM();
		sp.m = CMonomial::lcm(h0, h1);
		sp.m.tryDivide(h0, sp.s[0].multiplier);
		sp.m.tryDivide(h1, sp.s[1].multiplier);
		if (f5_globals.latex_log)
		{
			*f5_globals.latex_log << LatexNewLine << "CritPair: S-пара $("
			<<Latex(sp.s[0].multiplier)<<"\\cdot "<<LatexRI(sp.s[0].polyPtr, f5_globals)<<","
			<<Latex(sp.s[1].multiplier)<<"\\cdot "<<LatexRI(sp.s[1].polyPtr, f5_globals)<<")$ ";
		}
		
		if (RejectedByF5Criteria(sp.s[0].multiplier, sp.s[0].polyPtr, phi, f5_globals))
		{
			continue;
		}
		if (index_k == sp.s[1].polyPtr->label.index && RejectedByF5Criteria(sp.s[1].multiplier, sp.s[1].polyPtr, phi, f5_globals))
		{
			continue;
		}
		if (more_output) {
			cout << "SPAIR  " << sp.s[0].multiplier.toString() << " * ";
			sp.s[0].polyPtr->printTo(cout);
			cout << " and "  << sp.s[1].multiplier.toString() << " * ";
			sp.s[1].polyPtr->printTo(cout);
			cout << endl;
		}
		if (f5_globals.latex_log)
		{
			*f5_globals.latex_log << "добавляется в $P$.";
		}
		result.push_back(sp);
	}
	return result;
}

SortedSPols F5Spols(const vector<SP>& P, F5Globals& f5_globals)
{
	SortedSPols result;
	for(const SP& sp : P)
	{
		if (f5_globals.latex_log)
		{
			*f5_globals.latex_log << LatexNewLine << "Spol: S-многочлен $"
			<<Latex(sp.s[0].multiplier)<<"\\cdot "<<LatexRI(sp.s[0].polyPtr, f5_globals)<<" - "
			<<Latex(sp.s[1].multiplier)<<"\\cdot "<<LatexRI(sp.s[1].polyPtr, f5_globals)<<"$ ";
		}
		
		if (IsRewritten(sp.s[0].multiplier, sp.s[0].polyPtr, f5_globals) || IsRewritten(sp.s[1].multiplier, sp.s[1].polyPtr, f5_globals))
		{
			continue;
		}
		//добавление одного элемента к data.R
		auto label = sp.s[0].polyPtr->label;
		label *= sp.s[0].multiplier;
		auto monlabel2 = sp.s[1].polyPtr->label;
		monlabel2 *= sp.s[1].multiplier;
		if (monlabel2.index == label.index)
		{
			if (monlabel2.compareToByLabelOrder(label) > 0) label = monlabel2;
		}
		auto newPoly = sp.s[0].polyPtr->poly;
		auto p1 = sp.s[1].polyPtr->poly;
		p1 *= sp.s[1].multiplier;
		p1 *= newPoly.HC()*CModular::inverseMod(p1.HC());
		newPoly *= sp.s[0].multiplier;
		newPoly = newPoly - p1;

		if (newPoly.empty())
		{
			cerr << "REDUCTION TO ZERO" << endl;
			continue;
		}
		
		newPoly.normalize();
		auto lpoly_ptr = AddNewPolyAndRule(label, newPoly, f5_globals);
		result.insert(lpoly_ptr);
		if (more_output) {
			cout << "computed spol  " << sp.s[0].multiplier.toString() << " * ";
			sp.s[0].polyPtr->printTo(cout);
			cout << " - "  << sp.s[1].multiplier.toString() << " * ";
			sp.s[1].polyPtr->printTo(cout);
			cout << " = ";
			lpoly_ptr->printTo(cout);
			cout << endl;
		}
	}
	cout << "SPOLS return size " << result.size() << endl;
	return result;
}

LabeledPolyVector AlgorithmF5(const LabeledPolyVector &Gprev, const LPolyPtr ith_poly, F5Globals& f5_globals)
{
	//MEASURE_TIME_IN_BLOCK("AlgorithmF5");
	LabeledPolyVector G = Gprev;
	G.push_back(ith_poly);
	PolynomSet oldBasisAsPlainSet;
	for (auto labeled_poly : Gprev)
	{
		oldBasisAsPlainSet.push_back(labeled_poly->poly);
	}
	ReduceBySet phi(oldBasisAsPlainSet);
	int index_i = ith_poly->label.index;
	cout << "Processing index " << index_i << endl;
	if (f5_globals.latex_log)
	{
		*f5_globals.latex_log << LatexNewLine << "\\textbf{Шаг "<<f5_globals.m() - index_i<<":} AlgorithmF5 обрабатывает многочлен $"<< LatexRI(ith_poly, f5_globals)
		<<"$ и формирует по предыдущему базису $G_{"<<index_i<<"}=G_"<<index_i+1<<"\\cup \\{" << PosForLogging(ith_poly, f5_globals) << "\\}$.";
	}		
	SortedSPairs P;
	AddSPairs(F5CritPairs(ith_poly, Gprev, index_i, phi, f5_globals), P);
	while (!P.empty())
	{
		const auto d = P.cbegin()->first;//smallest deg
		if (f5_globals.latex_log)
		{
			*f5_globals.latex_log << LatexNewLine << "AlgorithmF5: выбор S-пар степени "<< d << ".";
		}
		vector<SP> Pd;
		for (auto sPairIt = P.cbegin(); sPairIt != P.cend() && sPairIt->first == d; ++sPairIt)
		{
			Pd.push_back(sPairIt->second);
		}
		P.erase(d);
		SortedSPols Fd = F5Spols(Pd, f5_globals);
				
		auto Rd = F5Reduction(G, f5_globals, phi, Fd);
		//iterate over new polynomials added to G
		for (auto& r : Rd)
		{
			AddSPairs(F5CritPairs(r, G, index_i, phi, f5_globals), P);
			G.push_back(r);
		}
		if (f5_globals.latex_log && !Rd.empty())
		{
			*f5_globals.latex_log << LatexNewLine << "AlgorithmF5: добавляет в $G_{"<<index_i<<"}$ следующее множество позиций многочленов в массиве $R: \\{";
			for(auto i = Rd.begin(); i != Rd.end(); ++i)
			{
				if (i != Rd.begin())
				{
					*f5_globals.latex_log << ", ";
				}
				*f5_globals.latex_log << PosForLogging(*i, f5_globals);
			}
			
			*f5_globals.latex_log << "\\}$.";
		}
		cout << "Computed basis till degree " << d << " G size = " << G.size() << endl;
		
		if (more_output) {

			for(auto j = G.cbegin(); j != G.cend(); ++j)
			{
				(*j)->printTo(cout);
				cout << endl;
				auto k = j;
				for(++k; k != G.cend(); ++k)
				{
					auto& pj = **j;
					auto& pk = **k;
					if (pj.poly.HM().divisibleBy(pk.poly.HM()))
					{
						pj.printTo(cout);
						cout << " is  divisible by ";
						pk.printTo(cout);
						cout << endl;
					}
				}
			}
		}
	}
	if (f5_globals.latex_log)
	{
		*f5_globals.latex_log << LatexNewLine << "Результирующим базисом шага оказываются многочлены R с позициями из $G_"<<index_i<<"=\\{";
		for(auto i = G.begin(); i != G.end(); ++i)
		{
			if (i != G.begin())
			{
				*f5_globals.latex_log << ", ";
			}
			*f5_globals.latex_log << PosForLogging(*i, f5_globals);
		}
		*f5_globals.latex_log << "\\}$.";
	}
	return G;
}

PolynomSet IncrementalF5(PolynomSet F, const F4AlgData* f4options){
	//MEASURE_TIME_IN_BLOCK("IncrementalF5");
	Normalize(F);
	F5Globals f5_globals;
	f5_globals.latex_log = f4options->latex_log;
	LabeledPolynom f_i;
	f_i.label.index = 1;
	for(auto poly = F.rbegin(); poly != F.rend(); ++poly, ++f_i.label.index)
	{
		f_i.poly = *poly;
		f5_globals.R.push_back(f_i);
	}
	f5_globals.Rule.resize(f5_globals.R.size());
	if (f5_globals.latex_log)
	{
		*f5_globals.latex_log << "Мы будем рассматривать работу алгоритма в кольце многочленов над конечным полем $\\mathbb{Z}_{"<<CModular::getMOD()<<"}$, "
		<< "с мономиальным порядком $" << LatexOrderDescr() << "$." << LatexNewLine
		<< MULTILINESTR((
\\textbf{Инициализация:} IncrementalF5 записывает в массив $R$ входные многочлены, добавляя к ним сигнатуры с единичным
мономом и индексом, равным порядковому номеру во входных данных:\\\\*
$R\\leftarrow[
								 ));
		bool first = true;
		for(int i = 0; i < f5_globals.R.size(); ++i)
		{
			auto lpoly = f5_globals.R[i];
			*f5_globals.latex_log << "["<<i+1<<"]="<<Latex(&lpoly);
			if (i != f5_globals.R.size() - 1)
			{
				*f5_globals.latex_log << ",$" << LatexNewLine << "$";
			}
		}
		*f5_globals.latex_log << "]$" << LatexNewLine
		<< "$G_"<<f5_globals.m()<<"$ присваивается значение $\\{"<<f5_globals.m()<<"\\}$.";
	}
	auto ith_poly = f5_globals.R.rbegin();
	LabeledPolyVector G;
	G.push_back(&*ith_poly);
	++ith_poly;
	for(; ith_poly != f5_globals.R.rend(); ++ith_poly)
	{
		G = AlgorithmF5(G, &*ith_poly, f5_globals);
	}
	PolynomSet result;
	if (f5_globals.latex_log)
	{
		*f5_globals.latex_log << LatexNewLine <<"После убирания сигнатур финальный нередуцированный базис принимает вид:";
	}
	for (auto labeled_poly : G)
	{
		if (f5_globals.latex_log)
		{
			*f5_globals.latex_log << LatexNewLine << "$" <<LatexRI(labeled_poly, f5_globals)<<"="<<Latex(labeled_poly->poly)<< "$";
		}
		result.push_back(labeled_poly->poly);
	}
	if (f5_globals.latex_log)
	{
		*f5_globals.latex_log << ".";
	}
	ReduceBySet reducedBasis;
	AutoReduceSetWithBasisFull(result, f4options, reducedBasis);
	cout << "size after reduction " << reducedBasis.allReducers().size() << endl;
	return result;
}
}
}
