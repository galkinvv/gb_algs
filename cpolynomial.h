#ifndef CPolynomial_h
#define CPolynomial_h
/**
\file
Представление многочленов.
Определяет простую рализацию многочленов CPlainPolynomial и обёртку над ней с подсчётом ссылок CRefPolynomial.
В качестве CPolynomial может использоваться любая из них.
*/

#include "cmonomial.h"
#include "cmodular.h"
#include "algs.h"

#include <vector>
#include <cassert>
#include <algorithm>
#include <ostream>

namespace F4MPI{
///Представление полиномов над полем CModular.
class CPlainPolynomial:public IntrusiveRefCountBase{
public:
	typedef CInternalMonomial MonomialInPoly;
private:	
	///тип содержащий массив коэффициентов
	typedef std::vector<CModular> CoeffContanier;
	/**тип содержащий массив мономов.
	Поскольку размер, занимаемый мономом не известен заранее,
	а использование хранения данных мономов во внешней памяти в виде CMonomial неэффективно,
	применяется PODvecSize<MonomialInPoly> - массив элементов, размер которых не известен на момент компиляции.
	*/
	typedef PODvecSize<MonomialInPoly> MonomialContanier;
	/**коэффициенты при мономах.
	Нулевые коэффициенты в массиве отсутсвуют.
	*/
	CoeffContanier coeffs;
	/**мономы, присутсвтующие в многочлене
	мономы многочлена упорядочены по убыванию в соответствии с порядком,
	поэтому старший моном всегда находится в начале массива.
	*/
	MonomialContanier mons;
public:
	///итераторы по коээфициентам при мономах
	typedef CoeffContanier::iterator c_iterator;
	typedef CoeffContanier::const_iterator c_const_iterator;
	//итераторы мономам
	typedef MonomialContanier::iterator m_iterator;
	typedef MonomialContanier::const_iterator m_const_iterator;
	
	CPlainPolynomial(){}

	///устанавливает размер (массивов коэффициентов и мономов)
	void resize(size_t sz){
		coeffs.resize(sz);
		mons.resize(sz);
	}
	void reserve(size_t sz)
	{
		coeffs.reserve(sz);
		mons.reserve(sz);
	}
	
	void clear()
	{
		coeffs.clear();
		mons.clear();
	}
	///число мономов с ненулевыми коэффициентами
	size_t size()const{
		return coeffs.size();
	}

	///проверка, имеется ли хоть один терм
	bool empty() const {
		return size()==0;
	}
	
	/**Делает единичным старший коэффициент.
	Нормализует строку, т.е. домножает все элементы строки на коэффициент обратный старшему.
	Таким образом старший коэффициент становится единичным.
	*/
	void normalize(){
		(*this)*=CModular::inverseMod(HC());
	}

	m_const_iterator m_end()const{
		return mons.end();
	}

	m_const_iterator m_begin()const{
		return mons.begin();
	}

	///возвращает итератор по мономам - позиция за последним
	m_iterator m_end(){
		return mons.end();
	}

	///возвращает итератор по мономам - начальная позиция
	m_iterator m_begin(){
		return mons.begin();
	}

	///возвращает итератор по коэффициентам - позиция за последним
	c_iterator c_end(){
		return coeffs.end();
	}

	///возвращает итератор по коэффициентам - начальная позиция
	c_iterator c_begin(){
		return coeffs.begin();
	}

	const MonomialInPoly& getMon(int i)const{
		return mons[i];
	}

	///Возвращает ссылку на i-й моном
	MonomialInPoly& getMon(int i){
		return mons[i];
	}

	CModular getCoeff(int i)const{
		return coeffs[i];
	}

	///Возвращает ссылку на i-й коэффициент
	CModular& getCoeff(int i){
		return coeffs[i];
	}

	/**добавляет терм к полиному.
	Терм добавляется таким образом, чтоб сохранилось упорядочение в массиве мономов
	Время работы - O(size()) из-за того, что приходится сдвигать все элементы массивов коэффициентов и мономов. 
	*/
	void addTerm(CModular c, const CMonomial& m);

	/**добавляет терм в конец полинома.
	Для того чтоб при добавлении терма сохранилось свойство упорядоченнсти мономов,
	необходимо, чтоб моном \a m, переданный в качестве аргумента был младше всех мономов, содержащихся в полиноме.
	Если это требование нарушено, многочлен будет некорректен и дальнейшее поведение не определено.
	*/
	void pushTermBack(CModular c, const CMonomial& m);

	/**
	Текстовое представление полинома.
	Выводит мнгочлен в текстовос виде.
	\param output поток, в который нужно произвести вывод
	\param names соответствие между индексами и текстовыми именами переменных.
	Если передан нулевой указатель (по умолчанию), используются стандартные имена x1 .. xN
	*/
	void printPolynomial(std::ostream& output,ParserVarNames* names=0) const;

	///возвращает старший моном полинома
	const MonomialInPoly& HM() const{
		return getMon(0);
	}

	///возвращает коэффициент при старшем мономе
	CModular HC() const{
		return getCoeff(0);
	}

	/**отношение порядка на полиномах.
	Сравнивает полиномы, используя лексикографическое упорядочивание, основанное на порядке на мономах.
	\retval число, определяющее результат сравнения:
	\arg -1: этот полином младше \a p2
	\arg 1: этот полином старше \a p2
	\arg 0: этот полином равен \a p2
	*/
	int compareTo(const CPlainPolynomial& p2)const{
		int res;
		for (int i=0,iend=std::min(size(),p2.size());i!=iend;++i){
			res=getMon(i).compareTo(p2.getMon(i));
			if (res) return res;
			res=getCoeff(i).compareTo(p2.getCoeff(i));
			if (res) return res;
		}
		res=size()-p2.size();
		if (res<0) return 1;
		if (res>0) return -1;
		return 0;
	}

	/**меняет местами многочлены.
	Эффективная реализация swap работает за O(1) просто обменивает массивы мономов и коэффициентов между многочленами.
	*/
	inline void swap(CPlainPolynomial& p2){
		mons.swap(p2.mons);
		coeffs.swap(p2.coeffs);
	}

	///Сравнение старшинства полиномов
	bool operator< (const CPlainPolynomial& p2)const{
		return compareTo(p2)<0;
	}

	///Сравнение полиномов на равенство
	bool operator== (const CPlainPolynomial& p2)const{
		return compareTo(p2)==0;
	}

	//Данные операторы неэффективны для CPlainPolynomial
	/*
	///вычитание полиномов
	CPlainPolynomial operator-(const CPlainPolynomial& p)const;
	///сложение полиномов
	CPlainPolynomial operator+(const CPlainPolynomial& p)const;
	*/
	///умножение
	CPlainPolynomial operator*(const CPlainPolynomial& p)const;

	/**вычисляет хеш-функцию от полинома
	При вычисления хеша от полинома он формируется на основе значений хешей лишь некоторых старших мономов.
	Вычисление на большем числе мономов займёт много времени, но не имеет смысла,
	поскольку в большинстве случаев многочлены отличаются по первым мономам
	*/
	int hash() const{
		int res = 0;
		int optimalhashsz=size()/25;
		if (optimalhashsz<4) optimalhashsz=size()/5;
		if (optimalhashsz<3) optimalhashsz=size();
		for (int i=0;i<optimalhashsz;++i){
			res=(res*2+getMon(i).hash())^384923;
		}
		return res;
	}

	///домножение на моном
	CPlainPolynomial& operator*= (const CMonomial& givenMonomial){
		if (!(givenMonomial.isOne()))
		{
			for(m_iterator i = mons.begin();i!=mons.end(); ++i)
			{
				(*i)*=givenMonomial;
			}
		}
		return *this;
	}

	/**присваивание полинома, умноженного на моном.
	Присваивает данному полинома \a poly умноженный на \a givenMonomial без линих копирований
	*/
	void AssignMultiply(const CPlainPolynomial& poly,const CMonomial& givenMonomial){
		coeffs=poly.coeffs;
		mons.resize(poly.mons.size());
		m_iterator To;
		m_const_iterator From;
		for(From = poly.m_begin(), To=m_begin(); To!=m_end(); ++To, ++From){
			To->assignmul(*From,givenMonomial);
		}
	}

	///домножение на коэффициент
	CPlainPolynomial& operator*= (CModular c){
		if (c!=CModular(1)){
			c_iterator i;
			for(i = coeffs.begin(); i!=coeffs.end(); ++i){
				(*i)*=c;
			}
		}			
		return *this;
	}

	static void addCPolynomialMultiplied(const CPlainPolynomial& p1, const CPlainPolynomial& p2, CModular addCoeff, CPlainPolynomial& Result);
	static void addCPolynomialMultipliedWithMonom(const CPlainPolynomial& p1, const CPlainPolynomial& p2, CModular addCoeff, const CMonomial& addMon, CPlainPolynomial& Result);
};

/**Представление полинома со счётчиком ссылок и COW.
Интерфейс и реализация этого класса в точности идентична интерфейс CPlainPolynomial,
за исключением того испльзуется счётчик ссылок и данные копируются лишь при необходимости их изменить.
Несколько многочленов разделяют одни и те же данные в памяти, пока один из них не модифицируется.
Это позволяет уменьшить затраты памяти на хранение S-пар.
*/
class CRefPolynomial{
	explicit CRefPolynomial (const CPlainPolynomial& p):poly(new CPlainPolynomial(p))
	{}

	const CPlainPolynomial& cpoly()const{
		return const_cast<const CPlainPolynomial&>(*poly);
	}

  public:
	IntrusivePtr<CPlainPolynomial> poly;
	CRefPolynomial(){
		poly.makeDefault();
	}

	typedef CPlainPolynomial::MonomialInPoly MonomialInPoly;
	typedef CPlainPolynomial::c_iterator c_iterator;
	typedef CPlainPolynomial::m_iterator m_iterator;
	typedef CPlainPolynomial::c_const_iterator c_const_iterator;
	typedef CPlainPolynomial::m_const_iterator m_const_iterator;

	void resize(size_t sz){
		poly.makeUnique();
		poly->resize(sz);
	}
	
	size_t size()const{
		return poly->size();
	}

	//проверка, имеется ли хоть один терм
	bool empty() const {
		return poly->empty();
	}

	void normalize(){
		if (HC() == CModular(1)) return;
		poly.makeUnique();
		poly->normalize();
	}
	
	m_const_iterator m_end()const{
		return cpoly().m_end();
	}

	m_const_iterator m_begin()const{
		return cpoly().m_begin();
	}

	m_iterator m_end(){
		poly.makeUnique();
		return poly->m_end();
	}

	m_iterator m_begin(){
		poly.makeUnique();
		return poly->m_begin();
	}

	c_iterator c_end(){
		poly.makeUnique();
		return poly->c_end();
	}

	c_iterator c_begin(){
		poly.makeUnique();
		return poly->c_begin();
	}

	const MonomialInPoly& getMon(int i)const{
		return poly->getMon(i);
	}

	MonomialInPoly& getMon(int i){
		poly.makeUnique();
		return poly->getMon(i);
	}

	CModular getCoeff(int i)const{
		return poly->getCoeff(i);
	}

	CModular& getCoeff(int i){
		poly.makeUnique();
		return poly->getCoeff(i);
	}

	void addTerm(CModular c, const CMonomial& m){
		poly.makeUnique();
		poly->addTerm(c,m);
	}

	//добавляет терм в конец полинома
	void pushTermBack(CModular c, const CMonomial& m){
		poly.makeUnique();
		poly->pushTermBack(c,m);
	}

	void printPolynomial(std::ostream& output,ParserVarNames* names=0) const{
		poly->printPolynomial(output, names);
	}

	//возвращает старший моном полинома
	const MonomialInPoly& HM() const{
		assert(!empty());
		return poly->HM();
	}

	//возвращает коэффициент старшего монома
	CModular HC() const{
		assert(!empty());
		return poly->HC();
	}

	//отношение порядка на полиномах - для реализации polynomialComparator
	int compareTo(const CRefPolynomial& p2)const{
		return poly->compareTo(*p2.poly);
	}

	inline void swap(CRefPolynomial& p2){
		poly.swap(p2.poly);
	}

	bool operator< (const CRefPolynomial& p2)const{
		return (*poly)<(*p2.poly);
	}

	bool operator== (const CRefPolynomial& p2)const{
		return (*poly)==(*p2.poly);
	}

	CRefPolynomial operator-(const CRefPolynomial& p)const{
		CRefPolynomial result;
		CPlainPolynomial::addCPolynomialMultiplied(*poly, *p.poly, CModular(-1), *result.poly);
		return result;
		//return CRefPolynomial((*poly)-(*p.poly));
	}

	CRefPolynomial operator+(const CRefPolynomial& p)const{
		CRefPolynomial result;
		CPlainPolynomial::addCPolynomialMultiplied(*poly, *p.poly, CModular(1), *result.poly);
		return result;
		//return CRefPolynomial((*poly)+(*p.poly));
	}

	CRefPolynomial operator*(const CRefPolynomial& p)const{
		return CRefPolynomial((*poly)*(*p.poly));
	}

	CRefPolynomial GetLinearComb(const CRefPolynomial& p, CModular coeff, const CMonomial& mon)const{
		CRefPolynomial result;
		CPlainPolynomial::addCPolynomialMultipliedWithMonom(*poly, *p.poly, coeff, mon, *result.poly);
		return result;	
	}

	int hash() const{
		return poly->hash();
	}

	//умножение на моном
	CRefPolynomial& operator*= (const CMonomial& givenMonomial){
		/*poly.makeUnique();
		CPlainPolynomial* result=new CPlainPolynomial();
		(*poly)*=givenMonomial;*/

		CPlainPolynomial* result=new CPlainPolynomial();
		result->AssignMultiply(*poly,givenMonomial);
		poly.reset(result);

		return *this;
	}

	//умножение на коэффициент
	CRefPolynomial& operator*= (CModular c){
		poly.makeUnique();
		(*poly)*=c;
		return *this;
	}
};
/*
namespace std{
	//swap specialization
	inline void swap(CPolynomial& a, CPolynomial& b){
		a.swap(b);
	}
}
*/

///Испльзуемая реализация CPolynomial
typedef CRefPolynomial CPolynomial;
//typedef CPlainPolynomial CPolynomial;

} //namespace F4MPI
#endif
