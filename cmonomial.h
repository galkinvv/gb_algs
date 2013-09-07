#ifndef CMonomial_h
#define CMonomial_h
/** \file
Представление и операции над мономом многочлена.
Определяет класс CMonomial, хранящий данные в отдельно выделенной памяти
и CInternalMonomial, хранящий данные в памяти непояредственно после *this
*/

#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cassert>

#include "globalf4.h"

namespace F4MPI{

typedef std::vector<std::string> ParserVarNames;

/**определение базовых операций над мономамами.
Этот класс лишь оределяет операции над мономами, которые ипользуют отнаследованные от него классы,
но не имеет никаких данных и нестатических функций.
Существенным является то, что предполагаемое представление данных монома
содержит не только степени по каждой из переменной, но и заранее вычисленную суммарную степень,
что позволяет ускорить многие операции.
*/
class CMonomialBase{
  public:
	enum Order{
		lexOrder = 1,      ///<код лекксикографического порядка на мономах (lex)
		deglexOrder = 2,     ///<порядок на мономах: по суммарной степени, потом лексикографически (deglex)
		degrevlexOrder = 3,    ///<порядок на мономах: по суммарной степени, потом обратоно-лексикографически (degrevlex)
		blklexOrder = 4    ///<порядок на мономах: два обратоно-лексикографических блока переменных упорядоченные между собой лексикографически
	};

	///код порядка на мономах (глобальный параметр)
	static Order order;       

	///параметризация порядка (глобальный параметр, не для всех порядков)
	static int orderParam;

	///число переменных (глобальный параметр)
	static int theNumberOfVariables;

	/**\details
	число байт, занимаемых данными монома.
	Данные сотоят из последовательно расположенной суммарной степени и набора степеней по отдельным переменным.
	Поэтому degreessize должно быть равно (\c #theNumberOfVariables + 1)*(число байт, занимаемое одной степенью)
	*/
	static int degreessize;
  
	static const int MAX_DEGREE = 100;
	
	///возвращает код порядка на мономах
	static Order getOrder();

	///устанавливает порядок с кодом \a ord на мономах.
	static void setOrder(Order ord, int orderParameter);
    
	///устанавливает число переменных равным \a n (и degreessize соотвественно)
	static void setNumberOfVariables(int n){
		//static const int varsInHash=sizeof(int)/sizeof(Deg);
		//warning! This assume little endian arch;
		//static const int[] hashMasks={0,0x000000FF,0x0000FFFF};

		theNumberOfVariables=n;
		degreessize=(n+1)*sizeof(Deg);
	}

	//typedef int Deg;
	///Тип для хранения степеней - как суммарной, так и каждой переменной
	typedef signed char Deg;//На очень больших примерах может понадобиться заменить signed char на signed short
  protected:

	///Инициализирует данные монома \a degrees значениями по умолчанию (нулевой степенью)
	static inline void initmondefault(Deg *degrees){
		for(int i = 0; i<theNumberOfVariables+1; ++i)
			degrees[i] = 0;
	}

	/**Инициализирует моном по набору степеней переменных
	Инициализирует данные монома \a degrees значениями по заданному набору степеней переменных \a ptr
	и вычисляет суммарную степень
	*/
	static inline void initmonfromptr(Deg *degrees,const Deg *ptr){
		for(int i = 0; i<theNumberOfVariables; ++i)
			degrees[i+1] = ptr[i];
		setDegree(degrees);
	}

	///Инициализирует моном \a degrees по заданному моному \a dgoriginal
	static inline void initmonfrom(Deg *degrees,const Deg *dgoriginal){
		for(int i = 0; i<theNumberOfVariables+1; ++i)
			degrees[i] = dgoriginal[i];
	}

	///корректиует моном \a degrees - устанавливает суммарную степень равной сумме степеней переменных
	static inline void setDegree(Deg *degrees)
	{
		degrees[0]=0;
		for(int i =1; i< theNumberOfVariables+1; ++i)
			degrees[0]+=degrees[i];
	}

	/** НОД мономов.
	записывает в \a degreesres данные наибольшего общего делителя \a degrees1 и \a degrees2
	*/
	static inline void gcd (const Deg *degrees1, const Deg *degrees2, Deg *degreesres)
	{
		for(int i = 1; i<theNumberOfVariables+1; ++i)
			degreesres[i]=std::min(degrees1[i],degrees2[i]);
		setDegree(degreesres);
	}

	static inline void checkOverflowInLcm(const Deg *degrees1, const Deg *degrees2)
	{
		if (int(degrees1[0]) + int(degrees2[0]) > MAX_DEGREE) throw std::runtime_error(std::string("possible degree overflow in monomial lcm: ") + toString(degrees1) + " * " + toString(degrees2));
	}
	
	/** НОК мономов.
	записывает в \a degreesres данные наименьшего общего кратного \a degrees1 и \a degrees2
	*/
	static inline void lcm (const Deg *degrees1, const Deg *degrees2, Deg *degreesres)
	{
		checkOverflowInLcm(degrees1, degrees2);
		for(int i = 1; i<theNumberOfVariables+1; ++i)
			degreesres[i]=std::max(degrees1[i],degrees2[i]);
		setDegree(degreesres);
	}

	static inline void checkOverflowInMul(const Deg *degrees1, const Deg *degrees2)
	{
		if (int(degrees1[0]) + int(degrees2[0]) > MAX_DEGREE) throw std::runtime_error(std::string("degree overflow in monomial multiplication: ") + toString(degrees1) + " * " + toString(degrees2));
	}
	
	/** умножение мономов.
	записывает в \a degreesres данные поризведения \a degrees1 и \a degrees2
	*/
	static inline void mul(const Deg *degrees1, const Deg *degrees2, Deg *degreesres){
		checkOverflowInMul(degrees1, degrees2);
		for(int i = 0; i<theNumberOfVariables+1; ++i)
			degreesres[i]=degrees1[i]+degrees2[i];
	}
	
	/** in-place доножение монома.
	записывает в \a degrees1 данные поризведения \a degrees1 и \a degrees2
	*/
	static inline void mulby(Deg *degrees1, const Deg *degrees2){
		checkOverflowInMul(degrees1, degrees2);
		for(int i = 0; i<theNumberOfVariables+1; ++i)
			degrees1[i]+=degrees2[i];
	}

	/** попытка деления мономов.
	записывает в \a degreesres данные результата деления \a degrees1 на \a degrees2.
	Если деление невозможно возвращает false, а записанный в degreesres результат неопределён
	*/
	static inline bool tryDiv(const Deg *degrees1,const Deg *degrees2, Deg *degreesres){
		for(int i = 0; i<theNumberOfVariables+1; ++i)
		{
			if ((degreesres[i]=degrees1[i]-degrees2[i]) < 0) return false;
		}
		return true;
	}
	
	/** деление мономов.
	записывает в \a degreesres данные результата деления \a degrees1 на \a degrees2.
	Процедура требует, чтоб деление было возможно, иначе полученный результат будет некорректен
	*/
	static inline void div(const Deg *degrees1,const Deg *degrees2, Deg *degreesres){
		for(int i = 0; i<theNumberOfVariables+1; ++i)
		{
			degreesres[i]=degrees1[i]-degrees2[i];
			assert(degreesres[i]>=0);
		}
	}
	
	/** in-place деление монома.
	записывает в \a degrees1 данные результата деления \a degrees1 на \a degrees2.
	Процедура требует, чтоб деление было возможно, иначе полученный результат будет некорректен
	*/
	static inline void divby(Deg *degrees1,const Deg *degrees2){
		for(int i = 0; i<theNumberOfVariables+1; ++i)
		{
			degrees1[i]-=degrees2[i];
			assert(degrees1[i]>=0);
		}
	}

	/**Сравнивает поднаборы monomFrom и monomTo степеней degrees1 и degrees2 по порядку degrevlex*/
	static int compareDegRevLex(const Deg *degrees1, const Deg *degrees2, int monomFrom, int monomTo){
		int s1 = 0, s2 = 0;
		for (int i = monomFrom; i != monomTo; ++i){
			s1 += degrees1[i];
			s2 += degrees2[i];
		}
		int degdif = s1 - s2;
		if (degdif>0) return 1;
		if (degdif<0) return -1;
		for (int i = monomTo -1; i >= monomFrom; --i) {
			int dif = degrees1[i] - degrees2[i];
			if (dif<0) return 1;//reversed
			if (dif>0) return -1;
		}
		return 0;
	}
	
	/**отношение порядка на мономах.
	Сравнение происходит в соотвествии установленным на данный момент порядком.
	\retval число, определяющее результат сравнения:
	\arg -1: \a degrees1 \< \a degrees2
	\arg 1: \a degrees1 \> \a degrees2
	\arg 0: \a degrees1 равно \a degrees2
	*/
	static inline int compareTo(const Deg *degrees1, const Deg *degrees2, Order orderToCompareWith){
		switch (orderToCompareWith){
			case lexOrder:
				{
					for (int i = 1; i < theNumberOfVariables+1; ++i) {
						int dif = degrees1[i] - degrees2[i];
						if (dif>0) return 1;
						if (dif<0) return -1;
					}
					break;
				}
			case deglexOrder:
				{ 
					int degdif = degrees1[0] - degrees2[0];
					if (degdif>0) return 1;
					if (degdif<0) return -1;
					for (int i = 1; i < theNumberOfVariables+1; ++i) {
						int dif = degrees1[i] - degrees2[i];
						if (dif>0) return 1;
						if (dif<0) return -1;
					}
					break;
				}
			case degrevlexOrder:
				{
					int degdif = degrees1[0] - degrees2[0];
					if (degdif>0) return 1;
					if (degdif<0) return -1;
					for (int i = theNumberOfVariables; i > 0; --i) {
						int dif = degrees1[i] - degrees2[i];
						if (dif<0) return 1;//reversed
						if (dif>0) return -1;
					}
					break;
				}
			case blklexOrder:
				{
					int res1 = compareDegRevLex(degrees1, degrees2, 1, orderParam+1);
					if (res1) return res1;
					return compareDegRevLex(degrees1, degrees2, orderParam+1, theNumberOfVariables+1);
					break;
				}
		}
		return 0;
	}

	/**возвращает хеш-функцию на мономе \a degrees
	хеш-функция выбрана так, чтоб минимизировать вероятность совпадения кодов для двух различных мономов \b небольшой степени с небольшим числом переменных.
	Это достинается за счёт того, что информация о степенях "размазана" по всему числу.
	Старшие биты более зависят почти только от степеней старших (первых) переменных, а младшие - от тех и других.
	Таким образом, если 2 монома впервые отличаются отличаются в степени i-й переменной,
	то, вероятно что значение хеш функции будет отличаться в битах, зависящих от перемеенных {1, ... , i}.
	*/
	static int hash(const Deg *degrees){
		int res = 0;
		//int res = 123;
		for(int i = 0; i<theNumberOfVariables; i++){//на последний элемент можно не смотреть, поскольку он является функцией от предыдущих
			res*=11;
			res^=degrees[i];
			//res = (res*11 + degrees[i])^142857151;
		}
		return res;
	}

	/**сравнение мономов на равенство.
	Сравнение на равенство может оказаться эффективнее вызова compareTo(), поскольку понятие равенства не зависит от порядка и проще проверяется.
	*/
	static bool isEqual(const Deg *degrees1, const Deg *degrees2)
	{
		//return memcmp(degrees1,degrees2,(1+theNumberOfVariables*sizeof(Deg)))==0; //странно, это оказалось медленнее
		for(int i = 0; i<=theNumberOfVariables; ++i)//на последний элемент можно не смотреть, поскольку он является функцией от предыдущих
			if (degrees1[i]!=degrees2[i]) return false;
		return true;
	}

	///Проверка на возможность деления
	static bool divisibleBy(const Deg *degrees1,const Deg *degrees2)
	{
		for(int i = 0; i<theNumberOfVariables+1; ++i)
			if(degrees1[i] < degrees2[i])
				return false;
		return true;
	}  

	/**возвращает строку, содержащую текстовое представление монома.
	\param degrees данные монома
	\param names соответствие между индексами и текстовыми именами переменных.
	Если передан нулевой указатель (по умолчанию), используются стандартные имена x1 .. xN
	*/
	static std::string toString(const Deg *degrees, ParserVarNames* names=0){
		std::ostringstream out;
		bool was = false;
		for(int i = 1; i<theNumberOfVariables+1; ++i){									
			if(degrees[i]>0){
				if(was)	out<<"*";
				out<<varName(i, names);
				if(degrees[i]>1){
					out<<"^"<<int(degrees[i]);
				}
				was = true;
			}
		}
		return out.str();
	}
  public:	
	static std::string varName(int varIdxFrom1, ParserVarNames* names);
};

/** хранит данные монома в памяти, выделенной отдельно.
В качестве менеджера памяти используется MemoryManager.
Память под данные выделяется при создании и удаляется при разрушении объекта этого класса.
Чтоб избежать утечек или повторного освобождения памяти, копирование этого класса запрещено.
Предполагает, что MemoryManager настроен так, чтоб выделять память нужного размера.
*/
class MonomialExternalPlacing{
	MonomialExternalPlacing(const MonomialExternalPlacing&);//запрет на явное копирование класса
	void operator=(const MonomialExternalPlacing&);//запрет на явное копирование класса
  protected:
	typedef CMonomialBase::Deg Deg;
	///указатель на внешнюю память для хранения данных монома.
	Deg* degrees;

	MonomialExternalPlacing(){
		degrees = (Deg*)globalF4MPI::MonomialAllocator.getMem();
		//degrees = new Deg[CMonomialBase::theNumberOfVariables+1];
	}

	~MonomialExternalPlacing()
	{
		globalF4MPI::MonomialAllocator.eraseMem(degrees);	
		//delete[] degrees;
	}
/*
  public:
	void swap(const MonomialExternalPlacing& other){
		std::swap(degrees,other);
	}
*/
};

/** хранит данные монома в памяти по адресу this.
Объекты этого класса не создаются, не копируются и не удаляются,
но в указатель на этот тип преобразуются указатели на реально выделенную внешнюю память
(это делается в классе PODvecSize).
*/
class MonomialInternalPlacing{
	MonomialInternalPlacing();//запрет на явное создание переменной класса
	MonomialInternalPlacing(const MonomialInternalPlacing&);//запрет на явное копирование класса
	void operator=(const MonomialInternalPlacing&);
  protected:
	typedef CMonomialBase::Deg Deg;
	/**\details
	массив "указывающий" на данные монома; 
	его адрес как раз совпадает с this (началом объекта POD-типа), реальный размер выделенной памяти равен
	CMonomialBase::degreessize, но за это отвечает процедура, преобразовавшая указатель на память в указатель на MonomialInternalPlacing.
	*/
	Deg degrees[2];//реальный размер массива равен theNumberOfVariables+1
  public:
	static int realSize(){
		return CMonomialBase::degreessize;
	}
};

template <class Placing> class MonomialWithPlacing;

/**Обычное представление монома.
Это представление монома использует внешнюю память для хранения данных и с переменными этого типа можно делать всё что угодно.
*/
typedef MonomialWithPlacing<MonomialExternalPlacing> CMonomial;

/**Представление монома в многочлене.
Переменная этого типа сущетствует \b только, как элемент контейнера PODvecSize<CInternalMonomial>, который и выделил для неё память. В остальных местах используются только указатели/ссылки на неё.
*/
typedef MonomialWithPlacing<MonomialInternalPlacing> CInternalMonomial;
/**Реализация интерфейса монома.
Класс наследует реализацию из CMonomialBase, представление из \a Placing и определяет общий интерфейс доступа к мономам.
*/
template <class Placing>
class MonomialWithPlacing:public Placing, public CMonomialBase{
	using CMonomialBase::Deg;
public:
	/**\details
	Placing должен иметь член-данное degrees, ковертируемое в тип Deg*.
	Полученный указатель рассматривается, как память нужного размера для хранения данных монома.
	*/
	using Placing::degrees;
	MonomialWithPlacing():Placing(){
		initmondefault(degrees);
	}

	MonomialWithPlacing(const MonomialWithPlacing& givenMonomial):Placing(){
		initmonfrom(degrees,givenMonomial.degrees);//копируются сами данные, а не указатели на них
	}
	/**\details
	Шаблонный "конструктор копирования" позволяет преобразовывать друг в друга мономы,
	использующие различные способы рзмещения данных.
	*/
	template <class Placing2> MonomialWithPlacing(const MonomialWithPlacing<Placing2>& givenMonomial)
	{	
		initmonfrom(degrees,givenMonomial.degrees);
	}

	void operator=(const MonomialWithPlacing& givenMonomial){
		initmonfrom(degrees,givenMonomial.degrees);
	}

	/**\details
	Шаблонный оператор присваивания позволяет присваивать друг другу мономы,
	использующие различные способы рзмещения данных.
	*/
	template <class Placing2> void operator=(const MonomialWithPlacing<Placing2>& givenMonomial){
		initmonfrom(degrees,givenMonomial.degrees);
	}	

	///создание монома на основе вектора степеней переменных
	MonomialWithPlacing(const std::vector<Deg>& initMonomial)
	{
		initmonfromptr(degrees,&initMonomial[0]);
	}
	
	bool isOne()const
	{
		return getDegree() == 0;
	}
	///возвращает суммарную степень монома
	int getDegree()const
	{
		return degrees[0];
	}

	///возвращает  степень переменной i-й переменной
	int getDegree(int i) const
	{
		return degrees[i+1];
	}

	///возвращает НОД 2х мономов
	template <class Placing3> static const CMonomial gcd (const MonomialWithPlacing<Placing>& firstMonomial, const MonomialWithPlacing<Placing3> & secondMonomial)
	{
		CMonomial res;
		CMonomialBase::gcd(firstMonomial.degrees,secondMonomial.degrees,res.degrees);
		return res;
	}

	///возвращает НОК 2х мономов
	template <class Placing3> static const CMonomial lcm (const MonomialWithPlacing<Placing>& firstMonomial, const MonomialWithPlacing<Placing3> & secondMonomial)
	{
		CMonomial res;
		CMonomialBase::lcm(firstMonomial.degrees,secondMonomial.degrees,res.degrees);
		return res;
	}

	///возвращает произведение 2х мономов
	template <class Placing3> friend CMonomial operator*(const MonomialWithPlacing<Placing>& firstMonomial, const MonomialWithPlacing<Placing3> & secondMonomial){
		CMonomial res;
		CMonomialBase::mul(firstMonomial.degrees,secondMonomial.degrees,res.degrees);
		return res;
	}
	
	/**присваивает поризведение 2х мономов.
	Единственное отличие от использования операторов = и * в том, что лишних копирований не происходит.
	*/
	template <class Placing2, class Placing3> void assignmul(const MonomialWithPlacing<Placing2>& firstMonomial, const MonomialWithPlacing<Placing3> & secondMonomial){
		CMonomialBase::mul(firstMonomial.degrees,secondMonomial.degrees,degrees);
	}
	
/*	///возвращает результат деления 2х мономов, не проверяя возможно ли оно
	template <class Placing3> friend CMonomial operator/(const MonomialWithPlacing<Placing>& firstMonomial, const MonomialWithPlacing<Placing3> & secondMonomial){
		CMonomial res;
		CMonomialBase::div(firstMonomial.degrees,secondMonomial.degrees,res.degrees);
		return res;
	}
*/
	///домножение на моном
	template <class Placing2> void operator*= (const MonomialWithPlacing<Placing2>& givenMonomial)
	{
		CMonomialBase::mulby(degrees,givenMonomial.degrees);
	}


	/*///деление на моном; возможность не проверяется, в случае невозможности результат некорректен.
	template <class Placing2> void operator/= (const MonomialWithPlacing<Placing2>& givenMonomial)
	{
		CMonomialBase::divby(degrees,givenMonomial.degrees);
	}*/

	///сравнение монома по упорядочиванию с \a monomialToCompare; возвращает -1 (меньше), 0 (равно) или 1(больше)
	template <class Placing2> inline int compareTo(const MonomialWithPlacing<Placing2>& monomialToCompare, Order orderToCompareWith)const{
		return CMonomialBase::compareTo(degrees, monomialToCompare.degrees, orderToCompareWith);
	}

	///Сравнение с применением глобального порядка
	template <class Placing2> inline int compareTo(const MonomialWithPlacing<Placing2>& monomialToCompare)const{
		return compareTo(monomialToCompare, CMonomialBase::order);
	}

	///возвращает хеш-код монома
	int hash() const {
		return CMonomialBase::hash(degrees);
	}

	///сравнения на равенство; эффективней compareTo()
	template <class Placing2> bool operator== (const MonomialWithPlacing<Placing2>& givenMonomial) const
	{
		return CMonomialBase::isEqual(degrees,givenMonomial.degrees);
	}

	///сравнения на неравенство; эффективней compareTo()
	template <class Placing2> bool operator!= (const MonomialWithPlacing<Placing2>& givenMonomial) const
	{
		return !(*this == givenMonomial);
	}

	///Проводит попытку деления данного монома на моном \a m и возвращает \c true в случае успеха. Результат успешного деления сохраняется в \a result, иначе значение \a result неопределено
	template <class Placing2, class Placing3> bool tryDivide(const MonomialWithPlacing<Placing2>& m, MonomialWithPlacing<Placing3>& result) const
	{
		return CMonomialBase::tryDiv(degrees, m.degrees, result.degrees);
	}  

	///возвращает \c true, если данный моном делится на моном \a m
	template <class Placing2> bool divisibleBy(const MonomialWithPlacing<Placing2>& m) const
	{
		return CMonomialBase::divisibleBy(degrees,m.degrees);
	}  

	/**возвращает строку, содержащую текстовое представление монома.
	\param names соответствие между индексами и текстовыми именами переменных.
	Если передан нулевой указатель (по умолчанию), используются стандартные имена x1 .. xN
	*/
	std::string toString(ParserVarNames* names=0) const {
		return CMonomialBase::toString(degrees,names);
	}
};

class CMonomialForLabel: private CMonomial
{
	static Order orderForLabel;
public:
	template <class Placing2> explicit CMonomialForLabel(const MonomialWithPlacing<Placing2>& givenMonomial)
		:CMonomial(givenMonomial)
	{}
	CMonomialForLabel(){}

	const CMonomial& mon()const
	{
		return *this;
	}
	
	template <class Placing2> void operator *=(const MonomialWithPlacing<Placing2>& givenMonomial)
	{
		(CMonomial&)(*this)*=(givenMonomial);
	}
	template <class Placing2> CMonomialForLabel operator *(const MonomialWithPlacing<Placing2>& givenMonomial)const
	{
		return CMonomialForLabel(mon() * givenMonomial);
	}
	int compareToByLabelOrder(const CMonomialForLabel& monomialToCompare)const{
		return compareTo(monomialToCompare, orderForLabel);
	}
	template <class Placing2> int compareToByLabelOrder(const MonomialWithPlacing<Placing2>& monomialToCompare)const{
		return compareTo(monomialToCompare, orderForLabel);
	}
};

/**Возвращает текстовое название упорядочивания \a orderID.
возможные значения: <tt>"lex", "deglex", "degrevlex"</tt>.
*/
const std::string getMonomialOrderName(CMonomialBase::Order orderID);
} //namespace F4MPI
#endif
