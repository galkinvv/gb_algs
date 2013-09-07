#ifndef MemoryManager_h
#define MemoryManager_h
/**\file
Работа с памятью для мономов.
Определяет аллокатор для хранения мономов во внешней памяти
и массив хранящий элементы неизвестной при компиляции длинны для хранения мономов в нём
*/

#include <stack>
#include <vector>
#include <cassert>
#include <stdlib.h>
#include <string.h>
namespace F4MPI{
/**vector, оптимизированный для POD-типов.
в точности повторяет vector, за исключением того что конструкторы/деструкторы элементов никогда не вызываются,
и для копирования памяти при необходимости используется memcpy.
*/
template <typename T> class PODvector{
  public:
	typedef T& reference;
	typedef const T& const_reference;
	typedef T* pointer;
	typedef const T* const_pointer;
	typedef pointer iterator;
	typedef const_pointer const_iterator;
	typedef T value_type;
  private:
	pointer first,last,memlast;
	
	char* charfirst(){
		return reinterpret_cast<char*>(first);
	}

	void del(){
		delete[] charfirst();
	}

	void doreserve(size_t n){
		char *newc=new char[n*sizeof(T)];
		if (!empty()) memcpy(newc,charfirst(),size()*sizeof(T));
		del();
		first=reinterpret_cast<pointer>(newc);
		memlast=first+n;
	}
	
 public:
	void operator=(const PODvector& v){
		del();
		if (v.empty()) first=last=0;
		else{
			size_t n=v.size();
			first=reinterpret_cast<pointer>(new char[n*sizeof(T)]);
			//Выполнение этой строчки следует считать ошибкой, поскольку копирование лишнее
			memcpy(first,v.first,n*sizeof(T));
			last=memlast=first+n;
		}
	}
	
	PODvector(const PODvector& v):first(0),last(0){
		(*this)=v;
	}
	
	PODvector():
		first(0),last(0)
	{}

	~PODvector(){
		del();
	}
	
	bool empty()const{
		return (!first) || first==last;
	}

	void reserve(size_t n){
		if (capacity()>=n) return;
		size_t os=size();
		doreserve(n);
		last=first+os;
	}

	void swap(PODvector& v){
		std::swap(first,v.first);
		std::swap(last,v.last);
		std::swap(memlast,v.memlast);
	}

	void clear(){
		del();
		first=0;
		last=0;
	}

	void resize(size_t n){
		if (capacity()<n) doreserve(n);
		last=first+n;
	}

	T operator[](unsigned i)const{
		return first[i];
	}

	T& operator[](unsigned i){
		return first[i];
	}
	
	const_reference front()const{
		return *first;
	}

	reference front(){
		return *first;
	}

	const_reference back()const{
		return *(last-1);
	}

	const_iterator begin()const{
		return first;
	}

	iterator begin(){
		return first;
	}

	const_iterator end()const{
		return last;
	}

	iterator end(){
		return last;
	}

	size_t capacity(){
		return first ? (memlast-first) : 0;
	}
	
	size_t size()const{
		return last-first;
	}
};

/**vector c элементами неизвестной при компиляции длинны.
Контейнер аналогичный PODvector, размер элемента которого определяется во время выполнения.
Конструкторы/деструкторы элементов никода не вызываются, при необходимости копирование выполняется 
Документировано только то что отличает его от vector.
\tparam T тип ("неопределённой" длины), который хранится в контейнере. Функции контейнера возврвщают ссылки/указатели на него.
*/
template <typename T> class PODvecSize{
  static size_t iterator_step;///<размер элемента в байтах (глобальный параметр)
	/**итератор PODvecSize.
	С точки зрения интерфейса ничем не отличается от итереатора по vector.
	В операциях существенно используется рзмер элемента, сохранённый в PODvecSize::iterator_step.
	*/
	template <class ref, class ptr, class val_t, class Cont> class PODvecSize_iterator{
		typedef ptr pointer;
		typedef ref reference;

		val_t val;
		pointer getptr()const{
			return reinterpret_cast<pointer>(val);
		}

	  public:
		PODvecSize_iterator():val(0){}
		PODvecSize_iterator(val_t v):val(v){}

		reference operator*()const{
			return *getptr();
		}

		pointer operator->()const{
			return getptr();
		}

		bool operator!=(PODvecSize_iterator v){
			return val!=v.val;
		}

		size_t operator-(PODvecSize_iterator v)const{
			return (val-v.val)/Cont::iterator_step;
		}

		PODvecSize_iterator& operator++(){
			val+=Cont::iterator_step;
			return *this;
		}

		PODvecSize_iterator& operator--(){
			val-=Cont::iterator_step;
			return *this;
		}

		PODvecSize_iterator operator++(int){
			PODvecSize_iterator tmp=*this;
			++*this;
			return tmp;
		}

		PODvecSize_iterator operator--(int){
			PODvecSize_iterator tmp=*this;
			--*this;
			return tmp;
		}
	};
  
	typedef char *ptrhandler;
	typedef const char *const_ptrhandler;
	ptrhandler first, last, memlast;

  public:
	typedef T& reference;
	typedef T* pointer;
	typedef const T& const_reference;
	typedef const T* const_pointer;
	///итератор по контейнеру, учитывающий установленный размер элемента
	typedef PODvecSize_iterator<reference, pointer, ptrhandler, PODvecSize> iterator;
	///константный итератор по контейнеру, учитывающий установленный размер элемента
	typedef PODvecSize_iterator<const_reference, const_pointer, const_ptrhandler, PODvecSize> const_iterator;
	typedef T value_type;

	/**фиксирует размер элемента.
	Эту функцию можно вызывать для установки/изменения объёма памяти необходимого для типа T,
	только если ни одного объекта типа PODvecSize<T>, ни итератора по нему в данный момент не существует,
	иначе эти объекты/итераторы станут некорректно себя вести.
	*/
	static void setvalsize(size_t sz){
		iterator_step=sz;
	}
  private:
	
	void del(){
		delete[] first;
	}

	void doreserve(size_t n){
		char *newc=new char[n*iterator_step];
		if (!empty()) memcpy(newc,first,size()*iterator_step);
		del();
		first=newc;
		memlast=first+n*iterator_step;
	}
	
	const T& unsafeAt(unsigned i)const{
		return *reinterpret_cast<pointer>(first+i*iterator_step);
	}

 public:
	void operator=(const PODvecSize& v){
		del();
		if (v.empty()) first=last=0;
		else{
			size_t n=v.size();
			first=new char[n*iterator_step];
			//Выполнение этой строчки следует считать ошибкой, поскольку копирование лишнее
			//Однако сейчас копирование выпоняется при запихивании в PolynomMap
			memcpy(first,v.first,n*iterator_step);
			last=memlast=first+n*iterator_step;
		}
	}
	
	PODvecSize(const PODvecSize& v):first(0),last(0){
		(*this)=v;
	}
	
	PODvecSize():
		first(0),last(0)
	{}

	~PODvecSize(){
		del();
	}
	
	bool empty()const{
		return (!first) || first==last;
	}

	void reserve(size_t n){
		if (capacity()>=n) return;
		size_t os=size();
		doreserve(n);
		last=first+os*iterator_step;
	}

	void swap(PODvecSize& v){
		std::swap(first,v.first);
		std::swap(last,v.last);
		std::swap(memlast,v.memlast);
	}

	void clear(){
		del();
		first=last=0;
	}

	void resize(size_t n){
		if (capacity()<n) doreserve(n);
		last=first+n*iterator_step;
	}

	const T& operator[](unsigned i)const{
		assert(i>=0 && i<size());
		return unsafeAt(i);
	}

	T& operator[](unsigned i){
		assert(i>=0 && i<size());
		return *reinterpret_cast<pointer>(first+i*iterator_step);
	}
	
	reference front(){
		assert(!empty());
		return *reinterpret_cast<pointer>(first);
	}

	reference back(){
		assert(!empty());
		return *reinterpret_cast<pointer>(last-iterator_step);
	}

	const_iterator begin()const{
		return first;
	}

	iterator begin(){
		return first;
	}

	const_iterator end()const{
		return last;
	}

	iterator end(){
		return last;
	}

	size_t capacity(){
		return first ? (memlast-first)/iterator_step : 0;
	}
	
	size_t size()const{
		return (last-first)/iterator_step;
	}

	template <class T2> void push_back(const T2& nv){
		resize(size()+1);
		back()=nv;
	}

	template <class T2> void insert(iterator pos,const T2& nv){
		int ind=pos-begin();
		resize(size()+1);
		size_t sizeToMove = last-reinterpret_cast<const char*>(&unsafeAt(ind+1));
		if (sizeToMove != 0) memmove(&((*this)[ind+1]), &((*this)[ind]), sizeToMove);
		(*this)[ind]=nv;
	}
};

///Аллокатор памяти для мономов
class MemoryManager{
	static const int NBLOCKS = 1000;
	int BLOCKSIZE;
public:	
	typedef signed char Data;
private:
	typedef std::stack<Data* /*,vector<Data*>*/ > PoolStack;
	PoolStack pool;
	std::vector<Data*> allocated;
	void allocate()
	{
		Data* a = new Data[BLOCKSIZE*NBLOCKS];
		for(int i = 0; i<NBLOCKS; i++)
		{
			Data* p = a + BLOCKSIZE*i; 
			pool.push(p);
		}
		allocated.push_back(a);
	}
public:		
	/**фиксирует размер элемента.
	Эту функцию можно вызывать для установки/изменения объёма памяти необходимого для монома,
	только если pool и allocated пусты, то есть после вызова reset или до первого выделения памяти.
	*/
	void setSize(int sz){
		BLOCKSIZE = sz;
	}
	
	MemoryManager()
	{
		allocated.reserve(1000);
	}

	/**выделяет память.
	Перед первом вызовом для выделения памяти следует установить размер выделямого блока через setSize()
	*/
	Data* getMem()
	{
		if(pool.empty())
			allocate();
		Data* ret = pool.top();
		pool.pop();
		return ret;
	}

	/**Отдаёт память назад.
	Системе память не возвращается, а складывается в пул для возврата при последующих вызовах getMem()
	Для возврата системе следует использовать reset() (после того, как память больше не используется)
	*/
	void eraseMem(Data* p)
	{
		pool.push(p);
	}

	/**Возвращает память назад системе.
	Вызов reset необходим перед изменением размера выделямого блока
	(для того чтоб очистить пул от старых блоков неподходящего размера).
	После вызова reset использование памяти выделенной аллокатором, но не отданной ему,
	приведёт к неопределённому поведению.
	*/
	void reset()
	{
		//Освобождение всей выделеннную память
		for(int i = 0; i<allocated.size(); i++)
			delete[] allocated[i];		
		allocated.clear();
		//удаление её из пула
		pool=PoolStack();
	}

	/**деструктор освобождает всю выделенную память.
	После вызова деструктора использование памяти выделенной аллокатором, но не отданной ему,
	приведёт к неопределённому поведению.
	*/
	~MemoryManager()
	{
		reset();		
	}
};
} //namespace F4MPI
#endif
