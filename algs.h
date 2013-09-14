#ifndef Algs_h
#define Algs_h
/**
\file
Содержит обобщённые процедуры и типы.
Процедуры *_heap предназаначены для работы с кучей и являются дополнением к функциям работы с кучей из STL.
Предполагается, что куча хранится в виде бинарного дерева, записанного по строкам.
Поэтому необходимо, чтоб используемая STL-реализация функций \c make_heap/\c push_heap/\c pop_heap работала с таким же представлением.


IntrusivePtr\<T\> представляет собой реализацию указателя со счётчиком ссылок,
работающего с типами T унаследованными от IntrusiveRefCountBase
*/
namespace F4MPI{
/**
Восстанавливает кучу после изменения головного элемента.
Процедура предполагает что промежуток [\a first; \a last) образует пправильную кучу,
за исключением головного элемента *\a first, который был изменён.
После выполнения процедуры элементы переставлены местами так,
что весь промежуток является правильной кучей.
Единственной операцией, которая проводится с элементами является \c swap
\param first начало кучи
\param last конец кучи
\param comp упорядочивание (operator\<) в онтошении которого определяется понятие кучи
*/
template <class randomIterator, class Comparator>
void adjust_heap(randomIterator first,randomIterator last,Comparator comp){
	randomIterator cur=first;
	for(;;){
		const int num=cur-first;
		randomIterator ch=cur+num+1;
		if (ch<last && comp(*cur,*ch)){
			if (ch+1<last && comp(*ch,*(ch+1))) ++ch;
			swap(*cur,*ch);
			cur=ch;
			continue;
		}
		++ch;
		if (ch<last && comp(*cur,*ch)){
			swap(*cur,*ch);
			cur=ch;
			continue;
		}
		break;
	}
}

/**
Заменеят в куче головной элемент и корректирует кучу.
Процедура предполагает что промежуток [\a first; \a last) образует пправильную кучу.
После выполнения процедуры во множестве элементов головной элемент *\a first заменён на \a newhead
и элементы переставлены местами так, что весь промежуток является правильной кучей.
Единственной операцией, которая проводится с элементами является присваивание. 
Процедура позволяет избежать лишних копирований по сравнению с последовательной заменой головного элемента и вызовом adjust_heap()
\param first начало кучи
\param last конец кучи
\param comp упорядочивание (operator\<) в онтошении которого определяется понятие кучи
\param newhead элемент, который должен заменить головной во множестве элементов.
Он передаётся НЕ по ссылке специально, чтоб процедура могла быть вызвана с элементом из той же кучи в качестве аргумента
*/
template <class randomIterator, class Comparator, class T>
void adjust_heap_head(randomIterator first,randomIterator last,Comparator comp, T newhead){
	randomIterator cur=first;
	for(;;){
		const int num=cur-first;
		randomIterator ch=cur+num+1;
		if (ch<last && comp(newhead,*ch)){
			if (ch+1<last && comp(*ch,*(ch+1))) ++ch;
			*cur=*ch;
			cur=ch;
			continue;
		}
		++ch;
		if (ch<last && comp(newhead,*ch)){
			*cur=*ch;
			cur=ch;
			continue;
		}
		break;
	}
	*cur=newhead;
}

/*
template <class T> void adjust_heap2(T *first, T *last){
	T *cur=first;
	for(;;){
		const int num=cur-first;
		T *ch=cur+num+1;
		if (ch<last && (cur->curit->column > ch->curit->column)){
			if (ch+1<last && (ch->curit->column > (ch+1)->curit->column)) ++ch;
			swap(*cur,*ch);
			cur=ch;
			continue;
		}
		++ch;
		if (ch<last && (cur->curit->column > ch->curit->column)){
			swap(*cur,*ch);
			cur=ch;
			continue;
		}
		break;
	}
}
*/

/**
Базовый класс объектов с внедрённым счётчиком ссылок.
Для использования указателя IntrusivePtr необходимо,
чтоб класс на который он указывает публично наследовал IntrusiveRefCountBase
*/
struct IntrusiveRefCountBase{
	int refCounter;///<число ссылок
	///\details Конструктор, обнуляющий число ссылок у вновь создаваемого объекта
	inline IntrusiveRefCountBase():
		refCounter(0)
	{}
	/**\details Конструктор копирования, \b обнуляющий число ссылок у создаваемого объекта.
	Несмотря на то что это копия, число ссылок на неё изначально равно 0
	*/
	inline IntrusiveRefCountBase(const IntrusiveRefCountBase&):
		refCounter(0)
	{}
	/**\details Оператор присваивания, <b>не меняющий</b> число ссылок у объекта.
	Он ничего не делает, поскольку в результате присваивания меняется содержимое объекта,
	но не число ссылок на него
	*/
	void operator=(const IntrusiveRefCountBase&){}
};

/**
"Умный" указатель со счётчиком ссылок.
Для подсчёта ссылок используется внедрённый счётчик, при обнулении которого объект автоматически удаляется.
\tparam T тип, на который указывает указатель. Должен быть унаследован от IntrusiveRefCountBase
*/
template <class T>
class IntrusivePtr{
  private:
	/**
	удаляет ссылку указателя на объект.
	Число ссылок на объект уменьшается на 1, и, если после этого оно оказывается нулевым, объект удалется опреатором delete
	*/
	void release(){
		if (ptr){
			if (! --(ptr->refCounter)){
				delete ptr;
			}
			ptr=0;
		}
	}
	T* ptr;///<Объект, на который ссылается указатель
	///устанавливает ссылку указателя на объект *\a newPtr
	void set(T* newPtr){
		ptr=newPtr;
		if (ptr) ++(ptr->refCounter);
	}

  public:
	///Конструктор копирования добавляет новую ссылку на тот же объект
	IntrusivePtr(const IntrusivePtr& other){
		set(other.ptr);
	}
	
	///Конструктор добавляет новую ссылку на переданный объект
	IntrusivePtr(T* ptr){
		set(ptr);
	}

	///По умолчанию создаётся "нулевой", никуда не указывающий указатель
	IntrusivePtr():ptr(0){}

	///При присваивании заменяется указываемый объект и корректируется число ссылок
	IntrusivePtr& operator=(const IntrusivePtr& other){
		reset(other.ptr);
		return *this;
	}

	///разыменование - симулирует обычный указатель
	T& operator*() const{
		return get();
	}
	
	///разыменование - симулирует обычный указатель
	T* operator->() const{
		return &get();
	}

	///заменяет указываемый объект на \a newPtr и корректирует число ссылок
	void reset(T* newPtr){
		if (ptr==newPtr) return;//это необходимо, чтоб объект не терялся при присваивании единственного указателя самому себе
		release();
		set(newPtr);
	}

	///Операция разименования указателя в виде именованного метода
	T& get()const{
		//assert(ptr);
		return *ptr;
	}

	///Возвращает число указателей на тот объект, на который ссылается этот указатель
	int getRefCount() const{
		if (!ptr) return 0;
		return ptr->refCounter;
	}
	
	///Возвращает \c true, если на объект не ссылается других указателей
	bool unique() const{
		return getRefCount()<=1;
	}

	/**
	Делает данный указатель единственным, имеющим ссылку на объект.
	Если указатель уже уникален, ничего не делается.
	Если же есть ещё указатели на этот объект, то с помощью конструктора копирования создаётся копия объекта,
	и указатель становится единственным, на него указывающим	
	*/  
	void makeUnique(){
		if (!unique()){
			reset(new T(*ptr));
		}
	}

	/**Создаёт новый объект типа \a T, и устанавливает указатель на него.
	Новый объект, на который устанавливается указтель, создаётся конструктором по умолчанию
	*/
	void makeDefault(){
		reset(new T);
	}
	
	///меняет местами 2 указателя (не объекты!) эффективным образом
	void swap(IntrusivePtr& other){
		T* newPtr=other.ptr;
		other.ptr=ptr;
		ptr=newPtr;
	}

	///Деструктор автоматически убирает ссылку этого указателя на объект (и, соотвественно, при необходимости удаляет сам объект)
	~IntrusivePtr(){
		release();
	}
};
} //namespace F4MPI	

class NoCopy
{
	void operator=(const NoCopy&);
	NoCopy(const NoCopy&);
protected:
	NoCopy(){}
};

#endif
