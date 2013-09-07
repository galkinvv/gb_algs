
#ifndef MonomialMap_h
#define MonomialMap_h
/**
\file
Работа с множеством мономов.
*/

#include "cpolynomial.h"

#include <string>
#include <ostream>
#include <algorithm>
#include <unordered_map>
namespace F4MPI{
///сравнение мономов по порядку на них
struct monomialComparator{
	bool operator()(const CMonomial& m1,const CMonomial& m2) const{
		return m1.compareTo(m2)<0;
	}
};

/**\details
содержит указатель на объект, для которого operator== переопределён так,
что сравнивает равенство указываемых объектов.*/
template <class T> struct ptrAsObject{
	ptrAsObject(){}
	explicit ptrAsObject(const T& obj):ptr(&obj){}
	const T* ptr;
	bool operator==(const ptrAsObject& other)const{
		return (*other.ptr)==*ptr;
	}
};

/**тип "указателя" на моном.
Используется для помещения в hash_map для того, чтоб избежать лишних копирований непосредственно объектов.
*/
typedef ptrAsObject<CMonomial> MonomialPtr;
/**вычисляет хеш по монома по #MonomialPtr.
Используется в качестве параметра шаблона hash_map, определяющего хеш-функцию
*/
struct monomialHasher{
	size_t operator()(const CMonomial& v)const{
		return v.hash();
	}
	///Возвращает хеш монома, на который уквзывает \a v
	size_t operator()(const MonomialPtr& v)const{
		return v.ptr->hash();
	}
	
	//needed only for MSVC
	bool operator()(const CMonomial& u,const CMonomial& v)const{
		return monomialComparator()(u,v);
	}

	/**\details
	Сравнивает на \< мономы, на которые указывают \a u и \a v.
	Этот оператор используется рализауией hash_map в MSVC
	*/
	bool operator()(const MonomialPtr& u, const MonomialPtr& v)const{
		return monomialComparator()(*u.ptr,*v.ptr);
	}
	static const size_t bucket_size = 1;
#ifdef _DEBUG
	static const size_t min_buckets = 64;
	//иначе страшно тормозит в отладке
#else
	static const size_t min_buckets = 262144;
#endif
};


/**множество мономов, сопоставленных целым числам.
Класс используется как для упорядоченного соотношения между мономами и номерами столбцов матрицы,
так и просто для хранения множества мономов не использующего соответствие с числами.
*/
class MonomialMap{
	///Тип, определяющий множество с отображением
	typedef std::unordered_map<MonomialPtr, int, monomialHasher> Container;
	///отображение мономов в целые числа
	Container M;//Возможно следует менять размер хеша в зависимости от реальной потребности

	///тип для отображения чисел на мономы.
	typedef std::vector<CMonomial> RevMStorgae;
	///отображение целых чисел в мономы
	RevMStorgae revM;
	/**\details
	Контейнер используемый для хранения мономов, указатели на которые кладут в hash_map.
	Дек используется потому что контейнер не должен перемещать элементы при добавлении новых,
	чтоб не портились указатели.
	*/
	std::deque<CMonomial> monomialStorage;

public:
	
	typedef Container::iterator iterator;
	
	Container& getMap(){
		return M;
	}

	///возвращает размер множества
	int size(){
		return int(M.size());
	}

	///возвращает \c true, если множество содержит моном \a m
	bool containsMonomial(const CMonomial& m)const{		
		return  M.find(MonomialPtr(m))!=M.end();
	}    

	///добавляет моном \a m во множество, если он там еще не содержится
	void storeMonomial(const CMonomial& m){
/*
		Container::value_type objToInsert(MonomialPtr(m),monomialStorage.size());
		pair<iterator,bool> findingResult=M.insert(objToInsert);
		if (findingResult.second){//Такого элемента ещё не было
			monomialStorage.push_back(m);
			const_cast<MonomialPtr&>(findingResult.first->first)=MonomialPtr(monomialStorage.back());
		}
		return findingResult.first->second;
*/
		iterator it = M.find(MonomialPtr(m));
		if(it == M.end()){
			int newID=monomialStorage.size();
			monomialStorage.push_back(m);
			M[MonomialPtr(monomialStorage.back())] = newID;
		}
	}

	/**возвращает число, сопоставленное моному.
	Упорядочение возвращаемых чисел соотвествуют порядку на мономах, только если после последнего вызова UpdateForUsingReversed()
	новых мономов не добавлялось.
	Передаваемый моном \a должен присутствовать в множестве, иначе поведение не определено.
	*/
	int getMonomialID(const CMonomial& m)const{
		return M.find(MonomialPtr(m))->second;
	}

	///возвращает \c true если множество пустое, \c false иначе
	bool empty() const {
		return M.empty();
	}

	///удаляет моном \a m из множества.
	void erase(const CMonomial& m){
		M.erase(MonomialPtr(m));
	}
	
	///возвращает некоторый моном из множества.
	const CMonomial& selectMonoimial(){		
		return *(M.begin()->first.ptr);
	}
	
	/**возвращает моном, соответствующий числу \a i.
	Для использзования этой функции необходимо чтоб \a i соответсвовало корректному номеру монома.
	То есть множество ен должно было меняться с тех пор, как номер был получен от getMonomialID(), вызванной после UpdateForUsingReversed().
	*/
	const CMonomial& getMonomialRev(int i)const{
		return revM[i];			
	}

	/**Подготавливает множество мономов к использованию отображений
	Перенумерует мономы по их порядку (в обратную сторону, чтоб столбцы в матрице шли по возростанию).
	Функции getMonomialID() и getMonomialRev() корректно работают только после вызова этой функции.
	Полученные числа становятся некорректны, если множество изменить после этого.
	Более того, функция в силу особенности реализации:
	\arg добавит назад в множество все удалённые элементы.
	\arg повторный вызов этой функции приведёт к полной очистке множества.
	*/
	void UpdateForUsingReversed(){
		revM.clear();
		M.clear();//все мономы будут вноситься заново
		revM.resize(monomialStorage.size());
		std::copy(monomialStorage.begin(),monomialStorage.end(),revM.begin());
		monomialStorage.clear();//теперь мономы хранятся только в revM

		sort(revM.rbegin(),revM.rend(),monomialHasher());
		int num=0;
		for (RevMStorgae::iterator i=revM.begin();i!=revM.end();++i,++num){
			M[MonomialPtr(*i)]=num;
		}
	}

	///выводит в \a output сопоставление между мономами и числами
	void PrintMap(std::ostream& output){
		for(iterator i = M.begin(); i!=M.end(); ++i){
			CMonomial Ml = *(i->first.ptr);
			std::string s = Ml.toString();
			output << "M[" << s << "] = " << i->second << "\n";
		}
	}

	///добавляет мономы полинома \a p, не содержащиеся в \a notThere, во множество \a M.
	friend void storeNotProcessedMonomialsFromPoly(MonomialMap &M, MonomialMap& notThere, const CPolynomial& p);

	///добавляет мономы полинома \a p во множество \a M.
	friend void storeMonomialsFromPoly(MonomialMap &M, const CPolynomial &p);

};
} //namespace F4MPI
#endif
