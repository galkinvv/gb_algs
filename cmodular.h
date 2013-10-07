#ifndef CModular_h
#define CModular_h
#include <string>

/**
\file
Вычисления над конечным полем.
*/

namespace F4MPI{
/**Представление вычетов по модулю.
Реализует хранение и арифметические операции элементов поля вычетов по модулю.
Простой модуль, по которому ведутся вычисления хранится в глобальной переменной.
Ограничение данной реализации состоит в том, что его квадрат должен помещаться в переменную типа int.
*/
class CModular{
	///численное значение вычета (число от 0 до \c MOD-1)
	int val;
	///характеристика конечного поля
	static int MOD;

	///максимальный модуль для которого используется 32битная арифметика
	static const int MAXMUL32 = 40000;

  public:

	static bool NEEDMUL64;

	///возвращает характеристику поля
	static int getMOD(){
		return CModular::MOD;
	}

	///устанавливает характеристику поля равной \a m
	static void setMOD(int m){
		CModular::MOD = m;
		NEEDMUL64 = CModular::MOD > CModular::MAXMUL32;
	}

	///Конструктор по умолчанию обнуляет элемент
	CModular(){	val = 0; }

	///Конструктор, устанавливающий вычет равным \a givenValue
	CModular(int givenValue){
		val = givenValue % MOD;
		if (val<0) val+=MOD;
	}

	///сложение двух элементов
	friend CModular operator+(const CModular& A, const CModular& B);

	///унарный минус	
	CModular operator-(){
		CModular tmp(-val);
		return tmp;
	}

	///умножение двух элементов
	friend CModular operator*(const CModular& A, const CModular& B);	
	
	///сравнение двух элементов
	friend bool operator==(const CModular& A, const CModular& B);

	///сравнение двух элементов на неравенство
	friend bool operator!=(const CModular& A, const CModular& B);



	static inline void doMultiply(CModular& result, const CModular& m1, const CModular& m2){
		typedef long long lltype;
		//typedef int lltype;
		if (NEEDMUL64){
			result.val = (lltype(m1.val) * lltype(m2.val)) % lltype(CModular::MOD);
		}else{
			result.val = (m1.val * m2.val) % CModular::MOD;
		}
	}

	///домножение элемента
	void operator*=(const CModular& A){
		doMultiply(*this, *this, A);
	}

	///добавление элемента
	void operator+=(const CModular& A){
		val+=A.val;
		val%=CModular::MOD;
	}

	///возвращает строку, содержащую текстовое представление числа (просто число от 0 до \c MOD-1)
	std::string toString() const;

	///устанавливает значение вычета равным числу, записанному в строке givenValue
	void setValue(std::string givenValue);	
	
	/**\details
	Доступ к закрытой части, используется в исключительных случаях для оптимизации вычислений -
	несколько арифметических операций подряд можно выполнить и без приведения по модулю \c MOD,
	и взять по модулю лишь в самом конце
	*/
	int& pureint(){
		return val;
	}

	///явное преобразование к типу int
	int toint()const{
		return val;
	}

	/**
	отношение порядка на вычетах.
	Рассматривает вычеты как числа от 0 до \c MOD-1
	и сравнивает их.
	\retval число, определяющее результат сравнения:
	\arg -1: это число меньше \a m2
	\arg 1: это число больше \a m2
	\arg 0: это число равно \a m2
	*/
	
	int compareTo(const CModular& m2) const {
		int dif = val - m2.val;
		if(dif < 0)return -1;
		else if(dif==0)return 0;
		else return 1;		
	}

	///возвращает число, равное обратному к \a m2 в поле вычетов по модулю \c MOD
	static inline CModular inverseMod(const CModular& m2)
    {
		int x, y;
        extendedEuclid(m2.val, MOD, x, y);
		return CModular(x);
    }

private:
	/**
	расширенный алгоритм Евклида.
	Находит такие \a x и \a y, что <tt>a*x + b*y = НОД(a, b)</tt>.
	\retval НОД(a, b).
	*/
	static int extendedEuclid(int a, int b, int& x, int& y);
};
inline CModular operator+(const CModular& A, const CModular& B){
	CModular tmp;
	tmp.val = (A.val + B.val) % CModular::MOD;
	return tmp;
}

inline CModular operator*(const CModular& A, const CModular& B){
	CModular result;
	CModular::doMultiply(result, A, B);
	return result;
}

inline bool operator==(const CModular& A, const CModular& B){
	return A.val == B.val;
}

inline bool operator!=(const CModular& A, const CModular& B){
	return !(A == B);
}

} //namespace F4MPI
#endif
