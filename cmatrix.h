#ifndef CMatrix_h
#define CMatrix_h
/**
\file
Представление матриц и операции над ними.
Определят тип CMatrix, представляющий матрицу, CRow,
представляющий отдельную строку и структуру RowElement соответствующую одному ненулевому элементу строки.
*/

#include "settings.h"
#include "monomialmap.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <ostream>

namespace F4MPI{
/**
Структура для хранения элемента строки.
Поскольку матрица является разреженной, хранятся только ненулевые элементы строк.
Поэтому вместе с численным значением элемента хранится номер столбца, в котором он расположен.
Является POD-типом.
*/
struct RowElement{
	///номер столбца матрицы, в котором находится элемент
	int column;

	///значение элемента из поля Z_p (всегда ненулевое)
	CModular value;
};

/**
Тип представляющий строку матрицы.
Строка матрицы состоит из массива упорядоченных по возрастанию номера столбца переменных типа RowElement.
Этот массив представлен, как базовый класс PODvector<RowElement>,
от которого строка наследует стандартные операции над STL-контейенерами и тип итератора,
используемый для обхода элементов строки.
*/
struct CRow:public PODvector<RowElement>{
//struct CRow:public vector<RowElement >{

	///добавляет к строке rowTo строку rowFrom, домноженную на multBy
	void addRowMultypliedBy(const CRow& rowFrom, CModular multBy);

	///Базовый контейнер (массив) для хранения отдельных элементов
	typedef PODvector<RowElement> BaseContainer;

	typedef BaseContainer::iterator iterator;
	typedef BaseContainer::const_iterator const_iterator;

	using BaseContainer::swap;
	using BaseContainer::resize;
	using BaseContainer::begin;
	using BaseContainer::end;
	using BaseContainer::size;
	using BaseContainer::empty;
	
	/**Делает единичным ведущий элемент строки.
	Нормализует строку, т.е. домножает все элементы строки на коэффициент обратный ведущему (первому ненулевому).
	Таким образом ведущий элемент становится единичным.
	*/
	void normalize(){
		CModular c  = CModular::inverseMod(HC());	
		for(unsigned i = 0; i<size(); i++){
			(*this)[i].value = (*this)[i].value*c;
		}
	}

	/**
	возвращает численное значение ведущего элемента строки.
	Численное значение ведущего элемента строки определено только для строк с по крайней мере одним ненулевым элементом,
	поэтому функцию нельзя вызывать для нулевых строк.
	*/
	CModular HC()const{
		return begin()->value;
	}
	
	/**
	Возвращает столбец, содержащий ведущий элемент строки.
	Столбец, содержащий ведущий элемент строки определён только для строк с по крайней мере одним ненулевым элементом,
	поэтому функцию нельзя вызывать для нулевых строк.	
	*/
	int HM()const{
		return begin()->column;
	}

	/**
	функция сравнения столбцов двух элементов строки.
	Сравнение элементов по номеру столбца соответствует упорядочиванию, соответсвующему хранению элементов строки в массиве.
	Операция сравнения используется в качестве operator< для функций бинарного поиска элемента по номеру столбца.
	*/
	static bool rowElementComparator(const RowElement& a, const RowElement& b){
		return a.column<b.column;
	}
	
	/**
	возвращает элемент строки  в заданном столбце.
	Выполняет поиск в массиве ненулевых элементов элемента с заданным номером столбца.
	Если элемент найден, возвращает соответсвующее численное значение, иначе возвращает 0,
	так как в столбцах, которым не соответствует элемент массива, находятся нули.
	В силу упорядоченности массива используется бинарный поиск и время работы процедуры равно
	O(log(число ненулевых элементов)).
	\param columnID столбец, значение элемента в котором нужно вернуть.
	*/
	CModular getCoefByMonom(int columnID){
		//Если известно, что элемент в начале, то линейный поиск будет лучше
		//но это не всегда известно, поэтому пусть пока будет бинарный
	
		//порядок индексов соответствует порядку элементов в векторе, поэтому можно применить бинарный поиск
		if (columnID < this->front().column) return 0;
		RowElement searchfor;
		searchfor.column=columnID;
		iterator i;
		i = std::lower_bound(begin(), end(), searchfor, rowElementComparator);
		if (i!=end() && i->column==columnID){
			return i->value;
		}
		return 0;
	}
};

/** тип, представляющий матрицу.
Матрица состоит из массива строк #СRow.
Этот массив представлен, как базовый класс vector<СRow>,
от которого матрица наследует стандартные операции над STL-контейенерами и тип итератора,
используемый для обхода строк матрицы.
*/
class CMatrix:public std::vector<CRow>{
  public:
	///Тип строки матрицы
	typedef CRow Row;
  private:
	///Базовый контейнер (массив) для хранения строк
	typedef std::vector<Row> BaseContainer;
  public:
	using BaseContainer::back;
	using BaseContainer::clear;
	using BaseContainer::resize;
	using BaseContainer::reserve;
	using BaseContainer::begin;
	using BaseContainer::end;
	using BaseContainer::size;
	using BaseContainer::empty;

	typedef Row::iterator Rowit;
	typedef Row::const_iterator ConstRowit;
	typedef BaseContainer::iterator iterator;
	typedef iterator MatrixIterator;

private:
	/**соответствие между мономами и номерами столбцов.
	На этапе формирования матрицы из многочленов в матрице сохраняется соответствие между мономами и номерами столбцов,
	которое позже используется для обратного преобразования строк матрицы в многочлены.
	*/
	MonomialMap M;

public:

	/**возвращает ссылку на меременную, в которой должна быть сохранена статистика для данной матрицы.
	Статистика для данной матрицы сохраняется в последний элемент массива со статистикой,
	поэтому возвращаемая значение верно только до тех пор, пока сбор статистики не начнётся для следующей матрицы
	*/
	MatrixInfo& getMyStats(const F4AlgData* f4options);

	///Записывает статистику по текущему состоянию матрицы в заданную переменную
	void makeMatrixStats(SingleMatrixInfo &);

	/**
	Работа со статистикой до приведения матрицы к ступенчатому виду.
	Добавляет в массив статистических данных новый элемент, соответсвующий текущей матрице,
	записывает статистику по текущему (до редукции) состоянию, также выводит его в файл.
	Если статистика отключена, ничего не делает.
	*/
	void doMatrixStatsPre(const F4AlgData* f4options);

	/**
	Работа со статистикой после приведения матрицы к ступенчатому виду.
	Записывает статистику по текущему (после редукции) состоянию, время затраченное на редукцию, также выводит это в файл.
	Если статистика отключена, ничего не делает.
	*/
	void doMatrixStatsPost(const F4AlgData* f4options);

	const Row& getRow(int i)const{
		return (*this)[i];
	}
	
	///возвращает ссылку на строку матрицы по заданному номеру
	Row& getRow(int i){
		return (*this)[i];
	}

	const MonomialMap& getMonomialMap()const{
		return M;
	}

	///возвращает ссылку на соответствие между мономами и номерами столбцов
	MonomialMap& getMonomialMap(){
		return M;
	}

	///Печатает матрицу в заданный поток вывода \a output
	void printMatrix(std::ostream& output);

	void printMatrixInfo(FILE *output,const char* name, int ID=0);

	/**
	Проводит полную авторедукцию матрицы \a m.
	После проведения полной авторедукции матрица \a mимеет сильно ступенчатый вид -
	в столбце, содержащем ведущий элемент какой-то строки других ненулей нет.
	*/
	static void fullAutoReduce(CMatrix& m);

	/**
	Проводит полную авторедукцию набора строк.
	После проведения полной авторедукции набор строк [\a from;\a to) имеет сильно ступенчатый вид -
	в столбце, содержащем ведущий элемент какой-то строки других ненулей нет. Последовательная версия.
	*/
	static void fullAutoReduce(MatrixIterator from, MatrixIterator to);

private:

	/**
	Переносит из матрицы набор строк для отсылки на заданный процессор в начале редукции.
	Выбор строк для каждого из процессоров осуществляется так, что матрица представляется в виде
	объединения непересекающихся множеств строк по всем процессорам.
	Распределение строк по процессорам относительно равномерно.
	Из исходной матрицы строки удаляются (чтоб не пришлось их копировать).
	\param res матрица, в которую записываются выбранные строки
	\param id номер процесса, для которого нужно выбрать строки
	\param N общее рассматриваемое число строк в матрице
	\param usefulProcesses общее число процессов, которое будет использоваться для этой матрицы
	\param f4options прочие параметры (размеры блоков, и т.д.)
	*/
	void selectRowsForProcessor(CMatrix& res, int id, int N, int usefulProcesses,const F4AlgData* f4options);

	/**
	редуцирует строку \a row по строке \a by c известным обратным к ведущему.
	При редукции строки по \a by, последняя вычитается из \a row с таким коэффициентом,
	чтоб обнулить элемент \a row в ведущем столбце \a by.
	\param row строка, которая модифицируется редукцией по \a by
	\param by строка, выступающая в роли редуктора
	\param mb обратный к ведущему элементу строки \a by.
	Для нескольких вызовов с разными \a row строка \a by одинакова,
	и чтоб не вычислять обратный каждый раз заново он вычисляется и передаётся отдельно.
	*/
	static void reduceRowByRow(Row& row,const Row& by, CModular mb);

	///удаляет из матрицы пустые строки, начиная поиск таких с \a firstRowToSearch
	void eraseEmptyRows(unsigned firstRowToSearch=0);

	///редуцирует каждую строку набора [\a mfrom;\a mto) по строке by
	static void reduceRangeByRow(MatrixIterator mfrom, MatrixIterator mto,const Row& by);

	//для fastReduceRangeByRange. Документировано в реализации
	struct PairCurEnd;
	struct RowCoeff;
	static bool pairCurEndComparer(const PairCurEnd& b, const PairCurEnd& a);

	/**
	редуцирует подматрицу по подматрице.
	редуцирует каждую строку набора [\a mfrom;\a mto) по каждой из строк [\a byfrom;\a byto).
	Функция пытается использовать блочность операции и эффективно использовать кеш процессора.
	Требуется, чтоб редуцирующий набор [\a byfrom;\a byto) был полностью авторедуцирован.
	*/
	static void fastReduceRangeByRange(MatrixIterator mfrom, MatrixIterator mto, MatrixIterator byfrom, MatrixIterator byto);
public:
	/**
	\details
	используется в реализации fastReduceRangeByRange() с зафиксированным параметром \a willNotOverflow
	\tparam willNotOverflow параметр, определяющий может ли произойти переполнение int
	при простом суммировании всех слагаемых, образующих новый элемент без взятия их по модулю.
	Если он равен \c true, то взятие по модулю происходит только после суммирования,
	иначе взятие по модулю производится для каждого слагаемого (являющегося произведением двух чисел).
	Вынесен в шаблонный параметр, чтоб быть известным во время компиляции, что позволит компилятору породить более оптимальный код в самом главном внутреннем цикле программы.
	*/
	template <bool willNotOverflow> static void fastReduceRangeByRangeConstOverflow(MatrixIterator mfrom, MatrixIterator mto, MatrixIterator byfrom, MatrixIterator byto);

	/**
	редуцирует подматрицу по матрице с заданным размером блока.
	Редуцирует каждую строку набора [\a mfrom;\a mto) по каждой из строк матрицы \a by.
	Если аргумент \a blocksize равен 1, то использзуется версия редуцирующая строки по одной через reduceRowByRow().
	При больших значениях размера блока используется блочный вариант fastReduceRangeByRange().
	Требуется, чтоб редуцирующая матрица \a by была полностью авторедуцирована.
	*/
	static void reduceRangeByMatrix(MatrixIterator mfrom, MatrixIterator mto, CMatrix& by, int blocksize);

	/**\details
	редуцирует подматрицу по подматрице (версия для обратного хода).
	Отличается от обычной тем, что есть дополнительное предусловие:
	после редукции гарантированно не могут получиться нулевые строки
	*/
	static void backReduceRangeByRange(MatrixIterator mfrom, MatrixIterator mto, MatrixIterator byfrom, MatrixIterator byto);

	/**Параллельная версия метода Гаусса.
	Метод Гаусса, распаралелленный с помощью MPI.
	\param f4options параметры, переданные алгоритму F4 (размеры блоков, опции статистики, и т.д.)
	\param doAutoReduce указывает на необходимость проведения авторедукции (доведения до сильно ступенчатого вида)
	*/
	void MPIDiagonalForm(const F4AlgData* f4options,int doAutoReduce=true);

	///Параллельно приводит матрицу к сильно ступенчатому виду (с помощью MPIDiagonalForm())
	void toDiagonalNormalForm(const F4AlgData* f4options);
	
	///Параллельно приводит матрицу к обычному ступенчатому виду (с помощью MPIDiagonalForm())
	void toRowEchelonForm(const F4AlgData* f4options);
};
} //namespace F4MPI
#endif
