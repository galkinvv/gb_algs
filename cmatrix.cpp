/**
\file
Реализация представления матриц и операций над ними*/

#include "cmatrix.h"
#include "mpimatrix.h"
#include "cmodular.h"
#include "f4main.h"
#include "algs.h"
#include "matrixinfoimpl.h"
#include <iostream>
#include <algorithm>
#include <limits>
#define MPICH_SKIP_MPICXX
#include <mpi.h>
using namespace std;
namespace F4MPI{
/**
описывает состояние итерации по строке.
текущий и конечный итераторы строки,
описывающие состояние итерируемой строки с возможностью проверки окончания итрерации.
используется в fastReduceRangeByRange()
*/
struct CMatrix::PairCurEnd{
	ConstRowit curit;///<текущее положение итерации по строке
	ConstRowit lastit;///<итератор, указывающий за конец строки

	/**Информация о номере строки, которой эта пара итераторов соответствует.
	Для редуцирующих строк: номер строки в наборе \n
	Для редуцируемых строк: (номер строки в ответе, куда записывается результат редукции этой строки) + (число редуцирующих строк).
	
	Таким образом, если значение этой переменной меньше, чем (число редуцирующих строк), то итераторы попрождены редуцирующей строкой, иначе редуцируемой.
	*/
	int row;
	///Инициализирует структуру для начала итерации по данной строке \a r
	void initfrom(Row& r){
		curit=r.begin();
		lastit=r.end();
	}
	///Определяет, является ли текущая позиция итерирования последней позицией со значащим элементом
	bool isprelast(){
		return curit+1==lastit;
	}
};


/**элемент матрицы коэффтцтентов прибавляемых строк.
описывает вхождение редуцирующей строки в ответ, используется в fastReduceRangeByRange()
*/
struct CMatrix::RowCoeff{
	CModular coeff;///<коэффициент
	Row::value_type** pRowit;///<указатель на итератор по строке ответа, в текущую позицию которого дописывается текущий элемент редуцирующей строки
};

/**
Оператор сравнения состояний итерации строк по их текущим столбцам.
Определяет упорядочивание на состояниях итерации строк, ориентируясь на текущий рассматриваемый столбец.
Используется в fastReduceRangeByRange() для создания кучи из элементов PairCurEnd и выбора среди них строки с наименьшим нерассмотренным столбцом
*/
bool CMatrix::pairCurEndComparer(const PairCurEnd& b, const PairCurEnd& a){
	return a.curit->column<b.curit->column;
};

 
void CMatrix::fastReduceRangeByRange(MatrixIterator mfrom, MatrixIterator mto, MatrixIterator byfrom, MatrixIterator byto){
	bool willNotOverflow;
	if (CModular::NEEDMUL64) willNotOverflow = false;
	else{
		int maxVal=CModular::getMOD()*CModular::getMOD();
		int maxInt=std::numeric_limits<int>::max();
		int bysz=byto-byfrom;
		willNotOverflow=maxVal<(maxInt/(1+bysz));
	}
	if (willNotOverflow){
		fastReduceRangeByRangeConstOverflow<true>(mfrom,mto,byfrom,byto);
	}else{
		fastReduceRangeByRangeConstOverflow<false>(mfrom,mto,byfrom,byto);
	}
}

template <bool willNotOverflow> 
void CMatrix::fastReduceRangeByRangeConstOverflow(MatrixIterator mfrom, MatrixIterator mto, MatrixIterator byfrom, MatrixIterator byto){
	
	const int bysz=byto-byfrom;
	const int msz=mto-mfrom;
	int maxresrowsize=0;//оценка на длину результирующих строк
	//Например максимальный элемент строки среди всех строк:
	for (MatrixIterator i=byfrom;i!=byto;++i){
		maxresrowsize=max(maxresrowsize,i->back().column+1);
	}
	for (MatrixIterator i=mfrom;i!=mto;++i){
		if(i->size()>0){
			maxresrowsize=max(maxresrowsize,i->back().column+1);
		}
	}
	//(возможно сумма размеров byfrom..byto и максимального из mfrom..mto будет лучше)

	const int res1size=maxresrowsize+1;
	//матрица для записи результатов
	char *const resultsmem=new char[msz*res1size*sizeof(typename Row::value_type)];
	typename Row::value_type *const results=reinterpret_cast<typename Row::value_type *const>(resultsmem);
	
	vector<typename Row::value_type*> res(msz);//итераторы, соответствующие текущему элементу в результатах

	for (int j=0;j<msz;++j){
		res[j]=results+j*res1size;//инициализируем res указателями на строки для записи результата
	}
	
	vector<bool> isModified(msz,false);//определяет, изменяется ли заданная строка

	vector<int> res2origRow;//соотвествие строк результата, и исходных
	res2origRow.reserve(msz);
	vector<int> origRow2res(msz,-1);//соотвествие строк результата, и исходных

	vector<bool> reducerUsed(bysz);
	//хранится в линейном виде чтоб эффективней использовать кеш
	char* const cmem=new char[bysz*(msz+1)*sizeof(RowCoeff)];//матрица коэффициентов
	RowCoeff* const c=reinterpret_cast<RowCoeff* const>(cmem);//матрица коэффициентов
	//для каждой строки byto..byfrom представлена набором пар (коэффициент,указатель на результат)
	//завершается парой с нулевым коэффициентом

	//Её инициализация
	for (int i=0;i<bysz;++i){
		CModular mb = CModular::inverseMod((byfrom+i)->HC());
		int monID = (byfrom+i)->HM();
		int baseindexinc=i*(msz+1);
		int curindexinc=baseindexinc;
		for (int j=0;j<msz;++j){
			if ((mfrom+j)->size()==0) continue;
			//коэффициент с которым i-ю редуцирующую строку прибавляют к j-й редуцируемой
			CModular k = -(mfrom+j)->getCoefByMonom(monID)*mb;
			if (k!=0){
				if (!isModified[j]){//Этой исходной строке ещё нет соотвествия в результате
					int nextFreeRow=res2origRow.size();
					//установим прямое и обратное соотвествие
					origRow2res[j]=nextFreeRow;
					res2origRow.push_back(j);
					isModified[j]=true;
				}

				c[curindexinc].pRowit=&res[origRow2res[j]];
				c[curindexinc].coeff=k;
				++curindexinc;
			}
		}
		reducerUsed[i]=curindexinc!=baseindexinc;
		c[curindexinc].coeff=0;//завершающая пара
	}
	
	const int usedResRows=res2origRow.size();
	PairCurEnd *const m=new PairCurEnd[bysz+msz];//итераторы, соответствующие текущему и последнему элементу в исходных строках и редуцирующих строках
	PairCurEnd *mlast=m;//Конец кучи
	for (int j=0;j<bysz;++j){
		if (!reducerUsed[j]) continue;
		mlast->initfrom(*(byfrom+j));//инициализируем первую часть m указателями на редуцирующие строки
		mlast->row=j;//определяет где в матрице c искать различные коэффициенты
		++mlast;
	}

	for (int j=0;j<msz;++j){
		if (!isModified[j]) continue;
		const int resRowNum=origRow2res[j];
		res[resRowNum]->value=0;//поставим нулевой коеффициент в начало
		mlast->initfrom(*(mfrom+j));//инициализируем вторую часть m указателями на исходные строки
		mlast->row=resRowNum+bysz;//Строка в которую следует копировать элементы
		++mlast;
	}

	static const int lastcolorignal=-1;
	int lastcol=lastcolorignal;
	make_heap(m,mlast,pairCurEndComparer);
	
	
	
	for(;;){
		int smallestcolumn=lastcolorignal;
		if (mlast!=m){
			smallestcolumn=m->curit->column;
		}
		if (smallestcolumn!=lastcol){//перешли в следующую колонку, нужно подвинуть итераторы результата
			for (int j=0;j<usedResRows;++j){
				if (res[j]->value.toint() && ((res[j]->value).pureint()%=CModular::getMOD())){//получили ненулевой элемент в строке
					res[j]->column=lastcol;//запишем номер столбца
					++(res[j]);//продвинем указатель
					res[j]->value=0;//обнулим элемент, чтоб на следующем шаге не принять за новый
				}
			}
			lastcol=smallestcolumn;
		}
		if (mlast==m) break;//больше ничего не осталось
		const CModular l=m->curit->value;//число в строке, которое нужно прибавить к результату с тем или иным коэффициентом
		const int resRow = m->row - bysz;
		if (resRow>=0){//это редуцируемая строка; её нужно просто прибавить
			res[resRow]->value+=l;
		}else{//это редуцирующая строка; её нужно прибавить к нескольким с разными коэффициентеми
			const RowCoeff *f=c+(m->row)*(msz+1);
			if (!willNotOverflow){//известно на момент компиляции
				while (f->coeff.toint()){//нультерменированный список коэффициентов и ссылок
					(*(f->pRowit))->value+=f->coeff*l;//Деление на MOD здесь необходимо лишь для больших модулей
					++f;
				}
			}else{
				while (f->coeff.toint()){//нультерменированный список коэффициентов и ссылок
					((*(f->pRowit))->value).pureint()+=f->coeff.toint()*l.toint();
					++f;
				}
			}
		}
		if (!m->isprelast()){
			++(m->curit);
			//добавим изменённый элемент вместо вершины
			adjust_heap_head(m,mlast,pairCurEndComparer,m[0]);
		}else{
			/*
			pop_heap(m,mlast,columncomparer2());
			--mlast;
			*/
			--mlast;
			//добавим последний выкинутый элемент вместо вершины
			adjust_heap_head(m,mlast,pairCurEndComparer,mlast[0]);
		}
	}
	delete[] cmem;
	delete[] m;

	for (int j=0;j<msz;++j){
		if (!isModified[j]) continue;
		const int resRow=origRow2res[j];
		(mfrom+j)->clear();
		//Установим реальный размер полученного результата
		(mfrom+j)->resize(res[resRow]-(results+resRow*res1size));
		//Запишем результат на место исходных данных
		memcpy(&(mfrom+j)->front(),results+resRow*res1size,(mfrom+j)->size()*sizeof((mfrom+j)->front()));
	}
	delete[] resultsmem;
}

 void CMatrix::printMatrixInfo(FILE *output,const char* name,int ID){
	SingleMatrixInfo info;
	makeMatrixStats(info);
	fprintf(output,"matrix %s   \tnrows: %d  columns: %d  elements: %d  fill percent: %f%%\n",
		name,
		info.rows,
		info.columns,
		info.elems,100.*info.filling());
	fflush(output);
}

void CMatrix::printMatrix(ostream& output){
	CMatrix& A=*this;
	int M = 0;
	for(unsigned i = 0; i<A.size(); i++)
		for(unsigned j = 0; j<A[i].size(); j++)
			M = max(M, A[i][j].column);
	
	for(unsigned i = 0; i<A.size(); i++){
		int prev  = -1;
		for(unsigned j = 0; j<A[i].size(); j++){
			int k = 1;
			while(A[i][j].column!=prev+k){
				output << " 0";
				k++;
			}
			prev = A[i][j].column;
			output <<" "<< A[i][j].value.toint();
		}
		while(prev++<M) output << " 0";
		output << "\n";
	}
	output<<endl;
}


void CMatrix::reduceRangeByMatrix(MatrixIterator mfrom, MatrixIterator mto,CMatrix& by, int blocksize){
	int sz=mto-mfrom;
	if (blocksize==1){
		for (MatrixIterator i=by.begin();i!=by.end();++i){
			reduceRangeByRow(mfrom,mto,*i);
		}
	}else{
		for (int i=0;i<sz;i+=blocksize){
			fastReduceRangeByRange(
					mfrom+i,
					mfrom+min(i+blocksize,sz),
					by.begin(),
					by.end()
					);
		}
	}
}

///Набор времён, замеренных разными способами
struct MPITimeMesurement{
	double mw;///<Время, замеренное MPI_Wtime
	//AbsoluteTime cl;///<Время, замеренное системными функциями
};

///Возвращает набор замеров текущего времени
const MPITimeMesurement getMPITimeMesurement(){
	MPITimeMesurement res={MPI_Wtime()/*,getTime()*/};
	return res;
}

//int MPI_PROCESS_CIRCULATE_ORDER=0;

/**порядок перебора процессов в MPI.
Определяет порядок перебора процессов при распределении нагрузки в методе Гаусса,
предоставляя функции перехода "к следующему" и "к предыдущему", использкющие внутренний счётчик текущего шага
*/
class ResIdCalculator{
	int nProc;///<общее число процессов
	int step;///<номер текущего шага
	int getID()const{
		int direction=0;
/*		switch(MPI_PROCESS_CIRCULATE_ORDER){
			case 0:{
				direction=0;
				break;
			}
			case 1:{
				//Порядок обхода процесооров такой: 1....N,N....1,1....N,N....1, и т.д.
				direction=(step%(2*nProc))/(nProc);
				break;
			}
			case 2:{
				//Порядок обхода процесооров такой: 1....N,N....1,N....1,1....N,1....N,N....1,N....1,1....N, и т.д.
				direction=((step+nProc)%(4*nProc))/(2*nProc);
				break;
			}
			case 3:{
				direction=1;
				break;
			}
		}
*/
		if (!direction) return step%nProc;
		else return (nProc-1-step%nProc);
	}
  public:
	///Конктруктор определяет числом процессов, которые надо перебирать
	ResIdCalculator(int numberOfProcesses):nProc(numberOfProcesses),step(-1){}
	///Перейти к следующему процессу и вернуть его номер
	int getNext(){
		++step;
		return getID();
	}
	///Перейти к предыдущему процессу и вернуть его номер
	int getPrev(){
		--step;
		return getID();
	}
	///Вернуть номер текущего процесса
	int getCur()const{
		return getID();
	}
};

///номер строки матрицы, с которой начинает блок \a n
int getNthBlockStart(int n,const F4AlgData* f4options){
	return n*f4options->MPIBlockSize;
}

///Число блоков на которое разбивается матрица с \a matrixSize строк при использовании \a usefulProcesses процессов
int getTotalBlocksOnAllProcessors(int matrixSize, int usefulProcesses,const F4AlgData* f4options){
	int fullBlocks=0;
	while (getNthBlockStart(fullBlocks+1,f4options)*usefulProcesses<=matrixSize) ++fullBlocks;
	int remBlock=matrixSize-getNthBlockStart(fullBlocks,f4options)*usefulProcesses;
	int remFullLines=remBlock/usefulProcesses;
	int remExtraLines=remBlock%usefulProcesses;
	int res=fullBlocks*usefulProcesses;
	if (remFullLines>0){
		res+=usefulProcesses;
	}else{
		res+=remExtraLines;
	}
	return res;
}

void CMatrix::selectRowsForProcessor(CMatrix& res, int id, int matrixSize, int usefulProcesses,const F4AlgData* f4options){
	res.clear();
	int fullBlocks=0;
	while (getNthBlockStart(fullBlocks+1,f4options)*usefulProcesses<=matrixSize) ++fullBlocks;
	int remBlock=matrixSize-getNthBlockStart(fullBlocks,f4options)*usefulProcesses;
	int remFullLines=remBlock/usefulProcesses;
	int remExtraLines=remBlock%usefulProcesses;
	bool hasExtraLine=(id<remExtraLines);
	res.resize(getNthBlockStart(fullBlocks,f4options)+remFullLines+hasExtraLine);
	for(int level=0;level<=fullBlocks;++level){
		int rowFrom=getNthBlockStart(level,f4options)*usefulProcesses;
		int thisBlockSize;
		if (level==fullBlocks){
			thisBlockSize=remFullLines;
			if (hasExtraLine){
				rowFrom+=id;
			}else{
				rowFrom+=remExtraLines;
			}
		}else{
			thisBlockSize=getNthBlockStart(level+1,f4options)-getNthBlockStart(level,f4options);
		}
		rowFrom+=thisBlockSize*id;
		int rowTo=getNthBlockStart(level,f4options);
		int rowToLast=min(getNthBlockStart(level+1,f4options),int(res.size()));
		for(;rowTo!=rowToLast;++rowTo,++rowFrom){
			res[rowTo].swap((*this)[rowFrom]);
		}
	}
}

/**
Набор аргументов для метода Гаусса.
Содержит набор аргументов метода Гаусса, которые при вызове параллельной рализации известны толко в главном процессе,
и должны быть разосланы на все остальные.
*/
struct MPIDiagonalFormOptions{
	int matrixSize;///<число строк в матрице
	int doAutoReduce;///<необходимость доведения до сильно ступенчатого вида 
	int usefulProcesses;///<число процессов, реально задействованных в операции
};

/**
Рассылает/принимает аргументы метода Гаусса из главного процесса в остальные.
Перед рассылкой параметров вызов из главного процесса выводит остальные процессы из спячки,
посылая им флаг продолжения вычислений.
После этого остальные процессы также заходят в wakeUpSyncWithOtherProcesses() и происходит пересылка аргументов всем процессам.
*/
void wakeUpSyncWithOtherProcesses(const F4AlgData* f4options, MPIDiagonalFormOptions& options){
	if (!f4options->thisProcessRank){
		//Нулевой процесс должен разбудить остальных
		int finish=false;//флаг окончания не установлен
		MPI_Bcast(&finish,1,MPI_INT,0,MPI_COMM_WORLD);//Выведем остальные процессы из ожидания
	}
	MPI_Bcast(&options,sizeof(options),MPI_CHAR,0,MPI_COMM_WORLD);//Раздадим всем опции
}

/**
Параллельный прямой ход метода Гаусса.
Прямой ход меотда Гаусса проводится для матрицы состоящей представленной в виде объединения \a mymat по всем процессорам.
Результат представляется в виде объединения \a resultmat по всем процессорам
\param mymat строки исходной матрицы. По окончанию неопределено (портится)
\param resultmat по окончании содержит строки полученной матрицы, хранящиеся на этом процессоре
\param f4options параметры F4 (размеры блоков, параметры коммуникаций)
\param reduceOptions исходные аргументы метода Гаусса
\param resultblocksizes размеры блоков строк результата, полученных из блоков строк исходной матрицы
\param resIdCalc отвечает за порядок перебора процессов
\param commUseful MPI-коммуникатор, объединяющий все процессы, реально участвующие в проведении редукции этой матрицы
\retval суммарное число строк в результате по всем процессорам
*/
int forwardGaussElimination(
		CMatrix& mymat, CMatrix& resultmat,
		const F4AlgData* f4options, const MPIDiagonalFormOptions& reduceOptions,
		vector<int>& resultblocksizes, ResIdCalculator& resIdCalc,
		MPI_Comm commUseful){
	int totalLinesDone=0;
	vector<int> rowBlockStarts;//содержит номера начал разосланных блоков строк в mymat
	rowBlockStarts.reserve(2+reduceOptions.matrixSize/f4options->MPIBlockSize);//число блоков + 1 на обозначение конца
	int lastNotProcessedBlock=0;//Блок, относительно которого в следующий раз будет происходить редукция
	for(;;){
		int nextBlockStart=getNthBlockStart(rowBlockStarts.size(),f4options);
		if (nextBlockStart>=mymat.size()){
			rowBlockStarts.push_back(mymat.size());//В конце добавляется номер "за последним блоком"
			break;
		}
		rowBlockStarts.push_back(nextBlockStart);
	}
	int totalSteps=getTotalBlocksOnAllProcessors(reduceOptions.matrixSize, reduceOptions.usefulProcesses,f4options);
	CMatrix myMatNew;//матрица, где формируется результат редукции mymat по extramat
	CMatrix extramat;//матрица-редуктор
	int resid;//Процессор на который нужно поместить редуцирующую строку в готовые.
	int getid;//Процессор на котором находится редуцирующая строка.
	for(int step=0;step<totalSteps;++step){
		getid=resid=resIdCalc.getNext();
		if (f4options->thisProcessRank==getid){//Разобьём mymat на 2 матрицы:
			//extramat - по которой будет происходить редукция
			//myMatNew - оставшиеся строки. Попутно выкинем нулевые
			extramat.reserve(rowBlockStarts[lastNotProcessedBlock+1]-rowBlockStarts[lastNotProcessedBlock]);
			myMatNew.reserve(mymat.size()-rowBlockStarts[lastNotProcessedBlock]);
			for (int i=lastNotProcessedBlock;i+1<rowBlockStarts.size();++i){
				int from=rowBlockStarts[i];
				int to=rowBlockStarts[i+1];
				//выкинуть нулевые строки
				CMatrix* destMat;
				if (i==lastNotProcessedBlock){//Определение того, в какую из матриц пойдёт блок
					destMat=&extramat;
				}else{
					destMat=&myMatNew;
					rowBlockStarts[i]=myMatNew.size();
					//Предыдущие блоки полностью уложены; Запомним начало текущего
				}
				for(int j=from;j!=to;++j){
					CMatrix::Row& r=*(mymat.begin()+j);
					if (r.empty()) continue;
					destMat->push_back(CMatrix::Row());
					r.swap(destMat->back());
				}
			}
			rowBlockStarts.back()=myMatNew.size();//За концом последнего блока
			mymat.swap(myMatNew);
			myMatNew.clear();
			++lastNotProcessedBlock;
			//приведём extramat к сильно ступенчатому виду
			CMatrix::fullAutoReduce(extramat);
			//PrintMatrix(extramat);
		}
		if (reduceOptions.usefulProcesses>1) bcastMat(getid,extramat,f4options->thisProcessRank,commUseful,f4options->MPIUseBigSends);//разослать extramat с getid на все процессоры

		if (f4options->thisProcessRank==resid){
			resultblocksizes.push_back(extramat.size());
		}
		if (extramat.size()==0){
			continue;
		}
		//Отредуцировать
		CMatrix::reduceRangeByMatrix(mymat.begin()+rowBlockStarts[lastNotProcessedBlock],mymat.end(),extramat,f4options->innerGaussBlockSize);

		totalLinesDone+=extramat.size();
		if (f4options->thisProcessRank==resid){
			//На этом процессоре положим результирующую строку в результаты
			for (CMatrix::iterator i=extramat.begin();i!=extramat.end();++i){
				resultmat.push_back(CMatrix::Row());
				resultmat.back().swap(*i);
			}
		}
		extramat.clear();
	}
	return totalLinesDone;
}

/**
Параллельный обратный ход метода Гаусса.
Обратный ход меотда Гаусса проводится для матрицы состоящей представленной в виде объединения \a resultmat по всем процессорам.
Результат представляется в виде matrix на основном процессе.
\param resultmat строки исходной матрицы. По окончанию неопределено (портится)
\param matrix по окончании на основном процессе содержит строки полученной матрицы, на остальных пусто
\param f4options параметры F4 (размеры блоков, параметры коммуникаций)
\param reduceOptions аргументы метода Гаусса с \c matrixSize равным значению числа строк в объединении resultmat по всем процессорам
\param resultblocksizes размеры блоков строк, из которых составлены данные resultmat
\param resIdCalc отвечает за порядок перебора процессов
\param commUseful MPI-коммуникатор, объединяющий все процессы, реально участвующие в проведении редукции этой матрицы
*/
void backGaussElimination(
		CMatrix& resultmat, CMatrix& matrix,
		const F4AlgData* f4options, const MPIDiagonalFormOptions& reduceOptions,
		vector<int>& resultblocksizes, ResIdCalculator& resIdCalc,
		MPI_Comm commUseful){
	int lastreducibleline=resultmat.size();//последняя неавторедуцированная строка в resultmat на текущем процессоре
	int getid=resIdCalc.getCur();
	//Теперь getid содержит номер последнего процесса, который записывал строки => с него и начинать
	int linesremains=reduceOptions.matrixSize;
	matrix.reserve(reduceOptions.matrixSize);
	while(linesremains>0){
		int bsize;
		if (f4options->thisProcessRank==getid){//редуцирующий блок строк находится на данном процессоре
			bsize=resultblocksizes.back();//размер этого блока
			resultblocksizes.pop_back();
			lastreducibleline-=bsize;//последний блок на данный момент авторредуцирован
			//место под строки результата
			matrix.resize(matrix.size()+bsize);
			CMatrix::iterator dest=matrix.end()-bsize;
			CMatrix::iterator src=resultmat.begin()+lastreducibleline;
			CMatrix::iterator srcend=src+bsize;
			//Перенос строк из промежуточного resultmat в matrix
			for(;src!=srcend;++dest,++src){
				(*src).swap(*dest);
				//Возможно это и вообще не нужно, а даже если нужно, то не здесь...
				//dest->normalize();
			}
			//рассылка блока на все процессоры, если это необходимо
			if (reduceOptions.usefulProcesses>1){
				bcastSendSubMatrix(getid,matrix.end()-bsize,matrix.end(),commUseful,f4options->MPIUseBigSends);
			}
		}else{
			//место под строки результата
			bsize=bcastRecvToMatrix(getid,matrix,commUseful,f4options->MPIUseBigSends);			
		}
		//авторедукция строк с текущего процесса по блоку
		CMatrix::backReduceRangeByRange(
				resultmat.begin(),
				resultmat.begin()+lastreducibleline,
				matrix.end()-bsize,
				matrix.end()
				);
		if (f4options->thisProcessRank!=0){
			//результат больше не нужен на этом процессоре
			matrix.clear();
		}
		//общее число неавторедуцированных строк
		linesremains-=bsize;
		//следующий процессор с которого будет браться блок
		getid=resIdCalc.getPrev();
	}
	
}

void CMatrix::MPIDiagonalForm(const F4AlgData* f4options,int doAutoReduce){
	MPITimeMesurement ttt[4];
	ttt[0]=getMPITimeMesurement();
	int ID=f4options->thisProcessRank;
	CMatrix &matrix=*this;
	MPIDiagonalFormOptions reduceOptions;
	reduceOptions.doAutoReduce=doAutoReduce;
	reduceOptions.matrixSize=matrix.size();
	reduceOptions.usefulProcesses=0;
	wakeUpSyncWithOtherProcesses(f4options,reduceOptions);
	doAutoReduce=reduceOptions.doAutoReduce;
	int& matrixSize = reduceOptions.matrixSize;
	int& usefulProcesses=reduceOptions.usefulProcesses;//Число процессов, которые имеет смысл использовать при данном размере матрицы
	if (matrixSize<200) usefulProcesses=min(1, f4options->numberOfProcs);//на маленьких матрицах лучше, чтоб вообще работал только 1 процесс
	else if (matrixSize<2000) usefulProcesses=min(2, f4options->numberOfProcs);//на средних 2 процесса
	else usefulProcesses=f4options->numberOfProcs;//иначе все
	ResIdCalculator resIdCalc(usefulProcesses);//отвечает за порядок циркулюции процессов
	CMatrix mymat;//Содержит строки из matrix, относящиеся к текущему процессору
	CMatrix resultmat;//Содержит готовые строки, относящиеся к текущему процессору
	vector<int> resultblocksizes;//размеры блоков строк, добавленных в resultmat
	mymat.reserve(1+matrixSize/usefulProcesses);//выделим сразу всю память, чтоб потом не копировать при push_back.
	resultmat.reserve(mymat.capacity());
	resultblocksizes.reserve(1+mymat.capacity()/f4options->MPIBlockSize);
	MPI_Comm commUseful;
	MPI_Group groupWorld, groupUseful;
	//Необходимость создания коммуникатора не по всем процессам
	bool needNonTrivialCommUseful = usefulProcesses!=f4options->numberOfProcs && usefulProcesses!=1;
	if (needNonTrivialCommUseful){
		vector<int> usefulRanks(usefulProcesses);
		for (int i=0;i<usefulRanks.size();++i){
			usefulRanks[i]=i;
		}
		MPI_Comm_group(MPI_COMM_WORLD, &groupWorld); 
		MPI_Group_incl(groupWorld, usefulRanks.size(), &(usefulRanks[0]), &groupUseful);
		//MPI_Group_excl(groupWorld, usefulRanks.size(), &(usefulRanks[0]), &groupUseful);  /* local */ 
		MPI_Comm_create(MPI_COMM_WORLD, groupUseful, &commUseful);
	}else{
		commUseful=MPI_COMM_WORLD;
	}
	if (ID<usefulProcesses){ //на этом процессоре нужно проводить рассчёт
		if (ID==0){
			//Выборка и рассылка матриц на все процессоры
			for (int id=1;id<usefulProcesses;++id){
				selectRowsForProcessor(mymat,id,matrixSize,usefulProcesses,f4options);
				sendSubMatrix(id,mymat.begin(),mymat.end(),f4options->MPIUseBigSends);
			}
			//Выборка своей матрицы
			selectRowsForProcessor(mymat,0,matrixSize,usefulProcesses,f4options);
		}else{
			recvToMatrix(0,mymat,f4options->MPIUseBigSends);
		}
		matrix.clear();
		ttt[1]=getMPITimeMesurement();
		reduceOptions.matrixSize=forwardGaussElimination(mymat,resultmat,f4options,reduceOptions,resultblocksizes,resIdCalc,commUseful);
		ttt[2]=getMPITimeMesurement();
		if (doAutoReduce){
			//Обратный ход совмещённый со сбором результата
			backGaussElimination(resultmat,matrix,f4options,reduceOptions,resultblocksizes,resIdCalc,commUseful);
		}else{
			//сбор результата
			if (ID){
				// Рассылаем матрицы 0му процессу
				// Строки в полученной матрице идут в произвольном порядке
				sendSubMatrix(0,resultmat.begin(),resultmat.end(),f4options->MPIUseBigSends);
			}else{
				matrix.reserve(reduceOptions.matrixSize);
				for (int id=1;id<usefulProcesses;++id){
					recvToMatrix(id,matrix,f4options->MPIUseBigSends);
				}
				int lastSize=matrix.size();
				matrix.resize(reduceOptions.matrixSize);
				iterator curline=matrix.begin()+lastSize;
				for (iterator i=resultmat.begin();i!=resultmat.end();++i,++curline){
					curline->swap(*i);
				}
			}
		}
		ttt[3]=getMPITimeMesurement();
	}
	if (needNonTrivialCommUseful){
		MPI_Group_free(&groupUseful); 
		MPI_Group_free(&groupWorld); 
		if (commUseful!=MPI_COMM_NULL) MPI_Comm_free(&commUseful);
	}
	/*if (!ID) {
		matrix.PrintMatrixInfo("reduced");
		fprintf(stats,"Detailed times:\n");
		const char *infos[]={
			"Initial bcast",
			"Forward Gauss",
			"Backward Gauss/Final bcast"
		};
		static const unsigned NUMBER_OF_TIMES=sizeof(ttt)/sizeof(ttt[0]);
		vector <MPITimeMesurement> othersTimes(f4options->numberOfProcs*NUMBER_OF_TIMES);
		MPI_Gather(ttt,sizeof(ttt),MPI_BYTE,&othersTimes.front(),sizeof(ttt),MPI_BYTE,0,MPI_COMM_WORLD);
		for (unsigned i=0;i+1<NUMBER_OF_TIMES;++i){
			fprintf(stats,"\t%s: ",infos[i]);
			fprintf(stats,"Absolute time=%f, cpu times in processes=[%f",ttt[i+1].mw-ttt[i].mw,(ttt[i+1].cl-ttt[i].cl)/CLOCKS_PER_SEC);
			for (unsigned id=1;id<nProc;++id){
				fprintf(stats," %f",(othersTimes[NUMBER_OF_TIMES*id+i+1].cl-othersTimes[NUMBER_OF_TIMES*id+i].cl)/CLOCKS_PER_SEC);
			}
			fprintf(stats,"]\n");
		}
	}else{
		MPI_Gather(ttt,sizeof(ttt),MPI_BYTE,0,0,MPI_BYTE,0,MPI_COMM_WORLD);
	}*/
}




void CMatrix::toDiagonalNormalForm(const F4AlgData* f4options){
	MPIDiagonalForm(f4options,true);
}


void CMatrix::toRowEchelonForm(const F4AlgData* f4options){
	MPIDiagonalForm(f4options,false);
}


void CMatrix::reduceRowByRow(Row& row,const Row& by, CModular mb){
	CModular c = row.getCoefByMonom(by.HM());
	if (c!=0){
		row.addRowMultypliedBy(by, -c*mb);
	}
}


void CMatrix::eraseEmptyRows(unsigned firstRowToSearch){
	CMatrix& m=*this;
	unsigned i=firstRowToSearch;
	while(i<m.size()){//Удалим пустые строки
		if (m[i].empty()){
			m[i].swap(m.back());
			m.pop_back();
			continue;
		}else ++i;
	}
}


void CRow::addRowMultypliedBy(const CRow& rowFrom, CModular multBy){
	CRow& rowTo=*this;
	CRow result;

	//выделим память под максимально возможное число элементов
	result.resize(rowTo.size()+rowFrom.size());

	
	const_iterator it1 = rowTo.begin();
	const_iterator it1Finish = rowTo.end();
	const_iterator it2 = rowFrom.begin();
	const_iterator it2Finish = rowFrom.end();		
	iterator itResult = result.begin();

	if(multBy!=CModular(1)){
		//Сложить с домножением
		while(it1!=it1Finish && it2!=it2Finish){
			if(it1->column==it2->column){				
				itResult->value = it1->value+multBy*it2->value;
				itResult->column = it1->column;
				++it1;
				++it2;
				if(itResult->value!=0){
					//result.push_back(resPair);
					++itResult;
				}
			}
			else if(it1->column<it2->column){
				*itResult=*it1;
				++it1;
				++itResult;
			}
			else{
				itResult->value = multBy*it2->value;				
				itResult->column = it2->column;
				++it2;
				++itResult;
			}
		}
		while(it1!=it1Finish){
			*itResult=*it1;
			++it1;
			++itResult;
		}
		while(it2!=it2Finish){
			itResult->value = multBy*it2->value;				
			itResult->column = it2->column;
			++it2;
			++itResult;
		}
	}else{
		//Просто сложить
		while(it1!=it1Finish && it2!=it2Finish){
			if(it1->column==it2->column){				
				itResult->value = it1->value+it2->value;
				itResult->column = it1->column;
				++it1;
				++it2;
				if(itResult->value!=0){
					//result.push_back(resPair);
					++itResult;
				}
			}
			else if(it1->column<it2->column){
				*itResult=*it1;
				++it1;
				++itResult;
			}else{
				*itResult=*it2;
				++it2;
				++itResult;
			}
		}
		while(it1!=it1Finish){
			*itResult=*it1;
			++it1;
			++itResult;
		}
		while(it2!=it2Finish){
			*itResult=*it2;
			++it2;
			++itResult;
		}
	}
	//установим size в фактически занятый размер. Реально память при этом не освобождается 
	result.resize(itResult-result.begin());
	rowTo.swap(result);
}	


/*
struct CMatrix::RowComparator{//оператор сравнения строк для выбора наиболее эффективной для редукции
	struct rowinfo{
		static const int VERY_BIG=2000000000;
		int size;
		int HMon;
		void getfrom(Row& r){
			size=r.size();
			HMon=r.HM();
		}
		void settoworst(){
			size=VERY_BIG;
			HMon=VERY_BIG;
		}
		bool isworst(){
			return size==VERY_BIG;
		}
	};
	struct rowinfowithprocess:rowinfo{
		int processor;
	};
	
	struct rowinfowithrownum{
		rowinfo info;
		int rownum;
		bool operator<(const rowinfowithrownum& r2)const{
			return comparator(info,r2.info)>0;
		}
	};
	

	static int comparator(const rowinfo& a,const rowinfo& b){
		//Функция должна быть строго антикоммутативна (для MPI)
		//поэтому если элементы несравнимы нужно возвращать 0
		
		//Сравнение числа элементов
		if (f4AlgOptions.useSizesForSelectingRow){
			int countres=a.size-b.size;
			//static const int ifl=2;
			//Если разница хоть сколько-нибудь заметна, то ориентируемся на неё
			//if (countres>ifl || countres<-ifl) return countres;
			//Оператор должен быть ассоциативен
			if (countres) return countres;
			
			//Иначе смотрим на положение ведущего элемента в строке
			if (a.HMon==b.HMon) return 0;
			return (b.HMon>a.HMon)?1:-1;
		}else return b.HMon-a.HMon;
	}
	
	//Функция должна быть строго коммутативна (для MPI)
	//результат записывается в b
	static void mpiopcomparator(void *a, void *b, int * ,MPI_Datatype*){
		//Длина пока игнорируется, сравниваются и обрабатываются только первые элементы
		rowinfowithprocess& ra=*(rowinfowithprocess*)a;
		rowinfowithprocess& rb=*(rowinfowithprocess*)b;
		int res=comparator(ra,rb);
		if (
				rb.isworst() ||(
					!(ra.isworst()) 
					&&(
						res<0
						||(res==0 && (ra.processor > (rb.processor)))
					)
				)
			){//если a лучше b
			rb=ra;
		}
	}
	bool operator()(Row& a,Row& b)const{
		//if (a.empty() || b.empty()) cout<<"Bug: empty line found"<<endl;
		rowinfo ra,rb;
		ra.getfrom(a);
		rb.getfrom(b);
		return comparator(ra,rb)>0;
	}
};*/


void CMatrix::backReduceRangeByRange(MatrixIterator mfrom, MatrixIterator mto, MatrixIterator byfrom, MatrixIterator byto){
	for (MatrixIterator i=byfrom;i!=byto;++i){
		CModular mb = CModular::inverseMod(i->HC());
		for (MatrixIterator j=mfrom;j!=mto;++j){
			reduceRowByRow(*j,*i,mb);
		}
	}
}


void CMatrix::fullAutoReduce(MatrixIterator from, MatrixIterator to){//вызывается только для маленьких матриц размера MPIblocksize
	//прямой ход
	for(MatrixIterator j=from;j!=to;++j){
		if (j->empty()) continue;
		CModular mb = CModular::inverseMod(j->HC());
		for(MatrixIterator i=j+1;i!=to;++i){
			if (i->empty()) continue;
			reduceRowByRow(*i,*j,mb);
		}
	}
	//обратный ход
	for(MatrixIterator i = to; i!=from; --i){
		MatrixIterator curi=i-1;
		if (curi->empty()) continue;
		CModular mb = CModular::inverseMod(curi->HC());
		for(MatrixIterator j = curi; j!=from; --j){
			MatrixIterator curj=j-1;
			if (curj->empty()) continue;
			reduceRowByRow(*curj,*curi,mb);
		}
	}
}


void CMatrix::fullAutoReduce(CMatrix& m){//вызывается только для маленьких матриц размера MPIblocksize
	//прямой ход
	for(unsigned j=0;j<m.size();++j){
		CModular mb = CModular::inverseMod(m[j].HC());
		for(unsigned i=j+1;i<m.size();++i){
			reduceRowByRow(m[i],m[j],mb);
		}
		m.eraseEmptyRows(j+1);//Удалим пустые строки, начиная с j+1
	}
	//обратный ход
	for(int i = m.size()-1; i>0; i--){
		CModular mb = CModular::inverseMod(m[i].HC());
		for(int j = i-1; j>=0; j--){
			reduceRowByRow(m[j],m[i],mb);
		}
	}
}


void CMatrix::reduceRangeByRow(MatrixIterator mfrom, MatrixIterator mto,const Row& by){
	CModular mb = CModular::inverseMod(by.HC());
	for(MatrixIterator m=mfrom;m!=mto;++m){//Редуцируем все строки
		if (m->empty()) continue;
		reduceRowByRow(*m,by,mb);
	}
}


MatrixInfo& CMatrix::getMyStats(const F4AlgData* f4options){
	return f4options->stats->matInfo.back();
}



void CMatrix::makeMatrixStats(SingleMatrixInfo & info){
	info.columns=getMonomialMap().size();
	info.rows=this->size();
	info.elems=0;
	for(int i = 0;i<info.rows; ++i){
		info.elems+=getRow(i).size();
	}

}


void CMatrix::doMatrixStatsPre(const F4AlgData* f4options){
	if (!f4options->detailedMatrixInfo) return;
	f4options->stats->matInfo.push_back(MatrixInfo());
	makeMatrixStats(getMyStats(f4options).mBefore);
	if (f4options->stats->matrixInfoFile){
		MatrixInfo& myInfo=getMyStats(f4options);
		char buf[1024];
		//sprintf - форматная строка должна гарантировать фиксированную длинну!
		sprintf(
			buf,
			"%3d: %6d x%6d (%8d, %5.2f%%)->",
			int(f4options->stats->matInfo.size()),
			myInfo.mBefore.rows,
			myInfo.mBefore.columns,
			myInfo.mBefore.elems,
			myInfo.mBefore.filling()*100
		);
		(*f4options->stats->matrixInfoFile)<<buf<<flush;
	}
	//getMyStats(f4options).beginTime=getTime();
}


void CMatrix::doMatrixStatsPost(const F4AlgData* f4options){
	if (!f4options->detailedMatrixInfo) return;
	makeMatrixStats(getMyStats(f4options).mAfter);
	//getMyStats(f4options).totalTime=getTimePassed(getMyStats(f4options).beginTime);
	if (f4options->stats->matrixInfoFile){
		MatrixInfo& myInfo=getMyStats(f4options);
		char buf[1024];
		//sprintf - форматная строка должна гарантировать фиксированную длинну!
		sprintf(
			buf,
			"%6d x%6d (%8d, %5.2f%%)     in %6.2f sec\n",
			myInfo.mAfter.rows,
			myInfo.mAfter.columns,
			myInfo.mAfter.elems,
			myInfo.mAfter.filling()*100,
			0./*myInfo.totalTime*/
		);
		(*f4options->stats->matrixInfoFile)<<buf<<flush;
	}
}
} //namespace F4MPI
