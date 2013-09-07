#ifndef MatrixInfoImpl_h
#define MatrixInfoImpl_h
/**\file
структуры, представляющие статистику по матрицам
*/
namespace F4MPI{
///Информация о матрицe в конкретный момент
struct SingleMatrixInfo{
	int columns;///<число столбцов
	int rows;///<число строк
	int elems;///<число ненулевых элементов
	///Возвращает заполненность - долю ненулевых элементов среди всех
	double filling(){
		double area=double(columns)*double(rows);
		if (area==0) return 1;
		else return double(elems)/double(area);
	}
};

///Информация о проведённом преобразовании одной матрицы
struct MatrixInfo{
	SingleMatrixInfo mBefore;///<Состояние матрицы до редукции
	SingleMatrixInfo mAfter;///<Состояние матрицы после редукции
	//AbsoluteTime beginTime;///<Время начала редукции
	//DifferenceTime totalTime;///<Время, затраченное ра редукцию
};
} //namespace F4MPI
#endif

