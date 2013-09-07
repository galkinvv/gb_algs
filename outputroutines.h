/** \file
Прототипы процедур для преобразование в текст.
*/


#ifndef OutputRoutines_h
#define OutputRoutines_h

#include "types.h"
namespace F4MPI{
/**печать множества полиномов.
\param output поток, в который нужно произвести вывод
\param polys множество полиномов, которое надо напечатать
\param names соответствие между индексами и текстовыми именами переменных.
Если передан нулевой указатель (по умолчанию), используются стандартные имена x1 .. xN
*/
void PrintPolynomSet(std::ostream& output,PolynomSet& polys,ParserVarNames* names=0);

void PrintPolynomSetRaw(std::ostream& output,PolynomSet& polys,ParserVarNames* names=0);


/**печать множества полиномов, содержащихся в матрице.
\param output поток, в который нужно произвести вывод
\param matrix матрица, представляющая множество полиномов, которое надо напечатать
\param names соответствие между индексами и текстовыми именами переменных.
Если передан нулевой указатель (по умолчанию), используются стандартные имена x1 .. xN
*/
void PrintMatrixAsPolynomials(std::ostream& output,CMatrix& matrix, ParserVarNames* names=0);

/**печать множества S-пар.
\param output поток, в который нужно произвести вывод
\param pairs множество S-пар, которое надо напечатать
*/
void PrintSPairSet(std::ostream& output,SPairSet& pairs);
} //namespace F4MPI
#endif
