#ifndef parse_tab_h
#define parse_tab_h
/**\file
Разбор задачи в текстовом виде
*/

#include "types.h"
#include <istream>
namespace F4MPI{
/**
Преобразует текстовую задачу в множество многочленов.
\param ins поток, содержащий текстовое представление задачи
\param readSet множество, в которое будет записано считанное множество
\param varNames указатель на переменную, в которую следует сохранить соответствие между номерами переменных в представлении монома и их текстовыми именами.
\retval успешность завершения. 0 - успешное, \<0 - произошла ошибка
*/
int ParseInput (std::istream& ins, PolynomSet& readSet, ParserVarNames* varNames);
} //namespace F4MPI
#endif

