#ifndef F4main_h
#define F4main_h
/**
\file
алгоритм F4, работающий с внутреннним представоением
*/

#include "settings.h"
#include "types.h"
///Содержит реализацию билиотеки libf4mpi
namespace F4MPI{
PolynomSet F4(PolynomSet F, const F4AlgData* f4options);
} //namespace F4MPI
#endif
