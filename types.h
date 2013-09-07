#ifndef Types_h
#define Types_h
/** \file
Определение основных типов для реализации F4()
*/

#include <vector>
#include <set>
#include <deque>
#include <utility>

#include "cmodular.h"
#include "cpolynomial.h"
#include "cmatrix.h"

namespace F4MPI{
///представление множества полиномов
typedef std::vector<CPolynomial> PolynomSet;
typedef std::deque<CPolynomial> SortedReducersSet;
///представление S-пар полиномов
typedef std::pair<CPolynomial, CPolynomial> SPair;
///представление множества S-пар
typedef std::vector<SPair> SPairSet;
typedef CMatrix::Row MatrixRow;
struct ReduceBySet;
} //namespace F4MPI
#endif
