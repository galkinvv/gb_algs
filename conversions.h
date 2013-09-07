#pragma once
/**\file
Преобразование многочленов в матрицу и обратно.
Содуржит процедуры преобразований строка матрицы \<=\> многочлен и,
реализованные на их основе, преобразования матрица \<=\> множество многочленов
Процедуры объявлены шаблонными, чтоб не требовать всех нужных заголовочных файлов для типов аргументов,
которые в любом случае будут определены в месте вызова.
*/

namespace F4MPI{
/**полином -> строка матрицы
преобразует полином \a p в строку матрицы, используя соотвествие между мономами и номерами столбцов \a M,
и записывает её на место \a result
*/
template <class CPolynomial, class Row, class MonomialMap>
void polynomialToRow(const CPolynomial &p, const MonomialMap &M, Row& result){
	result.resize(p.size());
	typename Row::iterator coef=result.begin();
	for(int i = 0; i!=result.size(); ++i){	
		coef->value = p.getCoeff(i);
		//Это условие надо убрать, когда в полиномах перестанут попадаться нули
		if (coef->value!=0){
			//к этому моменту в M уже должны быть все мономы из матрицы(при вызове из PolyToMatrix)
			coef->column = M.getMonomialID(p.getMon(i));
			++coef;
		}
	}
	result.resize(coef-result.begin());
}		

/**строка матрицы -> полином
преобразует строку \a r матрицы в полином, используя соотвествие между мономами и номерами столбцов \a M,
и записывает его на место \a result
*/
template <class CPolynomial, class Row, class MonomialMap>
void rowToPolynomial(const Row& r,const MonomialMap &M, CPolynomial& result){
	result.resize(r.size());
	typename Row::const_iterator copyFrom = r.begin();
	typename Row::const_iterator copyFromEnd = r.end();
	typename CPolynomial::m_iterator mcopyTo=result.m_begin();
	typename CPolynomial::c_iterator ccopyTo=result.c_begin();
	for( ; copyFrom!=copyFromEnd; ++copyFrom, ++ccopyTo, ++mcopyTo){
		*mcopyTo=M.getMonomialRev(copyFrom->column);
		*ccopyTo=copyFrom->value;	
	}		
}	

/**матрица -> набор полиномов
преобразует строки матрицы \a m в набор полиномов \a polys,
используя соотвествие между мономами и номерами столбцов сохранённое в матрице.
Строки, ведущие элементв которых находятся в столбцах соответствующих мономам из \a ignoreLines пропускаются.
*/
template <class CMatrix, class PolynomialSet, class MonomialMap>
void matrixToPoly(const CMatrix& m, PolynomialSet& polys, const MonomialMap& ignoreLines){
	for(int i = 0; i<m.size(); ++i){
		const typename CMatrix::Row& R = m.getRow(i);
		if(R.empty())
			continue;
		if (ignoreLines.containsMonomial(m.getMonomialMap().getMonomialRev(R.front().column))) continue;
		polys.resize(polys.size()+1);//Добавить пустой многочлен в множество
		rowToPolynomial(R,m.getMonomialMap(),polys.back());//записать на его место преобразованную строку
	}
}

/**набор полиномов -> матрица
Строит соответсвие между мономами и номерами столбцов создаваемой матрицы
и сохраняет его в матрицы для использования в последствии при обратоном преобразовании.
Далее преобразует набор полиномов \a polys в строки матрицы \a m, используя это соотвествие.
*/
template <class CMatrix, class PolynomialSet>
void polyToMatrix(const PolynomialSet& polys, CMatrix& m){
	m.clear();
	m.resize(polys.size());
	//добавим в M все мономы, которые будут в матрице
	for(typename PolynomialSet::const_iterator i = polys.begin(); i!=polys.end(); ++i){
		storeMonomialsFromPoly(m.getMonomialMap(), *i);
	}
	//переупорядочим мономы в M по возростанию
	m.getMonomialMap().UpdateForUsingReversed();
	typename CMatrix::iterator row=m.begin();
	for(typename PolynomialSet::const_iterator i = polys.begin(); i!=polys.end(); ++i,++row){
		polynomialToRow(*i,m.getMonomialMap(),*row);
		if(row->empty()) --row;///
	}
	m.resize(row-m.begin());///
}
} //namespace F4MPI
