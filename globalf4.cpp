/**
\file
Размещает глобальные переменные библиотеки
*/

#include "globalf4.h"
#include "cmodular.h"
#include "parse.tab.h"

namespace globalF4MPI{
	MemoryManager MonomialAllocator;
//	PolynomMap globalPolynomMap;

	GlobalOptions globalOptions;
	//Инициализация глобальных переменных после определения их в парсере
	void InitializeGlobalOptions(){
		CMonomial::setOrder((CMonomial::Order)globalOptions.monomOrder, globalOptions.monomOrderParam);
		CModular::setMOD(globalOptions.mod);
		CMonomial::setNumberOfVariables(globalOptions.numberOfVariables);
		globalF4MPI::MonomialAllocator.setSize(CMonomial::degreessize);
		PODvecSize<CInternalMonomial>::setvalsize(CMonomial::degreessize);
		MonomialAllocator.reset();
	}
	void Finalize(){
		MonomialAllocator.reset();
	}
}

