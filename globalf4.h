#ifndef globalF4_h
#define globalF4_h

#include "memorymanager.h"
/**\file
Объявляет глобальные переменные библиотеки
*/

///Глобальные данные алгоритма F4
namespace globalF4MPI{
	using namespace F4MPI;
	///аллокатор памяти для хранения мономов
	extern MemoryManager MonomialAllocator;
	
	/**основные глобальные пармаетры.
	Содержит глобальные параметры, определяемые в парсере и рассылаемые на все процессоры.
	*/
	struct GlobalOptions{
		int numberOfVariables;///<число переменных
		int mod;///<модуль по которому идут вычисления CModular
		int monomOrder;///<код порядка на мономах
		int monomOrderParam;///<параметр порядка на мономах
	};

	///центральная "точка доступа" к глобальным парметрам
	extern GlobalOptions globalOptions;

	/**Инициализация глобальных переменных.
	после определения в парсере или получения от главного процесса,
	глобальные опции нужно сообщить всем классам, поведение которых от них зависит.
	Особенно важно установить правильный размер памяти, занимаемый данными монома - это нужно сделать и как
	в аллокаторе MonomialAllocator для отдельных мономов, так и в классе PODvecSize\<CInternalMonomial\> который хранит масивы мономов в многочленах.
	*/

	void InitializeGlobalOptions();

	/**Деннициализация глобальных ресурсов.
	Освобождает всю память, выделенную аллокатором #MonomialAllocator.
	Это нужно делать отдельно, т.к. после удаления мономов память не освобождается, а складыывается в пул.
	*/
	void Finalize();

}
#endif
