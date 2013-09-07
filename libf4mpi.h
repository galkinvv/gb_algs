#ifndef libf4mpi_h
#define libf4mpi_h
/**
\file
Основной заголовочный файл библиотеки.
Определяет функции для вычисления базиса Грёбнера системы многочленов представленной в текстовом виде,
и структуру содержащую возможные параметры алгоритма.
*/

#include "f4_mpi_settings.h"

#ifdef __cplusplus
#include <string>
extern "C" {
#endif

/**установка параметров по умолчанию.
устанавливает набор параметров *opts в значения по умолчанию, близкие к оптимальным.
Перед передачей параметров в алгоритм следует инициализировать их этой функцией и после этого лишь откорректировать желаемые.
*/
void initDefaultF4Options(F4AlgOptions* opts);
#ifndef __cplusplus
LibF4ReturnCode runF4MPIFromFile(const char* inputName, const char* outputName, const F4AlgOptions* f4options);
#else

/**вычисление базиса из файла.
Считывает задачу из файла \a inputName, вычисляет базис с параметрами f4options и записывает полученный базис в файл \a outputName.
При указанни соотвествующих опций также в отдельные сохраняется статитстика.
\retval успешность выполнения алгоритма в соответсвии со значениями кодов возврата
*/
LibF4ReturnCode runF4MPIFromFile(const char* inputName, const char* outputName, const F4AlgOptions* f4options=0);
} //extern "C"

/**вычисление базиса из строки.
Считывает задачу из \a input, вычисляет базис с параметрами f4options и записывает полученный базис в строку \a output.
\retval успешность выполнения алгоритма в соответсвии со значениями кодов возврата
*/
LibF4ReturnCode runF4MPIFromString(const std::string& input, std::string& output, const F4AlgOptions* f4options=0);
#endif
#endif
