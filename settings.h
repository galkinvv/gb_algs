#ifndef settings_h
#define settings_h
#include "libf4mpi.h"
#include "mpi_start_info.h"
#include "matrixinfoimpl.h"

#include <iosfwd>
#include <vector>
/**\file
стуктуры статистики и параметров F4
*/
namespace F4MPI{
///статистика по матрицам
struct F4Stats{
	F4Stats():
		totalNumberOfReducedMatr(0),
		matrixInfoFile(0)
	{}

	///число редуцированных матриц
	int totalNumberOfReducedMatr;

	///статистика по редукции отдельных матриц
	std::vector<MatrixInfo> matInfo;

	///указатель на поток для записи статистики по матрицам, или NULL если сохранять не надо
	std::ostream* matrixInfoFile;
};

/**Внутренние параметры F4.
Расширенные внутренние параметры F4 состоят из переданных в библиотеку параметров #F4AlgOptions,
информации о MPI, место под статистику и передаются практически во все существенные шаги алгоритма.
*/
struct F4AlgData:F4AlgOptions{
	F4AlgData(const F4AlgOptions& options, F4Stats *f4stats, const MPIStartInfo &a_mpi_start_info, std::ostream* a_latex_log):
			F4AlgOptions(options),
			mpi_start_info(a_mpi_start_info),
			stats(f4stats),
			latex_log(a_latex_log)
	{}

	const MPIStartInfo &mpi_start_info;

	///указатель на место под статистику по матрицам
	F4Stats *stats;

	///stream for writing latex log
	std::ostream* latex_log;

};

} //namespace F4MPI
#endif

