/**
\file
Реализация внешнего интерфейса библиотеки.
Процедуры, подготавливающие полученные текстовые многочлены к передаче в F4()
и вызывающие инициализизаторы глобальных переменных во всех процессах MPI.
По завершении вывод статистики и вызов обратное преобразования в текст
*/
#if WITH_MPI
#include <mpi.h>
#endif
#include "gbimpl.h"
#include "outputroutines.h"
#include "matrixinfoimpl.h"
#include "settings.h"

#include "parse.tab.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cstring>
#include <cerrno>
#include <memory>
using namespace std;
namespace F4MPI{

/**Текстовые описания параметров F4.
Каждая переменная этого типа устанавливает соответствие между параметрами F4AlgData и их описаниями в файле статистики
*/
struct F4OptForStatsDescription{
	const char* statsname;///<текстовое описание для файла статистики
	int F4AlgData::* value;///<отдельный параметр (указатель на член-данное)
};

F4OptForStatsDescription optionsDescription[]={
//	{"Number of processes         ", &F4AlgData::numberOfProcs},
	{"Diagonal form each iteration", &F4AlgData::diagonalEachStep},
	{"Inner loop block size       ", &F4AlgData::innerGaussBlockSize},
	{"MPI block size              ", &F4AlgData::MPIBlockSize},
	{"MPI use big sends           ", &F4AlgData::MPIUseBigSends},
	{"Use sizes for selecting row ", &F4AlgData::useSizesForSelectingRow}
//	{"matrixSheduler", &CMatrix::matrixSheduler},
//	{"MPIProcessCirculation", &MPI_PROCESS_CIRCULATE_ORDER},
};

///Сообщает о коде возврата \a result, возникшем в процессе \a root всем остальным процессам, возвращая его
LibF4ReturnCode MPICheckResult(LibF4ReturnCode result=LIBF4_NO_ERROR, int root=0){
#if WITH_MPI
	MPI_Bcast(&result, sizeof(result), MPI_CHAR, root, MPI_COMM_WORLD);
#endif	
	return result;
}

/**вычисление базиса из потока.
Считывает задачу из \a input, вычисляет базис с параметрами f4givenOptions и записывает полученный базис в строку \a output.
При ненулевых параметрах статистика по матрицам собирается в \a f4stats, а по времени записывается в файл \a stats.
\retval успешность выполнения алгоритма в соответсвии со значениями кодов возврата
*/
LibF4ReturnCode runF4FromStream(istream& input, ostream& output, F4AlgOptions& f4givenOptions, const MPIStartInfo &mpi_start_info, F4Stats* f4stats=0, FILE* stats=0, ostream* latexLog = 0){
	F4Stats localf4stats;
	if (f4stats==0) f4stats=&localf4stats;
	F4AlgData f4data(f4givenOptions, f4stats, mpi_start_info, latexLog);
	PolynomSet givenSet;
	ParserVarNames varNames;
	LibF4ReturnCode parseSuccess=LIBF4_NO_ERROR;
	if (mpi_start_info.isMainProcess()){
		parseSuccess=LibF4ReturnCode(ParseInput(input, givenSet, &varNames));
		if(f4data.showInfoToStdout){
			if (parseSuccess>=0){
				string varDesc;
				for (int i=1;i<=globalF4MPI::globalOptions.numberOfVariables;++i)
				{
					if (i>1) varDesc += " > ";
					varDesc += CMonomialBase::varName(i, &varNames) + "=" +  CMonomialBase::varName(i, 0); 
				}
				printf("Task: order=%s, mod=%d, %d variables (%s)",
						getMonomialOrderName((CMonomial::Order)globalF4MPI::globalOptions.monomOrder).c_str(),
						globalF4MPI::globalOptions.mod,
						globalF4MPI::globalOptions.numberOfVariables,
						varDesc.c_str()
					);
				if (globalF4MPI::globalOptions.monomOrderParam){
					printf(", order_option=%d",
							globalF4MPI::globalOptions.monomOrderParam
						);
				}
				printf("\n");
				printf("Using %d processes\n",
						mpi_start_info.numberOfProcs
					);
			}
			fflush(stdout);
		}
	}
	//Разошлём всем успешность парсинга
	parseSuccess=MPICheckResult(parseSuccess);
	if (parseSuccess<0){
		if (mpi_start_info.isMainProcess()){
			globalF4MPI::Finalize();
		}
		return LIBF4_ERR_PARSE_FAILED;//неудачное завершение разбора
	}
#if WITH_MPI
	//Разошлём всем параметры, определённые в парсере
	MPI_Bcast(&globalF4MPI::globalOptions, sizeof(globalF4MPI::globalOptions), MPI_CHAR, 0, MPI_COMM_WORLD);
	if (!mpi_start_info.isMainProcess()){
		//во всех процессах, кроме главного нужно провести инициализацию полученных параметров
		globalF4MPI::InitializeGlobalOptions();
	}
#endif	
	if (stats){
		fprintf(stats, "\nOptions:\n");
		for (int i=0; i<sizeof(optionsDescription)/sizeof(optionsDescription[0]);++i){
			if (!optionsDescription[i].statsname) continue;
			fprintf(stats, "%s: %d\n", optionsDescription[i].statsname,f4data.*optionsDescription[i].value);
		}
	}

	PolynomSet basis;
	
	//Выполнение алгоритма F4
	if (mpi_start_info.isMainProcess()){
		basis = GB(givenSet, &f4data);
#if WITH_MPI
		int finished=true;
		MPI_Bcast(&finished,1,MPI_INT,0,MPI_COMM_WORLD);
	}else{
		CMatrix localmatrix;
		//Цикл вызовов вычислительной части (MPIDiagonalForm) в остальных процессах MPI
		for(;;){
			int finished=false;
			MPI_Bcast(&finished,1,MPI_INT,0,MPI_COMM_WORLD);
			if (finished) break;
			localmatrix.MPIDiagonalForm(&f4data);
		}
#endif
	}
	//Вывод результатов и сохранение статистики
	if (mpi_start_info.isMainProcess()){
		if (f4data.showInfoToStdout){
			printf("Computed basis - %d polynomials.\n", int(basis.size()));
			fflush(stdout);
		}
		PrintPolynomSet(output,basis,&varNames);
		output.flush();
		if (f4stats->matrixInfoFile){
			(*f4stats->matrixInfoFile)<<"Total matrices: "<<f4stats->matInfo.size()<<endl;
		}
		if (f4data.profileTime){
			if (stats){
				fflush(stats);
			}
		}
	}
	globalF4MPI::Finalize();
	return LIBF4_NO_ERROR;
}

/**Инициализирует параметры F4 во всех процессах.
\param givenF4Options параметры, переданные пользователем, возможно \c NULL (в таком случае используются значения по умолчанию).
\param localF4Options результирующие параметры, которые следует использовать (в основном процессе копируются переданные пользователем, остальные получают от него)
*/
void bcastF4Options(const F4AlgOptions* givenF4Options, F4AlgOptions& localF4Options, const MPIStartInfo &mpi_start_info){
	if (mpi_start_info.isMainProcess()){
		if (givenF4Options){
			localF4Options=*givenF4Options;
		}else{
			initDefaultF4Options(&localF4Options);
		}
	}
#if WITH_MPI
	MPI_Bcast(&localF4Options,sizeof(localF4Options),MPI_CHAR,0,MPI_COMM_WORLD);
#endif
}
} //namespace F4MPI
using namespace F4MPI;
LibF4ReturnCode runF4MPIFromString(const std::string& input, std::string& output, const F4AlgOptions* f4options, const MPIStartInfo &mpi_start_info){
	F4AlgOptions localF4Options;
	bcastF4Options(f4options, localF4Options, mpi_start_info);
	istringstream inputStream;
	ostringstream outputStream;
	inputStream.str(input);
	LibF4ReturnCode result=runF4FromStream(inputStream,outputStream,localF4Options,mpi_start_info);
	output=outputStream.str();
	return result;
}

LibF4ReturnCode runF4MPIFromFile(const char* inputName, const char* outputName, const F4AlgOptions* f4options, const MPIStartInfo &mpi_start_info){
	F4AlgOptions localF4Options;
	bcastF4Options(f4options, localF4Options, mpi_start_info);
	ifstream input;
	ofstream outputFile;
	ostream *outputPtr=&outputFile;
	ofstream matrixInfo;
	std::unique_ptr<ofstream> latexLog;
	FILE *stats=0;
	F4Stats f4stats;
	LibF4ReturnCode successCode;
	if (mpi_start_info.isMainProcess()){
		try{
			input.open(inputName);
			if (!input.is_open()){
				fprintf(stderr, "input (%s) open error: %s\n", inputName, strerror(errno));
				successCode=LIBF4_ERR_INPUT_OPEN_FAILED;
				throw runtime_error("opening input failed");
			}
			if(localF4Options.showInfoToStdout){
				printf("Opened file %s...\n",inputName);
				fflush(stdout);
			}
			if (outputName!=0 && strlen(outputName)==0) outputName=0;
			if (outputName != 0){
				outputFile.open(outputName);
				if (!outputFile.is_open()){
					fprintf(stderr, "output (%s) open error: %s\n", outputName, strerror(errno));
					successCode=LIBF4_ERR_OUTPUT_OPEN_FAILED;
					throw runtime_error("opening output failed");
				}
				string statsName=string(outputName)+".stats";
				if (localF4Options.profileTime){
					stats = fopen(statsName.c_str(), "w");
					if (!stats){
						//non-fatal error
						fprintf(stderr, "WARNING: stats (%s) open error: %s\n", statsName.c_str(), strerror(errno));
					}
				}
				if (localF4Options.detailedMatrixInfo){
					string matrInfoName=outputName;
					matrInfoName+=".matrinfo";
					matrixInfo.open(matrInfoName.c_str());
					if (!matrixInfo.is_open()){
						//non-fatal error
						fprintf(stderr, "WARNING: matrinfo (%s) open error: %s\n", matrInfoName.c_str(), strerror(errno));
					}
					f4stats.matrixInfoFile=&matrixInfo;
				}
				if (localF4Options.generateLatexLog){
					auto latex_log_name = string(outputName)+".tex";
					latexLog.reset(new ofstream(latex_log_name));
					if (!latexLog->is_open()){
						//non-fatal error
						fprintf(stderr, "WARNING: latex log (%s) open error: %s\n", latex_log_name.c_str(), strerror(errno));
					}
				}
				outputPtr=&outputFile;
			}else{
				outputPtr=&cout;
			}
			successCode=MPICheckResult(LIBF4_NO_ERROR);
		}catch(...){
			successCode=MPICheckResult(successCode);
		}
	}else{
		//оставшиеся процессы MPI получают успешность от нулевого
		successCode=MPICheckResult();
	}
	if (successCode!=0) return successCode;
	successCode=runF4FromStream(input,*outputPtr,localF4Options,mpi_start_info,&f4stats, stats,latexLog.get());
	if (stats) fclose(stats);
	return successCode;
}

void initDefaultF4Options(F4AlgOptions* opts){
	opts->detailedMatrixInfo=0;
	opts->diagonalEachStep=1;
	opts->useSizesForSelectingRow=0;
	opts->MPIUseBigSends=1;
	opts->autoReduceBasis=1;
	opts->MPIBlockSize=32;
	opts->innerGaussBlockSize=1;
	opts->profileTime=0;
	opts->minPercent=0;
	opts->showInfoToStdout=0;
	opts->generateLatexLog=0;
	opts->selectedAlgo=0;
}
