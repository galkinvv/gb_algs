#pragma once

#if WITH_MPI
#include <mpi.h>
#endif

struct MPIStartInfo
{
	///Число процессоров в MPI
	int numberOfProcs;

	///ранг текущего процесса в MPI
	int thisProcessRank;

	MPIStartInfo(int &argc, char ** &argv)
	{
#if WITH_MPI
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numberOfProcs);
	MPI_Comm_rank(MPI_COMM_WORLD,&thisProcessRank);
#else
	numberOfProcs = 1;
	thisProcessRank = 0;
#endif
	}
	bool isMainProcess()const
	{
		return !thisProcessRank;
	}
	bool isSingleProcess()const
	{
		return numberOfProcs == 1;
	}
	~MPIStartInfo()
	{
#if WITH_MPI
		MPI_Finalize();
#endif
	}
};
