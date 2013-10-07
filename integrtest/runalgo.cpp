#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <iterator>
using namespace std;
#include "mpi_start_info.h"
#include "libf4mpi.h"

struct ProgramOptions:F4AlgOptions{
	//Автоматическое именование результата
	int autoNameSuffix;
};
ProgramOptions localAlgOptions;

struct CMDLineOption{
	const char* cmdline;
	const char* fname;
	const char* helpcomment;
	int ProgramOptions::* value;
	enum cmdoptkinds {cmdopt_bool,cmdopt_int} kind;
};

CMDLineOption cmdlineoptions[]={
	{"--diag","DIAG" ,"use diagonal form (instead of rowechelon) each step", &ProgramOptions::diagonalEachStep, CMDLineOption::cmdopt_bool },
	{"--block","BLKS", "lines in reduced block", &ProgramOptions::innerGaussBlockSize, CMDLineOption::cmdopt_int},
	{"--MPIblock","MPIB", "lines in reducer block",  &ProgramOptions::MPIBlockSize, CMDLineOption::cmdopt_int},
	{"--MPIbig","MBIG", "minimize number of sends", &ProgramOptions::MPIUseBigSends, CMDLineOption::cmdopt_bool},
	{"--rowsz","SROW", "select reducing row by size", &ProgramOptions::useSizesForSelectingRow, CMDLineOption::cmdopt_bool},
	{"--time",0, "profile time", &ProgramOptions::profileTime, CMDLineOption::cmdopt_bool},
//	{"--shedul","SHED","use sheduler to select next reducer processor", &CMatrix::matrixSheduler, CMDLineOption::cmdopt_bool},
//	{"--circul","CIRC","method of process circulation (0-3)", &MPI_PROCESS_CIRCULATE_ORDER, CMDLineOption::cmdopt_int},
	{"--matrinfo",0, "generate detailed matrix info", &ProgramOptions::detailedMatrixInfo, CMDLineOption::cmdopt_bool},
	{"--autoname",0, "generate suffix for output filename", &ProgramOptions::autoNameSuffix, CMDLineOption::cmdopt_bool},
	{"--latex",0, "generate latex log", &ProgramOptions::generateLatexLog, CMDLineOption::cmdopt_bool},
	{"--algo",0, "Select algo. 0 = F4MPI, 1=F5, 2=F5C", &ProgramOptions::selectedAlgo, CMDLineOption::cmdopt_int}
};



void printUsage(const char* name){
	fprintf(stderr, "Usage: %s inputfile outputfile [OPTIONS] \n", name);
	fprintf(stderr, "OPTIONS (--option1 value1 --option2 value2 ... --optionN valueN):\n");
	for (const auto& cmdlineoption: cmdlineoptions){
		ostringstream hlp;
		hlp<<"  ";
		hlp.width(14);
		hlp<<left;
		hlp<<cmdlineoption.cmdline<<cmdlineoption.helpcomment<<"  Default: "<<localAlgOptions.*cmdlineoption.value;
		if (cmdlineoption.kind==CMDLineOption::cmdopt_bool) hlp<<", bool";
		hlp<<"\n";
		fprintf(stderr, "%s", hlp.str().c_str());
	}
}

int main(int argc, char * argv[]){
	MPIStartInfo mpi_info(argc, argv);
	initDefaultF4Options(&localAlgOptions);
	localAlgOptions.autoNameSuffix=0;
	localAlgOptions.showInfoToStdout=1;
	localAlgOptions.profileTime=1; //this is the default value for testf4mpi
	try{
		string outputname;
		string inputname;
		if (argc < 2){
			if (mpi_info.isMainProcess()){
				printUsage(argv[0]);
				fprintf(stderr, "Using input.dat by default...\n");
			}
			inputname="input.dat";
		}else{
			inputname=argv[1];
			string basename;
			string::size_type slashp=inputname.find_last_of("/\\");
			if (slashp==inputname.npos) basename=inputname;
			else{
				basename=inputname.substr(slashp+1);
			}
			if (basename.size()>4 && basename.substr(basename.size()-4)==".dat") basename.resize(basename.size()-4);
			int optstart=3;
			if (argc>optstart){
				int ci=optstart;
				while (ci<argc){
					string optname=argv[ci];
					const CMDLineOption* cur_option = nullptr;
					for	(const auto& cmdlineoption : cmdlineoptions){
						if (cmdlineoption.cmdline==optname) cur_option=&cmdlineoption;
					}
					
					if (!cur_option){
						fprintf(stderr, "\nERROR: Unknown option \"%s\"\n\n",argv[ci]);
						printUsage(argv[0]);
					}
					int *option=&(localAlgOptions.*(cur_option->value));
					++ci;
					char *err;
					if (ci<argc) *option=strtol(argv[ci],&err,0);

					if (ci>=argc || *err || (cur_option->kind==CMDLineOption::cmdopt_bool && (*option & ~1))){
						fprintf(stderr, "\nERROR: Option value expected after \"%s\"\n\n",argv[ci-1]);
						printUsage(argv[0]);
					}else ++ci;
				}
			}
			ostringstream optionsstream;
			if (!mpi_info.isSingleProcess()) optionsstream<<".mpi"<<mpi_info.numberOfProcs;
			for	(const auto& cmdlineoption : cmdlineoptions){
				if (!cmdlineoption.fname) continue;
				optionsstream<<'.'<<cmdlineoption.fname<<'_'<<localAlgOptions.*cmdlineoption.value;
			}
			if (argc>=3){
				outputname=argv[2];
				if (localAlgOptions.autoNameSuffix){
					outputname+=optionsstream.str();
				}
			}
		}
		int runningResult=runF4MPIFromFile(inputname.c_str(),outputname.c_str(),&localAlgOptions, mpi_info);
		if (runningResult<0 && mpi_info.isMainProcess()){
			const char* failReasons []={
				"NO error",
				"Parse failed",
				"Input open failed",
				"Output open failed"
			};
			fprintf(stderr,"runF4FromFile failed, reason: %s\n",failReasons[-runningResult]);
		}
	}catch(const std::exception& e){
		fprintf(stderr,"Exception: %s\n",e.what());
	}
	return 0;
}
