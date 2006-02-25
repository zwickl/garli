// GARLI version 0.94 source code
// Copyright  2005 by Derrick J. Zwickl
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//  Derrick Zwickl
//	Integrative Biology, UT
//	1 University Station, C0930
//	Austin, TX  78712
//  email: zwickl@mail.utexas.edu
//
//	Note: In 2006  moving to NESCENT (The National
//	Evolutionary Synthesis Center) for a postdoc

//	NOTE: Portions of this source adapted from GAML source, written by Paul O. Lewis


#ifdef WINDOWS
#include <conio.h>
#endif

#include "population.h"
#include "individual.h"
#include "parameters.h"
#include "adaptation.h"
#include "mlhky.h"
#include "defs.h"
#include "funcs.h"
#include "mpifuncs.h"
#include "hashdefines.h"
#include "tree.h"
#include "mlhky.h"
#include "errorexception.h"

//include Sioux console manipulators
#ifdef WINDOWS
#include "winsioux.h"
#else
#ifdef MAC
#include "sioux.h"
#endif
#endif

//allocation monitoring stuff from Paul, Mark and Dave
#define WRITE_MEM_REPORT_TO_FILE
#undef MONITORING_ALLOCATION
#ifdef MONITORING_ALLOCATION
	#define INSTANTIATE_MEMCHK
#endif
#include "memchk.h"

#ifdef PROFILE
#include "profiler.h"
#endif

#ifdef MPI_VERSION
#include "mpi.h"
#endif

char programName[81];

int main( int argc, char* argv[] )	{

	#ifdef MAC
	SIOUXSettings.columns = 90;
	SIOUXSettings.rows = 45;
	SIOUXSettings.asktosaveonclose = false;
	SIOUXSettings.autocloseonquit = true;
	SIOUXSettings.showstatusline = true;
	SIOUXSetTitle("\pGARLI 0.94");
	#endif


	bool poo=true;
//	while(poo);

	CREATE_MEMCHK{
	srand((unsigned)time(NULL));
	#ifdef MPI_VERSION
		MPIMain(argc, argv);
	#else
	// init some stuff

	char conf_name[20];
	strcpy(conf_name, "garli.conf");

    if (argc > 1) {
    	int curarg=1;
        while(curarg<argc){
#ifdef GANESH 
     		if(strlen(argv[curarg]) <= 2 && isdigit(argv[curarg][0])) {
	        	Tree::p_value       = atoi(argv[curarg]);
                if (Tree::p_value == 0) {
                    Tree::random_p = true;
                }
//                else {
  //                  Tree::ComputeRealCatalan();
    //            }
            }
	      	else
#endif
	        	strcpy(conf_name, argv[curarg]);
	        curarg++;
			}
		}
#ifdef GANESH
	if(Tree::random_p==false) Tree::ComputeRealCatalan();
#endif

		try{

			cout << "Running serial GARLI, version 0.94 (Feb 06)" << endl;
			MasterGamlConfig conf;
			cout << "Reading config file " << conf_name << endl;
			if (conf.Read(conf_name) < 0)	{
				cout << "Error in config file...aborting." << endl;
				return 0;
			}

			// Create the data object
			HKYData data;
			int err=ReadData(conf.datafname.c_str(), &data);
			if(err==-1) return 0;
			
			// create the parameters object
			Parameters params;
			params.SetParams(conf, data);
			
			Tree::meanBrlenMuts	= params.meanBrlenMuts;
			Tree::alpha		= params.gammaShapeBrlen;

			// create the population object
			Population pop;
		
			pop.Setup(params, &conf, 1, 0);

			#ifdef PROFILE
				ProfilerInit(collectDetailed, bestTimeBase, 100000, 1000);
			#endif
			if(pop.bootstrapReps == 0)
			pop.Run();
			else pop.Bootstrap();
			}catch(ErrorException err){
				err.Print(cout);
				return 0;
				}
			catch(int error){
				if(error==Population::nomem) cout << "not able to allocate enough memory!!!" << endl;
				return 0;
				}

#ifndef UNIX
	cout << "-Press enter to close program.-" << endl;
//	while( cin && (cin.get() != '\n') );
	char d=getchar();
#endif
	exit(0);


#endif
	}
	#if defined(MONITORING_ALLOCATION) && !defined(NDEBUG)
		#if defined(WRITE_MEM_REPORT_TO_FILE)
			char filename[50];
			#ifndef WIN32_VERSION
			int rank=0;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			sprintf(filename, "memecheck%d.txt", rank);
			#else
			strcpy(filename, "memcheck.txt");
			#endif
			ofstream memf(filename);
			MEMCHK_REPORT(memf)
			memf.close();
		#else
			MEMCHK_REPORT(cout)
		#endif
	#endif
	#ifdef PROFILE
	ProfilerDump("\phalfnewrescale.prof");
	ProfilerTerm();
	#endif

	return 0;
};

