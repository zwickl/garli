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


#ifdef WIN32
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
#include "outputman.h"

//include Sioux console manipulators
#ifdef MAC
	#include "sioux.h"
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

OutputManager outman;

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

#ifdef UNIX
	bool interactive=false;
#else
	bool interactive=true;
#endif

    if (argc > 1) {
    	int curarg=1;
        while(curarg<argc){
				if(argv[curarg][0]=='-'){
					//command line arguments with a dash
					if(argv[curarg][1]=='b') interactive=false;
					else {
						outman.UserMessage("Unknown command line option %s", argv[curarg]);
						exit(0);
						}
					}
				//if anything else appears, we'll assume that it's a config file
	        	else strcpy(conf_name, argv[curarg]);
	        curarg++;
			}
		}
#ifdef GANESH
	if(Tree::random_p==false) Tree::ComputeRealCatalan();
#endif

		try{
			MasterGamlConfig conf;
			bool confOK;
			confOK = ((conf.Read(conf_name) < 0) == false);
			//if(conf.Read(conf_name) < 0) throw ErrorException("Error in config file...aborting");
			
			char temp_buf[50];
			if(confOK == true) sprintf(temp_buf, "%s.screen.log", conf.ofprefix.c_str());
			else sprintf(temp_buf, "ERROR.log");
			outman.SetLogFile(temp_buf);
			
			outman.UserMessage("Running serial GARLI, version 0.95 BETA2 (July 2006)\nFlex Rates and Constraints testing version\n");
			outman.UserMessage("Reading config file %s", conf_name);
			if(confOK == false) throw ErrorException("Error in config file...aborting");

			// Create the data object
			HKYData data;
			int err=ReadData(conf.datafname.c_str(), &data);
			if(err==-1) return 0;
			
			// create the parameters object
			Parameters params;
			params.SetParams(conf, data);
			
			Tree::meanBrlenMuts	= params.meanBrlenMuts;
			Tree::alpha		= params.gammaShapeBrlen;

			//make sure that this is set before allocating any models	
			Model::noPinvInModel = conf.dontInferProportionInvariant;
			Model::nRateCats = conf.numratecats;
			Model::useFlexRates = conf.useflexrates;

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
//				err.Print(cout);
				outman.UserMessage("ERROR: %s\n\n", err.message);
				//err.Print(stderr);
				}
			catch(int error){
				if(error==Population::nomem) cout << "not able to allocate enough memory!!!" << endl;
				}

	if(interactive==true){
		outman.UserMessage("\n-Press enter to close program.-");
		char d=getchar();
		}
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

