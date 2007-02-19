// GARLI version 0.951 source code
// Copyright  2005-2006 by Derrick J. Zwickl
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//  Derrick Zwickl
//	National Evolutionary Synthesis Center
//	2024 W. Main Street, Suite A200
//	Durham, NC 27705
//  email: zwickl@nescent.org
//
//	NOTE: Portions of this source adapted from GAML source, written by Paul O. Lewis


//allocation monitoring stuff from Paul, Mark and Dave
#define WRITE_MEM_REPORT_TO_FILE
#undef MONITORING_ALLOCATION
#ifdef MONITORING_ALLOCATION
	#define INSTANTIATE_MEMCHK
#endif
#include "defs.h"

#ifdef WIN32
#include <conio.h>
#endif

#include "population.h"
#include "individual.h"
#include "adaptation.h"
#include "mlhky.h"

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

#ifdef MAC_FRONTEND
#import <Foundation/Foundation.h>
#import "MFEInterfaceClient.h"
#endif

#ifdef PROFILE
#include "profiler.h"
#endif

#ifdef MPI_VERSION
#include "mpi.h"
#endif

char programName[81];

OutputManager outman;

int CheckRestartNumber(const string str){
	int num=1;
	char temp_buf[100];
	sprintf(temp_buf, "%s.restart%d.best.tre", str.c_str(), num);
	while(FileExists(temp_buf)) {
		num ++;
		sprintf(temp_buf, "%s.restart%d.best.tre", str.c_str(), num);
		}
	return num;
	}

#ifndef SUBROUTINE_GARLI
int main( int argc, char* argv[] )	{
#else
int SubGarliMain(int rank)	{
int argc=1;
char **argv=NULL;
#endif

	#ifdef SIOUX
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

	#ifdef MPI_VERSION
		MPIMain(argc, argv);
	#else
	// init some stuff

	string conf_name;
#ifndef SUBROUTINE_GARLI
	conf_name = "garli.conf";
#else
	char temp[100];
	sprintf(temp, "run%d.conf", rank);
	conf_name = temp;
#endif

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
#ifdef MAC_FRONTEND
					else if (argv[curarg][1] == 'i') {
						NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
						NSString *arg = [[[NSProcessInfo processInfo] arguments] objectAtIndex:curarg];
						NSString *serverName = [[arg componentsSeparatedByString:@"="] objectAtIndex:1];
						BOOL success = [[MFEInterfaceClient sharedClient] connectToServerWithName:serverName];
						if (!success) {
							NSLog(@"Failed to connect to interface server");
							EXIT_FAILURE;
						}
						[pool release];
					}
#endif															
					else {
						outman.UserMessage("Unknown command line option %s", argv[curarg]);
						exit(0);
						}
					}
				//if anything else appears, we'll assume that it's a config file
				else conf_name = argv[curarg];
	        curarg++;
			}
		}
#ifdef GANESH
	if(Tree::random_p==false) Tree::ComputeRealCatalan();
#endif

		try{
			MasterGamlConfig conf;
			bool confOK;
			confOK = ((conf.Read(conf_name.c_str()) < 0) == false);

			// now set the random seed
			int randomSeed;
			if(conf.randseed < 1){
				srand((unsigned)time(NULL));
				randomSeed = RandomInt(1, 10000);
				}
			else randomSeed=conf.randseed;
			rnd.set_seed(randomSeed);
			
			char temp_buf[100];
			if(confOK == true){
				if(conf.restart == false) sprintf(temp_buf, "%s.screen.log", conf.ofprefix.c_str());
				else sprintf(temp_buf, "%s.restart%d.screen.log", conf.ofprefix.c_str(), CheckRestartNumber(conf.ofprefix));
				}
			else sprintf(temp_buf, "ERROR.log");
			outman.SetLogFile(temp_buf);
			
			outman.UserMessage("Running serial GARLI, version 0.951 (Oct 2006)\n");
			outman.UserMessage("Reading config file %s", conf_name.c_str());
			if(confOK == false) throw ErrorException("Error in config file...aborting");

			// Create the data object
			HKYData data;
			int err=ReadData(conf.datafname.c_str(), &data);
			if(err==-1) return 0;
			
			//create the population object
			Population pop;
		
			pop.Setup(&conf, &data, 1, 0);

			#ifdef PROFILE
				ProfilerInit(collectDetailed, bestTimeBase, 100000, 1000);
			#endif

			if(pop.conf->bootstrapReps == 0){
				pop.GetConstraints();
				if(conf.restart == false){
					outman.UserMessage("Running Genetic Algorithm with initial seed=%d\n", rnd.init_seed());
					pop.SeedPopulationWithStartingTree();
					//DEBUG - to look at effect of prec during init opt on score
/*					pop.InitializeOutputStreams();
					time_t repStart;
					ofstream res("optresults.log");
					for(FLOAT_TYPE prec=0.5;prec > 0.0001;){
						repStart = pop.stopwatch.SplitTime();
						for(int rep=0;rep<10;rep++){
							pop.adap->branchOptPrecision = prec;
							pop.SeedPopulationWithStartingTree();
							pop.AppendTreeToTreeLog(-1, -1);
							res << prec << "\t" << pop.indiv[0].Fitness() << endl;
							}
						res << "TIME: " << prec << "\t" << pop.stopwatch.SplitTime() - repStart << endl;
						if(prec > 0.051) prec -= 0.05;
						else prec -= 0.01;
						}
					return 1;
*/					}
				else{
					pop.ReadStateFiles();
					outman.UserMessage("Restarting Genetic Algorithm from checkpoint");
					outman.UserMessage("generation %d, seed %d, best lnL %.3f", pop.gen, rnd.init_seed(), pop.BestFitness());
					pop.adap->SetChangeableVariablesFromConfAfterReadingCheckpoint(&conf);
					}
				pop.InitializeOutputStreams();
				pop.AppendTreeToTreeLog(-1, -1);
				pop.Run();
				}
			else{
				pop.InitializeOutputStreams();
				if(conf.restart == true) throw(ErrorException("Restarting of bootstrap runs is not supported.\nYou should simply start a new bootstrap run and\ncombine all trees obtained into one bootstrap sample."));
				if(conf.inferInternalStateProbs == true) throw(ErrorException("You cannont infer internal states during a bootstrap run!"));
				pop.OutputModelReport();
				//if there are not mutable params in the model, remove any weight assigned to the model
				if(pop.indiv[0].mod->NumMutatableParams() == 0) {
					outman.UserMessage("NOTE: Model contains no mutable parameters!\nSetting model mutation weight to zero.\n");
					pop.adap->modelMutateProb=0.0;
					pop.adap->UpdateProbs();
					}
				pop.GetConstraints();
				pop.Bootstrap();
				}
			}catch(ErrorException err){
				outman.UserMessage("\nERROR: %s\n\n", err.message);
#ifdef MAC_FRONTEND
				NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
				NSString *messageForInterface = [NSString stringWithUTF8String:err.message];
				[[MFEInterfaceClient sharedClient] didEncounterError:messageForInterface];
				[pool release];
#endif				
				}
			catch(int error){
				if(error==Population::nomem) cout << "not able to allocate enough memory!!!" << endl;
				}

	if(interactive==true){
		outman.UserMessage("\n-Press enter to close program.-");
		char d=getchar();
		}
	
	outman.CloseLogFile();

#endif
	}
	#if defined(MONITORING_ALLOCATION) && !defined(NDEBUG)
		#if defined(WRITE_MEM_REPORT_TO_FILE)
			char filename[50];
			#ifndef WIN32
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

