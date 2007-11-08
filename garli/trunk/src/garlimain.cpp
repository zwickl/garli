// GARLI version 0.96b4 source code
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
#undef INSTANTIATE_MEMCHK
#include "defs.h"
#include "memchk.h"
#define WRITE_MEM_REPORT_TO_FILE

#ifdef MONITORING_ALLOCATION
	#define INSTANTIATE_MEMCHK
#endif

#ifdef WIN32
#include <conio.h>
#endif

#include "defs.h"
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

#ifdef MPI_VERSION
#include "mpi.h"
#endif

#undef RUN_TESTS

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
	
#ifdef BOINC
	int retval = boinc_init();
	if(retval){
		cout << "Problem initializing BOINC system!" << endl;
		return retval;
		}
	outman.SetNoOutput(true);

#endif

	CREATE_MEMCHK{

	#ifdef MPI_VERSION
		MPIMain(argc, argv);
	#endif

	string conf_name;
#ifndef SUBROUTINE_GARLI
	conf_name = "garli.conf";
#else
	char temp[100];
	sprintf(temp, "run%d.conf", rank);
	conf_name = temp;
#endif

#if defined(UNIX) || defined(BOINC)
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

		HKYData *data = NULL;
		//create the population object
		Population pop;

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

			string datafile = conf.datafname;

#ifdef BOINC
			//deal with stdout and stderr, although I don't think that anything is being
			//sent to them in BOINC mode
			char buffer[2048];
			boinc_resolve_filename("boinc_stdout", buffer, 2048);
			FILE *stdout_fp = freopen(buffer, "w", stdout);
			boinc_resolve_filename("boinc_stderr", buffer, 2048);
			FILE *stderr_fp = freopen(buffer, "w", stderr);

			//check for the presence of BOINC checkpoint files
			conf.restart = true;
			
			sprintf(temp_buf, "%s.adap.check", conf.ofprefix.c_str());
			boinc_resolve_filename(temp_buf, buffer, sizeof(buffer));
			if(FileExists(buffer) == false) conf.restart = false;

			sprintf(temp_buf, "%s.pop.check", conf.ofprefix.c_str());
			boinc_resolve_filename(temp_buf, buffer, sizeof(buffer));
			if(FileExists(buffer) == false) conf.restart = false;

			if(conf.uniqueSwapBias != ONE_POINT_ZERO){
				sprintf(temp_buf, "%s.swaps.check", conf.ofprefix.c_str());
				boinc_resolve_filename(temp_buf, buffer, sizeof(buffer));
				if(FileExists(buffer) == false) conf.restart = false;
				}

			if(confOK == true){
				sprintf(temp_buf, "%s.screen.log", conf.ofprefix.c_str());
				}
			else sprintf(temp_buf, "ERROR.log");

			if(conf.restart)
				outman.SetLogFileForAppend(temp_buf);
			else 
				outman.SetLogFile(temp_buf);

			outman.UserMessage("Running BOINC GARLI, version 0.96beta4 (Aug 2007)\n(->BOINC Amino acid and Codon testing version<-)\n");
			if(confOK && conf.restart == true) outman.UserMessage("Found BOINC checkpoint files.  Restarting....\n");

			boinc_resolve_filename(datafile.c_str(), buffer, 2048);
			datafile = buffer;
#else
			if(confOK == true){
				//changing this to always append to the .screen.log after a restart
				sprintf(temp_buf, "%s.screen.log", conf.ofprefix.c_str());
				//if(conf.restart == false) sprintf(temp_buf, "%s.screen.log", conf.ofprefix.c_str());
				//else sprintf(temp_buf, "%s.restart%d.screen.log", conf.ofprefix.c_str(), CheckRestartNumber(conf.ofprefix));
				}
			else sprintf(temp_buf, "ERROR.log");

			if(conf.restart) outman.SetLogFileForAppend(temp_buf);
			else outman.SetLogFile(temp_buf);

			outman.UserMessage("Running serial GARLI, version 0.96beta4 (Aug 2007)\n(->Amino acid and Codon testing version<-)\n");
#endif
			outman.UserMessage("Reading config file %s", conf_name.c_str());
			if(confOK == false) throw ErrorException("Error in config file...aborting");

			//set up the model specification
			modSpec.SetupModSpec(conf);

			// Create the data object
			if(modSpec.IsNucleotide() || modSpec.IsAminoAcid())
				data = new HKYData();
			else
				data = new CodonData();

			pop.usedNCL = ReadData(datafile.c_str(), data);
		
			pop.Setup(&conf, data, 1, 0);

			if(conf.runmode != 0){
				if(conf.runmode == 1)
					pop.ApplyNSwaps(10);
				else if(conf.runmode > 1)
					pop.SwapToCompletion(conf.startOptPrec);
				}
			else{
#ifdef RUN_TESTS
				pop.RunTests();
			
#endif
				if(pop.conf->restart) pop.ReadStateFiles();

				pop.SetOutputDetails();
				if(pop.conf->bootstrapReps == 0){//NOT bootstrapping
					pop.PerformSearch();
					}
				else pop.Bootstrap();
				pop.FinalizeOutputStreams(2);
				}
			}catch(ErrorException err){
				if(outman.IsLogSet() == false)
					outman.SetLogFile("ERROR.log");
				outman.UserMessage("\nERROR: %s\n\n", err.message);
#ifdef BOINC
				boinc_finish(1);
#endif

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
	if(data != NULL) delete data;
	
	if(interactive==true){
		outman.UserMessage("\n-Press enter to close program.-");
		char d=getchar();
		}

	outman.CloseLogFile();
#ifdef BOINC
	boinc_finish(0);
#endif

/* OLD WAY
			if(pop.conf->bootstrapReps == 0){
				pop.GetConstraints();
				if(conf.restart == false){
					outman.UserMessage("Running Genetic Algorithm with initial seed=%d\n", rnd.init_seed());
					pop.SeedPopulationWithStartingTree();
*/					//DEBUG - to look at effect of prec during init opt on score
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
/*				else{
					pop.ReadStateFiles();
					outman.UserMessage("Restarting Genetic Algorithm from checkpoint");
					outman.UserMessage("generation %d, seed %d, best lnL %.3f", pop.Gen(), rnd.init_seed(), pop.BestFitness());
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
/*
#ifdef BOINC
				boinc_finish(1);
#endif

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

	delete data;
	outman.CloseLogFile();
#ifdef BOINC
	boinc_finish(0);
#endif
	}
*/
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

#ifdef BOINC
	boinc_finish(0);
#endif

	return 0;
};

