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
#include "sequencedata.h"

#include "funcs.h"
#include "mpifuncs.h"
#include "hashdefines.h"
#include "tree.h"
#include "errorexception.h"
#include "outputman.h"

#ifdef MAC_FRONTEND
#import <Foundation/Foundation.h>
#import "MFEInterfaceClient.h"
#endif

#ifdef MPI_VERSION
#include "mpi.h"
#endif

OutputManager outman;
bool interactive;

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

#ifdef BOINC
int boinc_garli_main( int argc, char* argv[] );

int main( int argc, char* argv[] ){
	int retval = boinc_init();
	if(retval){
		cout << "Problem initializing BOINC system!" << endl;
		return retval;
		}
	retval = boinc_garli_main( argc, argv );
	if(retval) boinc_finish(retval);
	else boinc_finish(0);
	}

int boinc_garli_main( int argc, char* argv[] )	{
	outman.SetNoOutput(true);

#elif defined( SUBROUTINE_GARLI )
int SubGarliMain(int rank)	{
	int argc=1;
	char **argv=NULL;
#else
int main( int argc, char* argv[] )	{
#endif

	CREATE_MEMCHK{//memory lead detecting trick - no overhead when turned off

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
	interactive=false;
#else
	interactive=true;
#endif

	bool runTests = false;
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
					else if(argv[curarg][1]=='t') runTests = true;
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


		//create the population object
		Population pop;
		SequenceData *data = NULL;
		
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
#ifdef OPEN_MP
			outman.UserMessage("OpenMP multithreaded version for multiple processors/cores"); 
#endif
			outman.UserMessage("Compiled %s %s\n", __DATE__, __TIME__); 

			outman.UserMessage("Reading config file %s", conf_name.c_str());
			if(confOK == false) throw ErrorException("Error in config file...aborting");

			//set up the model specification
			modSpec.SetupModSpec(conf);

			// Create the data object
			if(modSpec.IsAminoAcid() && modSpec.IsCodonAminoAcid() == false)
				data = new AminoacidData();
			else //all data besides AA will be read into a DNA matrix and
				//then converted if necessary
				data = new NucleotideData();

			pop.usedNCL = ReadData(datafile.c_str(), data);
			
			if(modSpec.IsCodon()){
				CodonData *d = new CodonData(dynamic_cast<NucleotideData *>(data), modSpec.geneticCode);
				pop.rawData = data;
				data = d;
				//this probably shouldn't go here, but...
				if(modSpec.IsF1x4StateFrequencies()) d->SetEmpType(1);
				if(modSpec.IsF3x4StateFrequencies()) d->SetEmpType(2);
				}
			else if(modSpec.IsCodonAminoAcid()){
				AminoacidData *d = new AminoacidData(dynamic_cast<NucleotideData *>(data), modSpec.geneticCode);
				pop.rawData = data;
				data = d;				
				}

			data->Summarize();
			outman.UserMessage("\nSummary of dataset:");
			outman.UserMessage(" %d sequences.", data->NTax());
			outman.UserMessage(" %d constant characters.", data->NConstant());
			outman.UserMessage(" %d parsimony-informative characters.", data->NInformative());
			outman.UserMessage(" %d autapomorphic characters.", data->NAutapomorphic());
			int total = data->NConstant() + data->NInformative() + data->NAutapomorphic();
			outman.UserMessage(" %d total characters.", total);

			if(data->NMissing() > 0)
				outman.UserMessage(" %d characters were completely missing or ambiguous (removed).", data->NMissing());

			outman.flush();

			data->Collapse();
			outman.UserMessage("%d columns in compressed data matrix.\n", data->NChar());

			//DJZ 1/11/07 do this here now, so bootstrapped weights aren't accidentally stored as orig
			data->ReserveOriginalCounts();

			if(modSpec.includeInvariantSites && modSpec.IsCodon() == false)
				data->DetermineConstantSites();

			pop.Setup(&conf, data, 1, 0);

			if(runTests){
				outman.UserMessage("starting internal tests...");
				pop.RunTests();
				throw ErrorException("(not actually an error!) Successfully completed tests.");
				}
			
			if(conf.runmode != 0){
				if(conf.runmode == 1)
					pop.ApplyNSwaps(10);
				if(conf.runmode == 7)
					pop.VariableStartingTreeOptimization();
				else if(conf.runmode > 1)
					pop.SwapToCompletion(conf.startOptPrec);
				}
			else{
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
				pop.FinalizeOutputStreams(0);
				pop.FinalizeOutputStreams(1);
				pop.FinalizeOutputStreams(2);

#ifdef MAC_FRONTEND
				NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
				NSString *messageForInterface = [NSString stringWithUTF8String:err.message];
				[[MFEInterfaceClient sharedClient] didEncounterError:messageForInterface];
				[pool release];
#endif				

#ifdef BOINC
				return 1;
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

	return 0;
};

