// GARLI version 0.96b8 source code
// Copyright 2005-2008 Derrick J. Zwickl
// email: zwickl@nescent.org
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//	NOTE: Portions of this source adapted from GAML source, written by Paul O. Lewis

#define PROGRAM_NAME "GARLI"
#define MAJOR_VERSION 0.96
#define MINOR_VERSION "1b"
//DON'T mess with the following 2 lines!.  They are auto substituted by svn.
#define SVN_REV "$Rev$"
#define SVN_DATE "$Date$"
 
//allocation monitoring stuff from Paul, Mark and Dave
#define WRITE_MEM_REPORT_TO_FILE
#define INSTANTIATE_MEMCHK

#ifdef WIN32
#include <conio.h>
#endif

#include "defs.h"
#include "population.h"
#include "individual.h"
#include "adaptation.h"
#include "sequencedata.h"
#include "garlireader.h"

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

#ifdef CUDA_GPU
#include "cudaman.h"
CudaManager *cudaman;
int cuda_device_number=0;
#endif

OutputManager outman;
bool interactive;

//This is annoying, but the substituted rev and date from svn are in crappy format.
//Get what we need from them
//revision string looks like this: $Rev$
std::string GetSvnRev(){
	string temp = SVN_REV;
	string ret;
	for(int i=0;i<temp.length();i++){
		char c = temp[i];
		if(isdigit(c)) {
			ret += c;
			}
		}
	return ret;
	}
//date string looks like this: $Date$
std::string GetSvnDate(){
	string temp = SVN_DATE;
	string ret;
	int i=0;
	int len = temp.length();
	while(i < len && temp[i] != ',') i++;
	i++;
	while(i < len && !isdigit(temp[i])) i++;
	while(i < len && temp[i] != ')') ret += temp[i++];
	return ret;
	}

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

void UsageMessage(char *execName){
#ifdef SUBROUTINE_GARLI
	outman.UserMessage("This MPI version is for doing a large number of search replicates or bootstrap");
	outman.UserMessage("replicates, each using the SAME config file.  The results will be exactly");
 	outman.UserMessage("identical to those obtained by executing the config file a comparable number");
 	outman.UserMessage("of times with the serial version of the program.");
	outman.UserMessage("\nUsage: The syntax for launching MPI jobs varies between systems");
	outman.UserMessage("Most likely it will look something like the following:");
	outman.UserMessage("  mpirun [MPI OPTIONS] %s -[# of times to execute config file]", execName);
	outman.UserMessage("Specifying the number of times to execute the config file is mandatory.");
	outman.UserMessage("This version will expect a config file named \"garli.conf\".");
	outman.UserMessage("Consult your cluster documentation for details on running MPI jobs\n");
#elif defined (OLD_SUBROUTINE_GARLI)
	outman.UserMessage("This MPI version is for doing a large number of independent jobs in batch, each");
	outman.UserMessage("using a DIFFERENT config file.  This might be useful for analyzing a large");
	outman.UserMessage("number of simulated datasets or for analyzing a single dataset under a variety");
	outman.UserMessage("of models or search settings.  The results will be exactly the same as if each");
	outman.UserMessage("config file were executed separately by a serial version of GARLI.");
	outman.UserMessage("\nUsage: The syntax for launching MPI jobs varies between systems");
	outman.UserMessage("Most likely it will look something like the following:");
	outman.UserMessage("  mpirun [MPI OPTIONS] %s [# of provided config files]", execName);
	outman.UserMessage("This version will expect config files named \"run0.conf\", \"run1.conf\", etc.");
	outman.UserMessage("Consult your cluster documentation for details on running MPI jobs\n");
#else
	outman.UserMessage("Usage: %s [OPTION] [config filename]", execName);
	outman.UserMessage("Options:");
	outman.UserMessage("  -i, --interactive	interactive mode (allow and/or expect user feedback)");
	if(interactive) outman.UserMessage("        (interactive is the default for the version you are running)");
	outman.UserMessage("  -b, --batch		batch mode (do not expect user input)");
	if(!interactive) outman.UserMessage("        (batch is the default for the version you are running)");
	outman.UserMessage("  -v, --version		print version information and exit");
	outman.UserMessage("  -h, --help		print this help and exit");
	outman.UserMessage("  -t			run internal tests (requires dataset and config file)");
#ifdef CUDA_GPU
	outman.UserMessage("  --device d_number	use specified CUDA device");
#endif
	outman.UserMessage("NOTE: If no config filename is passed on the command line the program\n   will look in the current directory for a file named \"garli.conf\"\n");
#endif
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

#elif defined( SUBROUTINE_GARLI ) || defined(OLD_SUBROUTINE_GARLI)
int SubGarliMain(int rank)
	{
	int argc=1;
	char **argv=NULL;
	//clear out whatever is in the reader already - it might be full if a single
	//process has called SubGarliMain multiple times
	GarliReader &reader = GarliReader::GetInstance();
	reader.ClearContent();
#else
int main( int argc, char* argv[] )	{
#endif

	CREATE_MEMCHK{//memory leak detecting trick - no overhead when turned off

	#ifdef MPI_VERSION
		MPIMain(argc, argv);
	#endif

	string conf_name;

	string svnRev = GetSvnRev();
	string svnDate = GetSvnDate();

#ifdef OLD_SUBROUTINE_GARLI
	char name[12];
	sprintf(name, "run%d.conf", rank);
	conf_name = name;
#elif defined(SUBROUTINE_GARLI)
	conf_name = "garli.conf";
#else
	conf_name = "garli.conf";
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
					if(!_stricmp(argv[curarg], "-b") || !_stricmp(argv[curarg], "--batch")) interactive=false;
#ifndef MAC_FRONTEND
					else if(!_stricmp(argv[curarg], "-i") || !_stricmp(argv[curarg], "--interactive")) interactive=true;
#else
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
#endif
					else if(argv[curarg][1]=='t') runTests = true;
					else if(!_stricmp(argv[curarg], "-v") || !_stricmp(argv[curarg], "--version")){
						outman.UserMessage("%s Version %.2f%s.r%s", PROGRAM_NAME, MAJOR_VERSION, MINOR_VERSION, svnRev.c_str());
#ifdef SUBROUTINE_GARLI
						outman.UserMessage("MPI run distributing version");
#endif
#ifdef OPEN_MP
						outman.UserMessage("OpenMP multithreaded version");
#endif
						outman.UserMessage("Copyright Derrick J. Zwickl 2005-2009");
						outman.UserMessage("zwickl@nescent.org");
						exit(0);
						}
					else if(!_stricmp(argv[curarg], "-h") || !_stricmp(argv[curarg], "--help")){
						UsageMessage(argv[0]);
						exit(0);
						}
#ifdef CUDA_GPU
					else if(!_stricmp(argv[curarg], "--device")) cuda_device_number = atoi(argv[++curarg]);
#endif
					else {
						outman.UserMessage("Unknown command line option %s", argv[curarg]);
						UsageMessage(argv[0]);
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

#ifdef SUBROUTINE_GARLI
			//override the ofprefix here, tacking .runXX onto it
			char temp[10];
			if(rank < 10) sprintf(temp, ".run0%d", rank);
			else sprintf(temp, ".run%d", rank);
			conf.ofprefix += temp;
#endif

			// now set the random seed
			int randomSeed;
			if(conf.randseed < 1){
				srand((unsigned)time(NULL));
				randomSeed = RandomInt(1, 100000);
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

			outman.UserMessage("Running BOINC %s Version %.2f%s.r%s (%s)\n", PROGRAM_NAME, MAJOR_VERSION, MINOR_VERSION, svnRev.c_str(), svnDate.c_str());
			if(confOK && conf.restart == true) outman.UserMessage("Found BOINC checkpoint files.  Restarting....\n");

			boinc_resolve_filename(datafile.c_str(), buffer, 2048);
			datafile = buffer;
#else	//not BOINC
			if(confOK == true)
				sprintf(temp_buf, "%s.screen.log", conf.ofprefix.c_str());
			else
				sprintf(temp_buf, "ERROR.log");

			if(conf.restart) outman.SetLogFileForAppend(temp_buf);
			else outman.SetLogFile(temp_buf);

			outman.UserMessage("Running %s Version %.2f%s.r%s (%s)\n", PROGRAM_NAME, MAJOR_VERSION, MINOR_VERSION, svnRev.c_str(), svnDate.c_str());

#endif

#ifdef SUBROUTINE_GARLI //MPI versions
			outman.UserMessage("->MPI Parallel Version<-\nNote: this version divides a number of independent runs across processors.");
			outman.UserMessage("It is not the multipopulation parallel Garli algorithm.\n(but is generally a better use of resources)");
#endif
#if defined(OPEN_MP)
			outman.UserMessage("->OpenMP multithreaded version for multiple processors/cores<-");
#elif !defined(SUBROUTINE_GARLI)
			outman.UserMessage("->Single processor version<-\n", svnRev.c_str(), svnDate.c_str());
#endif

#ifdef SINGLE_PRECISION_FLOATS
			outman.UserMessage("->Single precision floating point version<-\n");
#endif

#ifdef CUDA_GPU
			outman.UserMessage("->CUDA GPU version<-\n");
#endif

			outman.UserMessage("This version has undergone much testing, but is still a BETA VERSION.\n   - Please check results carefully! -");

			outman.UserMessageNoCR("Compiled %s %s", __DATE__, __TIME__);

#if defined (_MSC_VER)
			outman.UserMessage(" using Microsoft C++ compiler version %.2f", _MSC_VER/100.0);
#elif defined(__INTEL_COMPILER)
			outman.UserMessage(" using Intel icc compiler version %.2f", __INTEL_COMPILER/100.0);
#elif defined(__GNUC__)
			outman.UserMessage(" using GNU gcc compiler version %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#else
			outman.UserMessage("");
#endif

#ifdef NCL_NAME_AND_VERSION
			outman.UserMessage("Using %s\n", NCL_NAME_AND_VERSION);
#endif

			outman.UserMessage("Reading config file %s", conf_name.c_str());
			if(confOK == false) throw ErrorException("Error in config file...aborting");

#ifdef SUBROUTINE_GARLI
			if(conf.randseed != -1)
				throw ErrorException("You cannot specify a random number seed with the MPI version.  This would cause all of the\n\tindependent MPI processes to give exactly identical results.  Set randomseed to -1");
//DEBUG
//			if(conf.checkpoint || conf.restart)
//				throw  ErrorException("Checkpointing cannot be used with the MPI version.  Set \"writecheckpoints\" and \"restart\" to 0");
#endif

			//set up the model specification
			modSpec.SetupModSpec(conf);

			// Create the data object
			if(modSpec.IsAminoAcid() && modSpec.IsCodonAminoAcid() == false)
				data = new AminoacidData();
			else //all data besides AA will be read into a DNA matrix and
				//then converted if necessary
				data = new NucleotideData();

			GarliReader &reader = GarliReader::GetInstance();
			pop.usedNCL = reader.ReadData(datafile.c_str(), modSpec);
			if(! pop.usedNCL) throw ErrorException("There was a problem reading the data file.");
			const NxsCharactersBlock *charblock = reader.CheckBlocksAndGetCorrectCharblock(modSpec);

			//using 2 argument form of CreateMatrix (taken from paritition branch) but passing dummy charset for
			//now.  Exsets will be taken care of in CreateMatrix
			NxsUnsignedSet charSet;
			data->CreateMatrixFromNCL(charblock, charSet);

			if(modSpec.IsCodon()){
				CodonData *d = new CodonData(dynamic_cast<NucleotideData *>(data), modSpec.geneticCode);
				pop.rawData = data;
				data = d;
				//this probably shouldn't go here, but...
				if(modSpec.IsF1x4StateFrequencies()) d->SetF1X4Freqs();
				else if(modSpec.IsF3x4StateFrequencies()) d->SetF3X4Freqs();
				else if(modSpec.IsEmpiricalStateFrequencies()) d->SetCodonTableFreqs();
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
			outman.UserMessage(" %d uninformative variable characters.", data->NVarUninform());
			int total = data->NConstant() + data->NInformative() + data->NVarUninform();
			if(data->NMissing() > 0){
				outman.UserMessage(" %d characters were completely missing or ambiguous (removed).", data->NMissing());
				outman.UserMessage(" %d total characters (%d before removing empty columns).", total, data->GapsIncludedNChar());
				}
			else outman.UserMessage(" %d total characters.", total);

			outman.flush();

			data->Collapse();
			outman.UserMessage("%d unique patterns in compressed data matrix.\n", data->NChar());

#ifdef CUDA_GPU
			cudaman = new CudaManager(modSpec.nstates, modSpec.numRateCats, data->NChar(),
					cuda_device_number, CUDA_TEST_ITERATIONS, CUDA_PRINT_TESTS, CUDA_PRINT_DEVICE_QUERY);
#endif

			//DJZ 1/11/07 do this here now, so bootstrapped weights aren't accidentally stored as orig
			data->ReserveOriginalCounts();

			data->DetermineConstantSites();

			pop.Setup(&conf, data, 1, 0);
			pop.SetOutputDetails();

			if(runTests){
				outman.UserMessage("starting internal tests...");
				pop.RunTests();
				outman.UserMessage("******Successfully completed tests.******");
				return 0;
				}

			if(conf.runmode != 0){
				if(conf.runmode == 1)
					pop.ApplyNSwaps(10);
				else if(conf.runmode == 7)
					pop.VariableStartingTreeOptimization(false);
				else if(conf.runmode == 9)
					pop.VariableStartingTreeOptimization(true);
				else if(conf.runmode == 8){
#ifdef OPEN_MP
					throw ErrorException("can't estimate site rates in openmp version!");
#endif
#ifndef ALLOW_SINGLE_SITE
					throw ErrorException("the program must be compiled with ALLOW_SINGLE_SITE defined in defs.h to use site rate estimation (runmode = 8)!");
#endif
					pop.OptimizeSiteRates();
					}
				else if(conf.runmode > 20){
					pop.GenerateTreesOnly(conf.runmode);
					}
				else if(conf.runmode > 1)
					pop.SwapToCompletion(conf.startOptPrec);
				}
			else{
				//if no checkpoint files are actually found conf->restart will be set to zero
				if(pop.conf->restart) pop.conf->restart = pop.ReadStateFiles();

				pop.SetOutputDetails();
				if(pop.conf->bootstrapReps == 0){//NOT bootstrapping
					pop.PerformSearch();
					}
				else pop.Bootstrap();
				pop.FinalizeOutputStreams(2);
				}
			}catch(ErrorException err){
				if(outman.IsLogSet() == false){
					outman.SetLogFile("ERROR.log");
					if(interactive == false) UsageMessage(argv[0]);
					}
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
				if(interactive==true){
					outman.UserMessage("\n-Press enter to close program.-");
					char d=getchar();
					}
				return 1;
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

