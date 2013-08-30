// GARLI version 2.0 source code
// Copyright 2005-2011 Derrick J. Zwickl
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
#define MAJOR_VERSION "2"
#define MINOR_VERSION "01"
//DON'T mess with the following 2 lines!.  They are auto substituted by svn.
#define SVN_REV "$Rev$"
#define SVN_DATE "$Date$"

//allocation monitoring stuff from Paul, Mark and Dave
#define WRITE_MEM_REPORT_TO_FILE
#define INSTANTIATE_MEMCHK

#ifdef WIN32
#include <conio.h>
#endif

#ifdef MPI_VERSION
#include "mpi.h"
#endif

#include "defs.h"
#include "population.h"
#include "individual.h"
#include "adaptation.h"
#include "sequencedata.h"
#include "garlireader.h"

#include "funcs.h"
#include "tree.h"
#include "errorexception.h"
#include "outputman.h"

#ifdef WIN32
#include <process.h>
#define PID_FUNC() _getpid()
typedef int pid_type;
#else
#include <sys/types.h>
#include <unistd.h>
#define PID_FUNC() getpid()
typedef pid_t pid_type;
#endif

#ifdef MAC_FRONTEND
#import <Foundation/Foundation.h>
#import "MFEInterfaceClient.h"
#endif

#ifdef CUDA_GPU
#include "cudaman.h"
CudaManager *cudaman;
int cuda_device_number=0;
#endif

OutputManager outman;
bool interactive;
bool is64bit = false;
 
vector<ClaSpecifier> claSpecs;
vector<DataSubsetInfo> dataSubInfo;
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

void OutputVersion(){
	if(is64bit)
		outman.UserMessage("%s Version %s.%s.%s (64-bit)", PROGRAM_NAME, MAJOR_VERSION, MINOR_VERSION, GetSvnRev().c_str());
	else
		outman.UserMessage("%s Version %s.%s.%s (32-bit)", PROGRAM_NAME, MAJOR_VERSION, MINOR_VERSION, GetSvnRev().c_str());
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
	OutputVersion();
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
	OutputVersion();
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
	outman.UserMessage    ("Usage: %s [OPTION] [config filename]", execName);
	outman.UserMessage    ("Options:");
	outman.UserMessage                 ("  -i, --interactive	interactive mode (allow and/or expect user feedback)");
	if(interactive) 
		outman.UserMessage("                    (interactive is the default for the version you are running)");
	outman.UserMessage                 ("  -b, --batch		batch mode (do not expect user input)");
	if(!interactive) 
		outman.UserMessage("                    (batch is the default for the version you are running)");
	outman.UserMessage                 ("  -v, --version		print version information and exit");
	outman.UserMessage                 ("  -h, --help		print this help and exit");
	outman.UserMessage                 ("  -t			run internal tests (requires dataset and config file)");
	outman.UserMessage                 ("  -V			validate: load config file and data, validate config file, data, starting trees"); 
	outman.UserMessage                 ("				and constraint files, print required memory and selected model, then exit");
#ifdef CUDA_GPU
	outman.UserMessage    ("  --device d_number	use specified CUDA device");
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
	//ditto for these other globals
	modSpecSet.Delete();	
	claSpecs.clear();
	dataSubInfo.clear();
#else
int main( int argc, char* argv[] )	{
#endif

	CREATE_MEMCHK{//memory leak detecting trick - no overhead when turned off

	#ifdef MPI_VERSION
		MPIMain(argc, argv);
	#endif

	//I'm not sure that this is dependable or portable, but it is only for screen output, so isn't that important
	int ptrSize = sizeof(int *);
	if(ptrSize == 8)
		is64bit = true;
	else
		is64bit = false;

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
	bool validateMode = false;
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
					else if(!strcmp(argv[curarg], "-v") || !_stricmp(argv[curarg], "--version")){
						OutputVersion();
#ifdef SUBROUTINE_GARLI
						outman.UserMessage("MPI run distributing version");
#endif
#ifdef OPEN_MP
						outman.UserMessage("OpenMP multithreaded version");
#endif
#ifdef SINGLE_PRECISION_FLOATS
						outman.UserMessage("Single precision floating point version");
#endif
						outman.UserMessage("(DNA, AA, codon, morphology and partitioned models)");
						outman.UserMessage("Copyright Derrick J. Zwickl 2005-2011");
						outman.UserMessage("http://www.nescent.org/wg/garli/");
						outman.UserMessage("garli.support@gmail.com");
						exit(0);
						}
					else if(!_stricmp(argv[curarg], "-h") || !_stricmp(argv[curarg], "--help")){
						UsageMessage(argv[0]);
						exit(0);
						}

					else if(!strcmp(argv[curarg], "-V"))
						//validate mode skips some allocation in pop::Setup, and then executes pop::ValidateInput,
						//which is essentially a stripped down version of pop::SeedPopWithStartingTree
						validateMode = true;
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

		//population is defined here, but not allocated until much later
		Population *pop = NULL;

		DataPartition dataPart;
		DataPartition rawPart;
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
				//Add in the pid with the time to get the seed.  Otherwise forking a bunch
				//of runs simultaneously has a good chance of giving identical seeds
				//I believe unsigned overflow is guaranteed to wrap around safely
				pid_type pid = PID_FUNC();
				srand((unsigned)time(NULL) + (unsigned)pid);
				randomSeed = RandomInt(1, 1000000);
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
			
			outman.UserMessage("Running BOINC %s Version %s.%s.%s (%s)", PROGRAM_NAME, MAJOR_VERSION, MINOR_VERSION, svnRev.c_str(), svnDate.c_str());
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

			outman.UserMessage("Running %s Version %s.%s.%s (%s)", PROGRAM_NAME, MAJOR_VERSION, MINOR_VERSION, svnRev.c_str(), svnDate.c_str());

#endif

#ifdef SUBROUTINE_GARLI //MPI versions
			outman.UserMessage("->MPI Parallel Version<-\nNote: this version divides a number of independent runs across processors.");
			outman.UserMessage("It is not the multipopulation parallel Garli algorithm.\n(but is generally a better use of resources)"); 
#endif
#if defined(OPEN_MP)
			outman.UserMessageNoCR("->OpenMP multithreaded version for multiple processors/cores");
#elif !defined(SUBROUTINE_GARLI)
			outman.UserMessageNoCR("->Single processor version");
#endif

			if(is64bit)
				outman.UserMessage(" for 64-bit OS<-");
			else
				outman.UserMessage(" for 32-bit OS<-");

#ifdef SINGLE_PRECISION_FLOATS
			outman.UserMessage("->Single precision floating point version<-\n");
#endif

#ifdef CUDA_GPU
			outman.UserMessage("->CUDA GPU version<-\n");
#endif
			outman.UserMessage("##############################################################");
			outman.UserMessage(" This is GARLI 2.0, the first \"official\" release including ");
			outman.UserMessage("          partitioned models.  It is a merging of"); 
			outman.UserMessage("   official release 1.0 and beta version GARLI-PART 0.97");
			outman.UserMessage("  Briefly, it includes models for nucleotides, amino acids,");
			outman.UserMessage(" codons, and morphology-like characters, any of which can be ");
			outman.UserMessage("  mixed together and applied to different subsets of data.\n"); 
			outman.UserMessage("    General program usage is extensively documented here:");
			outman.UserMessage("            http://www.nescent.org/wg/garli/");
			outman.UserMessage("      see this page for details on partitioned usage:");
			outman.UserMessage("  http://www.nescent.org/wg/garli/Partition_testing_version");
			outman.UserMessage("   and this page for details on Mkv mophology model usage:");
			outman.UserMessage("    http://www.nescent.org/wg/garli/Mkv_morphology_model");
			outman.UserMessage("         PLEASE LET ME KNOW OF ANY PROBLEMS AT:");
			outman.UserMessage("                garli.support@gmail.com");
			outman.UserMessage("##############################################################");

			//outman.UserMessage("This version has undergone much testing, but is still a BETA VERSION.\n   - Please check results carefully! -");

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
			outman.UserMessage("Using %s", NCL_NAME_AND_VERSION);
#endif

			OutputImportantDefines();
			outman.UserMessage("\n#######################################################");
			outman.UserMessage("Reading config file %s", conf_name.c_str());
			if(confOK == false) throw ErrorException("Error in config file...aborting");

#ifdef SUBROUTINE_GARLI
			if(conf.randseed != -1)
				throw ErrorException("You cannot specify a random number seed with the MPI version.  This would cause all of the\n\tindependent MPI processes to give exactly identical results.  Set randomseed to -1");
#endif
			//This is pretty hacky.  Create one modSpec now because it is needed
			//to read the data (to identify the expected type of sequence for phylip and
			//fasta files), then add more later if there are multiple char blocks or CharPartitions
			//in the case of a Nexus datafile, we actually don't need to know the datatype
			//in advance, and can handle data subsets with different datatypes
			modSpecSet.AddModSpec(conf.configModelSets[0]);

			//read the datafile with the NCL-based GarliReader - should allow nexus, phylip and fasta
			outman.UserMessage("###################################################\nREADING OF DATA");
			GarliReader &reader = GarliReader::GetInstance();
			bool usedNCL = reader.ReadData(datafile.c_str(), *modSpecSet.GetModSpec(0));
			if(! usedNCL) 
				throw ErrorException("There was a problem reading the data file.");
			
			//assuming a single taxa block
			if(reader.GetNumTaxaBlocks() > 1) throw ErrorException("Expecting only one taxa block in datafile");
			NxsTaxaBlock *taxblock = reader.GetTaxaBlock(0);

			//currently data subsets will be created for each separate characters block, and/or for each
			//part of a char partition within a characters block
			int numCharBlocks = reader.GetNumCharactersBlocks(taxblock);
			if(numCharBlocks == 0) throw ErrorException("No character data (in characters/data blocks) found in datafile");
			vector<pair<NxsCharactersBlock *, NxsUnsignedSet> > effectiveMatrices;

			outman.UserMessage("\n###################################################\nPARTITIONING OF DATA AND MODELS");
			//loop over characters blocks
			for(int c = 0;c < numCharBlocks;c++){
				NxsCharactersBlock *charblock = reader.GetCharactersBlock(taxblock, c);
				string cbName = charblock->GetTitle();
				NxsAssumptionsBlock *assblock = NULL;
				NxsUnsignedSet charSet;
				bool foundCharPart = false;

				int numAssBlocks = reader.GetNumAssumptionsBlocks(charblock);
				if(numAssBlocks > 0){
					//loop over assumptions blocks for this charblock
					for(int a = 0;a < numAssBlocks;a++){
						assblock = reader.GetAssumptionsBlock(charblock, a);
						int numParts = assblock->GetNumCharPartitions();
						if(numParts > 1) 
							throw ErrorException("Found more than one CHARPARTITION referring to CHARACTERS block %s in a single ASSUMPTIONS or SETS blocks", charblock->GetTitle().c_str());
						else if(numParts == 1){
							if(foundCharPart == true)
								throw ErrorException("Found more than one CHARPARTITION referring to CHARACTERS block %s in multiple ASSUMPTIONS or SETS blocks", charblock->GetTitle().c_str());\
							else 
								foundCharPart = true;
							//get the name of the charpartition
							vector<std::string> charPartNames;
							assblock->GetCharPartitionNames(charPartNames);
							const NxsPartition *part = assblock->GetCharPartition(charPartNames[0]);
							int numSubsets = part->size();
							int subsetNum = 0;
							//loop over the partition subsets, each of which creates a data subset in GARLI
							for(NxsPartition::const_iterator subit = part->begin();subit != part->end();subit++){
								charSet = (*subit).second;
								dataSubInfo.push_back(DataSubsetInfo(effectiveMatrices.size(), c, cbName, subsetNum, (*subit).first, DataSubsetInfo::NUCLEOTIDE, DataSubsetInfo::NUCLEOTIDE));
								effectiveMatrices.push_back(make_pair(charblock, charSet));
								subsetNum++;
								}
							}
						}
					if(foundCharPart == false){//no charpart found 
						dataSubInfo.push_back(DataSubsetInfo(effectiveMatrices.size(), c, cbName, -1, "", DataSubsetInfo::NUCLEOTIDE, DataSubsetInfo::NUCLEOTIDE));
						effectiveMatrices.push_back(make_pair(charblock, charSet));
						}
					}
				else{ //no assumptions block found,  
					dataSubInfo.push_back(DataSubsetInfo(effectiveMatrices.size(), c, cbName, -1, "", DataSubsetInfo::NUCLEOTIDE, DataSubsetInfo::NUCLEOTIDE));
					effectiveMatrices.push_back(make_pair(charblock, charSet));
					}
				}

			//report on how data and models line up, and deal with a few unsupported possibilites
			if(conf.linkModels && conf.configModelSets.size() > 1)
				throw ErrorException("Multiple model subsets specified, but linkmodels = 1");
			if(effectiveMatrices.size() > 1){
				if(conf.configModelSets.size() == 1){//only one model description found
					if(conf.linkModels)
						outman.UserMessage("\nCHECK: ONE MODEL APPLIES TO ALL DATA SUBSETS\n\t(full linkage, all parameters shared)\n");
					else
						outman.UserMessage("\nCHECK: ONE MODEL TYPE APPLIES TO ALL DATA SUBSETS,\n\tBUT WITH INDEPENDENT MODEL PARAMETERS (no linkage)\n");
					}
				else{//mulitple model descriptions found
					if(conf.configModelSets.size() != effectiveMatrices.size())
						throw ErrorException("Multiple data subsets and model subsets specified, but numbers don't match");
					else outman.UserMessage("\nCHECK: DIFFERENT MODEL TYPES AND MODEL PARAMETERS APPLY\n\tTO EACH DATA SUBSET (no linkage)\n");
					}
				}
			else if(conf.configModelSets.size() != 1)
				throw ErrorException("Multiple models specified, but only one data subset found");

			//set this
			modSpecSet.SetInferSubsetRates(conf.subsetSpecificRates && effectiveMatrices.size() > 1);

			//now create a datamatrix object for each effective matrix
			//because of exsets some subsets of a charpart could contain no characters,
			//but I'm not going to deal with that right now, and will crap out

			//EFFECTIVE matrices ( = datasubsets) are the actual chunks of data specified by separate charblocks and/or charpartitions
			//There is a one to one matching between effective matrices and CLAs.  EXCEPT in the case of Nstate data,
			//in which case multiple matrices will be spawned 

			for(int dataChunk = 0;dataChunk < effectiveMatrices.size();dataChunk++){
				ModelSpecification *modSpec = NULL;
				//for Mk type data the number of actual matrices created can be > the number of actual data chunks
				//e.g., a first characters block might spawn 3 matrices for 2 state, 3 state and 4 state characters
				//sucessive char blocks then need to take that into account
				//a dataSubset is equivalent to a matrix in this respect
				int nextMatrixNum = dataPart.NumSubsets();

				if(dataChunk > 0){
					if(conf.linkModels){//linkage means that all clas/matrices point to the same model/modSpec
						//EXCEPT in the case of Mk type data with different numbers of states.  That is taken 
						//care of below.
						claSpecs.push_back(ClaSpecifier(dataChunk,0,dataChunk));
						modSpec = modSpecSet.GetModSpec(0);
						}
					else{//models are not linked ...
						if(conf.configModelSets.size() == 1)//but are all described by the same settings in the config file
							modSpecSet.AddModSpec(conf.configModelSets[0]);
						else{ //each has its own description in the config
							modSpecSet.AddModSpec(conf.configModelSets[dataChunk]);
							}
						claSpecs.push_back(ClaSpecifier(nextMatrixNum, nextMatrixNum, nextMatrixNum));
						modSpec = modSpecSet.GetModSpec(nextMatrixNum);
						}
					}
				else{ //if this is the first model, it must correspond to the first modSpec
					modSpec = modSpecSet.GetModSpec(0);
					claSpecs.push_back(ClaSpecifier(0,0,0));
					}
				if(conf.linkModels && (modSpec->IsMkTypeModel() || modSpec->IsOrientedGap()))
					throw ErrorException("Model linkage cannot be used with Mk/Mkv models (nor does it\n\tneed to be, since there are no estimated parameters).\n\tSet linkmodels = 0");

				//defaults here are NUCLEOTIDE, so make changes as necessary
				if(modSpec->IsCodon()) 
					dataSubInfo[dataChunk].usedAs = DataSubsetInfo::CODON;
				else if(modSpec->IsCodonAminoAcid())
					dataSubInfo[dataChunk].usedAs = DataSubsetInfo::AMINOACID;
				else if(modSpec->IsAminoAcid())
					dataSubInfo[dataChunk].readAs = dataSubInfo[dataChunk].usedAs = DataSubsetInfo::AMINOACID;
				else if(modSpec->IsNState())
					dataSubInfo[dataChunk].readAs = dataSubInfo[dataChunk].usedAs = DataSubsetInfo::NSTATE;
				else if(modSpec->IsNStateV())

					dataSubInfo[dataChunk].readAs = dataSubInfo[dataChunk].usedAs = DataSubsetInfo::NSTATEV;
				else if(modSpec->IsOrderedNState())
					dataSubInfo[dataChunk].readAs = dataSubInfo[dataChunk].usedAs = DataSubsetInfo::ORDNSTATE;
				else if(modSpec->IsOrderedNStateV())
					dataSubInfo[dataChunk].readAs = dataSubInfo[dataChunk].usedAs = DataSubsetInfo::ORDNSTATEV;
				else if(modSpec->IsOrientedGap())
					dataSubInfo[dataChunk].readAs = dataSubInfo[dataChunk].usedAs = DataSubsetInfo::ORIENTEDGAP;
				else if(modSpec->IsBinaryNotAllZeros())
					dataSubInfo[dataChunk].readAs = dataSubInfo[dataChunk].usedAs = DataSubsetInfo::BINARY_NOT_ALL_ZEROS;
				else if(modSpec->IsBinary())
					dataSubInfo[dataChunk].readAs = dataSubInfo[dataChunk].usedAs = DataSubsetInfo::BINARY;

				dataSubInfo[dataChunk].Report();
				//outman.UserMessage("");

				// Create the data object
				//for nstate data the effective matrices will be further broken up into implied matrices that each have the same number of observed states
				//the implied matrix number will be that number of states
				int actuallyUsedImpliedMatrixIndex = 0;
				int maxObservedStates = effectiveMatrices[dataChunk].first->GetMaxObsNumStates(false);
				//for Mk the impliedMatrix number is the number of states
				for(int impliedMatrix = 2;impliedMatrix < (modSpec->IsMkTypeModel() ? maxObservedStates + 1 : 3);impliedMatrix++){
					if(modSpec->IsMkTypeModel() && !modSpec->IsOrientedGap()){
						bool isOrdered = (modSpec->IsOrderedNState() || modSpec->IsOrderedNStateV());
						bool isBinary = modSpec->IsBinary() || modSpec->IsBinaryNotAllZeros();
						bool isConditioned =  (modSpec->IsNStateV() || modSpec->IsOrderedNStateV() || modSpec->IsBinaryNotAllZeros());
						//data = new NStateData(impliedMatrix, (modSpec->IsNStateV() || modSpec->IsOrderedNStateV()), (modSpec->IsOrderedNState() || modSpec->IsOrderedNStateV()));
						data = new NStateData(impliedMatrix, isOrdered, isBinary, isConditioned);
						}
					else if(modSpec->IsOrientedGap())
						data = new OrientedGapData();
					else if(modSpec->IsAminoAcid() && modSpec->IsCodonAminoAcid() == false)
						data = new AminoacidData();
					else //all other data will be read into a DNA matrix and
						//then converted if necessary
						data = new NucleotideData();

					//it really shouldn't be necessary to use the PatternManager on non-sequence data
					if(modSpec->IsNucleotide() || modSpec->IsAminoAcid() || modSpec->IsCodon())
						data->SetUsePatternManager(conf.usePatternManager);
					else
						data->SetUsePatternManager(0);
					
					//if no charpart was specified, the second argument here will be empty
					data->CreateMatrixFromNCL(effectiveMatrices[dataChunk].first, effectiveMatrices[dataChunk].second);

#ifdef SINGLE_PRECISION_FLOATS
					if(modSpec->IsMkTypeModel() || modSpec->IsOrientedGap()) throw ErrorException("Sorry, Mk/Mkv type models have not yet been tested with single precision.");
#endif
				
					if(data->NChar() == 0){
						//if there weren't any characters with a certain number of states,
						//just get rid of the matrix.  This could in theory also work for
						//totally excluded subsets, but that gets complicated because it
						//isn't clear how the indexing of models specified in the config
						//file should work
						assert(modSpec->IsMkTypeModel());
						outman.UserMessage("NOTE: No characters found with %d observed states.", impliedMatrix);
						delete data;
						}
					else{//now we have a data matrix object created, already filtered for the correct sites or number of states
						if(modSpec->IsMkTypeModel()){
							if(modSpec->IsGammaRateHet() || modSpec->IsFlexRateHet())
								throw ErrorException("Sorry, rate heterogeneity cannot be used with Mk/Mkv models yet.\n\tSet ratehetmodel = none.");
							if(actuallyUsedImpliedMatrixIndex > 0){
								//the specs are being added as we read and create subsets, so we can add them for the implied matrices
								//as we go
								claSpecs.push_back(ClaSpecifier(nextMatrixNum + actuallyUsedImpliedMatrixIndex, nextMatrixNum + actuallyUsedImpliedMatrixIndex, nextMatrixNum + actuallyUsedImpliedMatrixIndex));
								//clone the current datasubset info, which applies to all of the implied matrices within this effective matrix
								dataSubInfo.push_back(dataSubInfo[dataChunk]);
								//also clone the modspec.  This isn't really necessary (or good) except that the number of states is stored by the modspecs
								if(conf.linkModels)
									modSpecSet.AddModSpec(conf.configModelSets[0]);
								else//there may be only a single model set specified, but no linkage
									modSpecSet.AddModSpec(conf.configModelSets[conf.configModelSets.size() > 1 ? dataChunk : 0]);
								modSpec = modSpecSet.GetModSpec(modSpecSet.NumSpecs() - 1);
								}
							modSpec->SetNStates(impliedMatrix);
							}
						else if(modSpec->IsOrientedGap()){
							if(modSpec->IsGammaRateHet() || modSpec->IsFlexRateHet())
								throw ErrorException("Sorry, rate heterogeneity cannot be used with gap models yet.\n\tSet ratehetmodel = none.");
							if(modSpecSet.InferSubsetRates())
								outman.UserMessage("WARNING - YOU SHOULD TURN OFF SUBSET SPECIFIC RATE ESTIMATION WHEN USING GAP MODELS");
							if(actuallyUsedImpliedMatrixIndex > 0){
								//the specs are being added as we read and create subsets, so we can add them for the implied matrices
								//as we go
								claSpecs.push_back(ClaSpecifier(nextMatrixNum + actuallyUsedImpliedMatrixIndex, nextMatrixNum + actuallyUsedImpliedMatrixIndex, nextMatrixNum + actuallyUsedImpliedMatrixIndex));
								//clone the current datasubset info, which applies to all of the implied matrices within this effective matrix
								dataSubInfo.push_back(dataSubInfo[dataChunk]);
								//also clone the modspec.  This isn't really necessary (or good) except that the number of states is stored by the modspecs
								if(conf.linkModels)
									modSpecSet.AddModSpec(conf.configModelSets[0]);
								else//there may be only a single model set specified, but no linkage
									modSpecSet.AddModSpec(conf.configModelSets[conf.configModelSets.size() > 1 ? dataChunk : 0]);
								modSpec = modSpecSet.GetModSpec(modSpecSet.NumSpecs() - 1);
								}
							}
						else if(modSpec->IsCodon()){
							rawPart.AddSubset(data);
							const NucleotideData *nuc = dynamic_cast<NucleotideData *>(data);
							CodonData *dat;
							if(nuc != NULL)
								dat = new CodonData(nuc, modSpec->geneticCode);
							else throw ErrorException("Attempted to create codon matrix from non-nucleotide data");

							//this probably shouldn't go here, but...
							if(modSpec->IsF1x4StateFrequencies()) dat->SetF1X4Freqs();
							else if(modSpec->IsF3x4StateFrequencies()) dat->SetF3X4Freqs();
							else if(modSpec->IsEmpiricalStateFrequencies()) dat->SetCodonTableFreqs();

							data = dat;
							}
						else if(modSpec->IsCodonAminoAcid()){
							rawPart.AddSubset(data);
							const NucleotideData *nuc = dynamic_cast<NucleotideData *>(data);
							AminoacidData *dat;
							if(nuc != NULL)
								dat = new AminoacidData(nuc, modSpec->geneticCode);
							else throw ErrorException("Attempted to translate to amino acids from non-nucleotide data");

							data = dat;
							}
				
						dataPart.AddSubset(data);

						if(modSpec->IsMkTypeModel()){
							outman.UserMessage("\tSubset of data with %d states:", impliedMatrix);
							string chars;
							data->GetStringOfOrigDataColumns(chars);
							outman.UserMessage("\t  chars%s", chars.c_str());
							}

						if(conf.combineAdjacentIdenticalGapPatterns && (modSpec->IsOrientedGap() || modSpec->IsBinaryNotAllZeros())){
							if(conf.usePatternManager)
								throw ErrorException("Sorry, the pattern manager can't be used with gap collapsing currently");
							data->EliminateAdjacentIdenticalColumns();
							}

						data->ProcessPatterns();

						dataSubInfo[dataChunk + actuallyUsedImpliedMatrixIndex].totalCharacters = data->TotalNChar();
						dataSubInfo[dataChunk + actuallyUsedImpliedMatrixIndex].uniqueCharacters = data->NChar();
						actuallyUsedImpliedMatrixIndex++;
						}
					}
					//subset specific rates will be set if:
					//1. subsetspecificrates = 1 in the conf
					// and 
					//2a. a partition is actually specified via multiple char blocks or a charpart
					// and/or
					//2b. nstate (Mk) model is specified, characters have different numbers of observed states
					//(2b. is what needs to be taken care of here because we don't know whether there will
					//be implied blocks in advance)
					if(conf.subsetSpecificRates && modSpecSet.InferSubsetRates() == false)
						if(actuallyUsedImpliedMatrixIndex > 1)
							modSpecSet.SetInferSubsetRates(true);
				}
			
			//this depends on the fact that an extra taxon slot was allocated by not yet used
			if(modSpecSet.AnyOrientedGap()){
				NxsTaxaBlock *tax = reader.GetTaxaBlock(0);
				if(!tax->IsAlreadyDefined("ROOT"))
					dataPart.AddDummyRoots();
				}
	
			outman.UserMessage("\n###################################################");
			//allocate the population
			pop = new Population();
			pop->usedNCL = usedNCL;
			pop->Setup(&conf, &dataPart, &rawPart, 1, (validateMode == true ? -1 : 0));
			pop->SetOutputDetails();

			outman.UserMessage("STARTING RUN");
			if(runTests){
				outman.UserMessage("starting internal tests...");
					pop->RunTests();
				outman.UserMessage("******Successfully completed tests.******");
				return 0;
				}

			if(conf.optimizeInputOnly)
				conf.runmode = 11;

			if(validateMode){
				//validate mode skips some allocation in pop::Setup, and then executes pop::ValidateInput,
				//which is essentially a stripped down version of pop::SeedPopWithStartingTree
					pop->ValidateInput(1);
				outman.UserMessage("VALIDATION COMPLETE. Check output above for information and possible errors.");
				}
				//the runmodes are essentially a hidden way of causing different (often very different) program
				//behavior at runtime.  not really for user consumption
			else if(conf.runmode != 0){
				if(conf.runmode == 1)
						pop->ApplyNSwaps(10);
				else if(conf.runmode == 7)
						pop->VariableStartingTreeOptimization(false);
				else if(conf.runmode == 9)
						pop->VariableStartingTreeOptimization(true);
				else if(conf.runmode == 8){
					throw ErrorException("Sorry, site rate estimation is not yet implemented in this version.");
#ifdef OPEN_MP
					throw ErrorException("can't estimate site rates in openmp version!");
#endif
#ifndef ALLOW_SINGLE_SITE
					throw ErrorException("the program must be compiled with ALLOW_SINGLE_SITE defined in defs.h to use site rate estimation (runmode = 8)!");
#endif
					pop->OptimizeSiteRates();
					}
				else if(conf.runmode == 11){
					pop->OptimizeInputAndWriteSitelikelihoods();
					}
				else if(conf.runmode == 12){
					pop->OptimizeInputAndWriteSitelikelihoodsAndTryRootings();
					}
				else if(conf.runmode > 20){
						pop->GenerateTreesOnly(conf.runmode);
					}
					else if(conf.runmode > 1) //this is runmodes 2-6
						pop->SwapToCompletion(conf.startOptPrec);
				}
			else{
					//if no checkpoint files are actually found conf->restart will be set to zero
					if(pop->conf->restart) pop->conf->restart = pop->ReadStateFiles();

					pop->SetOutputDetails();
					if(pop->conf->bootstrapReps == 0){//NOT bootstrapping
						pop->PerformSearch();
						}
					else pop->Bootstrap();
					pop->FinalizeOutputStreams(2);
					}
				dataPart.Delete();
				if(pop != NULL){
					delete pop;
					pop = NULL;
				}
			}catch(ErrorException &err){
				if(outman.IsLogSet() == false){
					outman.SetLogFile("ERROR.log");
					if(interactive == false) UsageMessage(argv[0]);
					}
				outman.UserMessage("\nERROR: %s\n\n", err.message);
				if(pop != NULL){
					pop->FinalizeOutputStreams(0);
					pop->FinalizeOutputStreams(1);
					pop->FinalizeOutputStreams(2);
					}

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
		dataPart.Delete();
		modSpecSet.Delete();
	
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
	//		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			sprintf(filename, "memcheck%d.txt", rank);
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

