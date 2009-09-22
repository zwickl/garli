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
#define MINOR_VERSION "r583"

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

#ifdef MPI_VERSION
#include "mpi.h"
#endif

OutputManager outman;
bool interactive;
 
vector<ClaSpecifier> claSpecs;
vector<DataSubsetInfo> dataSubInfo;

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
#ifndef SUBROUTINE_GARLI	
	outman.UserMessage("Usage: %s [OPTION] [config filename]", execName);
	outman.UserMessage("Options:");
	outman.UserMessage("  -i, --interactive	interactive mode (allow and/or expect user feedback)");
	if(interactive) outman.UserMessage("        (interactive is the default for the version you are running)");
	outman.UserMessage("  -b, --batch		batch mode (do not expect user input)");
	if(!interactive) outman.UserMessage("        (batch is the default for the version you are running)");
	outman.UserMessage("  -v, --version		print version information and exit");
	outman.UserMessage("  -h, --help		print this help and exit");
	outman.UserMessage("  -t			run internal tests (requires dataset and config file)");
	outman.UserMessage("NOTE: If no config filename is passed on the command line the program\n   will look in the current directory for a file named \"garli.conf\"");
#else
	outman.UserMessage("Usage: The syntax for launching MPI jobs varies between systems");
	outman.UserMessage("Most likely it will look something like the following:");
	outman.UserMessage("  mpirun [MPI OPTIONS] %s -[# of times to execute config file]", execName);
	outman.UserMessage("Specifying the number of times to execute the config file is mandatory.");
	outman.UserMessage("This version will expect a config file named \"garli.conf\".");
	outman.UserMessage("Consult your cluster documentation for details on running MPI jobs");
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

#elif defined( SUBROUTINE_GARLI )
int SubGarliMain(int rank)	
	{
	int argc=1;
	char **argv=NULL;
	//clear out whatever is in the reader already - it might be full if a single
	//process has called SubGarliMain multiple times
	GarliReader &reader = GarliReader::GetInstance();
	reader.ResetReader();
#else
int main( int argc, char* argv[] )	{
#endif

	CREATE_MEMCHK{//memory leak detecting trick - no overhead when turned off

	#ifdef MPI_VERSION
		MPIMain(argc, argv);
	#endif

	string conf_name;
#ifndef SUBROUTINE_GARLI
	conf_name = "garli.conf";
#else
//	char temp[100];
//	sprintf(temp, "run%d.conf", rank);
//	conf_name = temp;
	//use the same config here too
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
						outman.UserMessage("%s PARTITION Version %.2f.%s", PROGRAM_NAME, MAJOR_VERSION, MINOR_VERSION);
#ifdef SUBROUTINE_GARLI
						outman.UserMessage("MPI run distributing version");
#endif
#ifdef OPEN_MP
						outman.UserMessage("OpenMP multithreaded version");
#endif
#ifdef SINGLE_PRECISION_FLOATS
						outman.UserMessage("Single precision floating point version");
#endif
						outman.UserMessage("Partitioned/Mkv model testing version");
						outman.UserMessage("Copyright Derrick J. Zwickl 2005-2008");
						outman.UserMessage("zwickl@nescent.org");
						exit(0);
						}
					else if(!_stricmp(argv[curarg], "-h") || !_stricmp(argv[curarg], "--help")){
						UsageMessage(argv[0]);
						exit(0);
						}
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
		//PARTITION
		DataPartition dataPart;
		DataPartition rawPart;
		SequenceData *data = NULL;

		try{
			MasterGamlConfig conf;
			bool confOK;
			outman.UserMessage("Reading config file %s", conf_name.c_str());
			confOK = ((conf.Read(conf_name.c_str()) < 0) == false);
			if(confOK == false) throw ErrorException("Error in config file...aborting");

			string datafile = conf.datafname;

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
			
			outman.UserMessage("Running BOINC GARLI, PARTITION version %.2f.%s (Partitioning Branch) (July 2009)\n", MAJOR_VERSION, MINOR_VERSION);
			if(confOK && conf.restart == true) outman.UserMessage("Found BOINC checkpoint files.  Restarting....\n");

			boinc_resolve_filename(datafile.c_str(), buffer, 2048);
			datafile = buffer;
#else	//not BOINC
			//set up the screen log file
			if(confOK == true){
				//changing this to always append to the .screen.log after a restart
				sprintf(temp_buf, "%s.screen.log", conf.ofprefix.c_str());
				//if(conf.restart == false) sprintf(temp_buf, "%s.screen.log", conf.ofprefix.c_str());
				//else sprintf(temp_buf, "%s.restart%d.screen.log", conf.ofprefix.c_str(), CheckRestartNumber(conf.ofprefix));
				}
			else sprintf(temp_buf, "ERROR.log");

			if(conf.restart) outman.SetLogFileForAppend(temp_buf);
			else outman.SetLogFile(temp_buf);
	#ifdef SUBROUTINE_GARLI
			//MPI search forking version
			outman.UserMessage("Running GARLI, PARTITION version %.2f.%s (Partitioning Branch) (July 2009)\n", MAJOR_VERSION, MINOR_VERSION);
			outman.UserMessage("->MPI Parallel Version<-\nNote: this version divides a number of independent runs across processors.");
			outman.UserMessage("It is not the multipopulation parallel Garli algorithm.\n(but is generally a better use of resources)"); 

	#else	//nonMPI version
			outman.UserMessage("Running GARLI, PARTITION version %.2f.%s (Partitioning Branch) (July 2009)\n", MAJOR_VERSION, MINOR_VERSION);
	#endif

#endif //end of BOINC / nonBOINC

#ifdef OPEN_MP
			outman.UserMessage("OpenMP multithreaded version for multiple processors/cores"); 
#endif
			outman.UserMessage("###################################################");
			outman.UserMessage("THIS IS A TESTING VERSION FOR PARTITIONED MODELS AND THE MKV MORPHOLOGY MODEL.");
			outman.UserMessage("IT APPEARS TO BE WORKING PROPERLY, BUT PLEASE LET ME KNOW OF ANY PROBLEMS AT:");
			outman.UserMessage("		          garli.support@gmail.com");
			outman.UserMessage("        see this hidden page for details on partitioned usage:");
			outman.UserMessage("      https://www.nescent.org/wg_garli/Partition_testing_version");
			outman.UserMessage("       and this page for details on Mkv mophology model usage:");
			outman.UserMessage("      https://www.nescent.org/wg_garli/Mkv_morphology_model");
			outman.UserMessage("             CHECK WITH ME BEFORE PUBLISHING WITH IT");
			outman.UserMessage("                         !!!!THANKS!!!!\n");
			outman.UserMessage("###################################################");

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

#ifdef SINGLE_PRECISION_FLOATS
			outman.UserMessage("Single precision floating point version\n");
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
			pop.usedNCL = reader.ReadData(datafile.c_str(), *modSpecSet.GetModSpec(0));
			
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
				if(conf.linkModels && (modSpec->IsNState() || modSpec->IsNState()))
					throw ErrorException("Model linkage cannot be used with Mk/Mkv models (nor does it\n\tneed to be, since there are no estimated parameters.\n\tSet linkmodels = 0");

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

				dataSubInfo[dataChunk].Report();
				outman.UserMessage("");

				// Create the data object
				//for nstate data the effective matrices will be further broken up into implied matrices that each have the same number of observed states
				//the implied matrix number will be that number of states
				int actuallyUsedImpliedMatrixIndex = 0;
				int maxObservedStates = effectiveMatrices[dataChunk].first->GetMaxObsNumStates(false);
				//for Mk the impliedMatrix number is the number of states
				for(int impliedMatrix = 2;impliedMatrix < ((modSpec->IsNState() || modSpec->IsNStateV()) ? maxObservedStates + 1 : 3);impliedMatrix++){
					if(modSpec->IsNState() || modSpec->IsNStateV())
						data = new NStateData(impliedMatrix, modSpec->IsNStateV());
					else if(modSpec->IsAminoAcid() && modSpec->IsCodonAminoAcid() == false)
						data = new AminoacidData();
					else //all other data will be read into a DNA matrix and
						//then converted if necessary
						data = new NucleotideData();
					//if no charpart was specified, the second argument here will be empty
					data->CreateMatrixFromNCL(effectiveMatrices[dataChunk].first, effectiveMatrices[dataChunk].second);

#ifdef SINGLE_PRECISION_FLOATS
					if(modSpec->NState() || modSpec->NStateV()) throw ErrorException("Sorry, Mk/Mkv type models have not yet been tested with single precision.");
#endif
				
					if(data->NChar() == 0){
						//if there weren't any characters with a certain number of states,
						//just get rid of the matrix.  This could in theory also work for
						//totally excluded subsets, but that gets complicated because it
						//isn't clear how the indexing of models specified in the config
						//file should work
						assert(modSpec->IsNState() || modSpec->IsNStateV());
						outman.UserMessage("NOTE: No characters found with %d observed states.", impliedMatrix);
						delete data;
						}
					else{//now we have a data matrix object created, already filtered for the correct sites or number of states
						if(modSpec->IsNState() || modSpec->IsNStateV()){
#ifdef OPEN_MP
							throw ErrorException("Sorry, discrete Mk type models cannot currently be used with the OpenMP version");
#endif
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
								else
									modSpecSet.AddModSpec(conf.configModelSets[dataChunk]);
								modSpec = modSpecSet.GetModSpec(modSpecSet.NumSpecs() - 1);
								}
							modSpec->SetNStates(impliedMatrix);
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

						if(modSpec->IsNState() || modSpec->IsNStateV())
							outman.UserMessage("Subset of data with %d states:", impliedMatrix);

						//this accounts for the dummy character stuck into each data subset
						//for Mkv.  We don't want the screen output to include it.
						int mkvDiff = 0;
						if(modSpec->IsNStateV()) mkvDiff = 1;

						data->Summarize();
						outman.UserMessage("\tSummary of data or data subset:");
						outman.UserMessage("\t%5d sequences.", data->NTax());
						outman.UserMessage("\t%5d constant characters.", data->NConstant() - mkvDiff);
						outman.UserMessage("\t%5d parsimony-informative characters.", data->NInformative());
						outman.UserMessage("\t%5d autapomorphic characters.", data->NAutapomorphic());
						int total = data->NConstant() + data->NInformative() + data->NAutapomorphic();
						if(data->NMissing() > 0){
							outman.UserMessage("\t%5d characters were completely missing or ambiguous (removed).", data->NMissing());
							outman.UserMessage("\t%5d total characters (%d before removing empty columns).", total, data->GapsIncludedNChar() - mkvDiff);
							}
						else outman.UserMessage("\t%5d total characters.", total - mkvDiff);
						
						outman.flush();
						
						data->Collapse();
						outman.UserMessage("\t%5d unique patterns in compressed data matrix.\n", data->NChar() - mkvDiff);

						dataSubInfo[dataChunk + actuallyUsedImpliedMatrixIndex].totalCharacters = data->TotalNChar();
						dataSubInfo[dataChunk + actuallyUsedImpliedMatrixIndex].uniqueCharacters = data->NChar();
						actuallyUsedImpliedMatrixIndex++;

						//DJZ 1/11/07 do this here now, so bootstrapped weights aren't accidentally stored as orig
						data->ReserveOriginalCounts();
						
						if(!modSpec->IsNState() && !modSpec->IsNStateV())
							data->DetermineConstantSites();
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

			outman.UserMessage("\n###################################################");
			pop.Setup(&conf, &dataPart, &rawPart, 1, 0);
			pop.SetOutputDetails();

			outman.UserMessage("###################################################\nSTARTING RUN");
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
				if(pop.conf->restart) pop.ReadStateFiles();

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

