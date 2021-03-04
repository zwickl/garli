// GARLI version 2.1 source code
// Copyright 2005-2014 Derrick J. Zwickl
// email garli.support@gmail.com
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

#ifdef MPI_VERSION
	#include <mpi.h>
#endif
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <signal.h>
#include <stdint.h>

using namespace std;

#ifdef WIN32
#include <conio.h>
#include <windows.h>
#include <winbase.h>
#endif

#ifdef MAC_FRONTEND
#import <Foundation/Foundation.h>
#import "MFEInterfaceClient.h"
#endif

#include "defs.h"
#include "population.h"
#include "individual.h"
#include "sequencedata.h"
#include "tree.h"
#include "funcs.h"
#include "clamanager.h"
#ifdef USE_BEAGLE
#include "calculationmanager.h"
#endif
#include "stopwatch.h"
#include "bipartition.h"
#include "adaptation.h"
#include "errorexception.h"
#include "outputman.h"
#include "model.h"
#include "garlireader.h"

#ifdef ENABLE_CUSTOM_PROFILER
#include "utility.h"
extern Profiler ProfIntInt;
extern Profiler ProfIntTerm;
extern Profiler ProfTermTerm;
extern Profiler ProfRescale;
extern Profiler ProfScoreInt;
extern Profiler ProfScoreTerm;
extern Profiler ProfIntDeriv;
extern Profiler ProfTermDeriv;
extern Profiler ProfCalcPmat;
extern Profiler ProfCalcEigen;
extern Profiler ProfModDeriv;
extern Profiler ProfNewton;
extern Profiler ProfEQVectors;
#endif

extern OutputManager outman;
extern bool interactive;
bool swapBasedTerm = false;

int memLevel;
int calcCount=0;
int optCalcs;

ModelSpecificationSet modSpecSet;

ofstream outf, paupf;
int tempGlobal=1;
bool uniqueSwapTried;

FLOAT_TYPE globalBest;

#undef PERIODIC_SCORE_DEBUG

#undef NNI_SPECTRUM

#undef MASTER_DOES_SUBTREE

bool output_tree=false;

int CheckRestartNumber(const string str);
int debug_mpi(const char* fmt, ...);
int QuitNow();
void InterruptMessage( int );
void ClearDebugLogs();

int askQuitNow = 0;

int QuitNow()
{
	char ch = '?';
	cerr << endl << "Quit? (y/n) -->";
	do {
		cin.get(ch);
	} while( ch != 'y' && ch != 'n' );

        if( ch == 'n' ) askQuitNow = 0;

	return ( ch == 'y' ? 1 : 0 );
}

#ifdef WIN32
// A function to get a single character from the windows console.  Provided by POL.
// Prompts user with string s, then returns the first character typed. If any problems arise
// (e.g. cannot obtain handle to console input buffer), bails out by returning the null
// character (i.e. '\0').
char AskUser(std::string s)
   {
   HANDLE h;
   DWORD num_chars_read, new_console_mode, prev_console_mode;
   char char_buffer[2]; // may be able to get away with [1]

   // Output the prompt string
   std::cerr << s << std::endl;

   // Get handle to console input buffer
   h = GetStdHandle(STD_INPUT_HANDLE);
   if (h == INVALID_HANDLE_VALUE)
       return '\0';

   // Save the current input mode (will restore it before we leave this function)
   if (!GetConsoleMode(h, &prev_console_mode) )
       return '\0';

   // Set new console mode. There are five mode flags defined in wincon.h (ENABLE_LINE_INPUT, ENABLE_ECHO_INPUT,
   // ENABLE_PROCESSED_INPUT, ENABLE_WINDOW_INPUT and ENABLE_MOUSE_INPUT), only ENABLE_PROCESSED_INPUT is useful
   // to us, and we specifically want to avoid ENABLE_LINE_INPUT because it requires the user to press the enter
   // key before ReadConsole returns (much better to have this function return the instant the user presses any
   // key).
   new_console_mode = ENABLE_PROCESSED_INPUT;
   if (!SetConsoleMode(h, new_console_mode))
       return '\0';

   // Read 1 character and place it in char_buffer. num_chars_read should be 1 afterwards. Note that
   // the last argument is reserved and must be NULL.
   if (!ReadConsole(h, char_buffer, 1, &num_chars_read, NULL))
       return '\0';

   // Be nice and return console mode to its previous value
   if (!SetConsoleMode(h, prev_console_mode))
       return '\0';

   return char_buffer[0];
   }
#endif

void InterruptMessage( int )
{
	askQuitNow = 1;
}

void TurnOnSignalCatching()
{//if SIGINT (generally Ctrl-C) isn't already set to be ignored, set it to the custom handler 
	if( signal( SIGINT, SIG_IGN ) != SIG_IGN ){
		signal( SIGINT, InterruptMessage );
		}
}

void TurnOffSignalCatching()
{//if SIGINT (generally Ctrl-C) isn't already set to be ignored, set it back to the default
	if( signal( SIGINT, SIG_IGN ) != SIG_IGN ){
		signal( SIGINT, SIG_DFL );
		}
}

bool CheckForUserSignal(){
	//this will be set if the user raises a signal with ctrl-C
	if(askQuitNow == 1){
		char c;
		if(interactive == false){
			//The run will begin terminating gracefully after this returns, but turn off further catching
			//in case the user wants to fully kill the run fully
			TurnOffSignalCatching();
			return true;
			}
		else{
#if defined (WIN32)
			c = AskUser("Perform final branch-length optimization and terminate now? (y/n)");
	#else
			outman.UserMessage("Perform final branch-length optimization and terminate now? (y/n)");
	#ifdef MAC_FRONTEND
			NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
			BOOL shouldQuit = [[MFEInterfaceClient sharedClient] programShouldTerminate];
			c = shouldQuit ? 'y' : 'n';
			[pool release];
	#else
			c = getchar();
	#endif
	#endif
			if(c=='y'){
				//as above, give up further catching
				TurnOffSignalCatching();
	#ifdef MAC
				cin.get();
	#endif   
				return true;
				}
			else{
				//the user changed their mind
				askQuitNow = 0;
				TurnOnSignalCatching();
				outman.UserMessage("continuing ...");
	#ifndef MAC_FRONTEND
	#ifndef WIN32
				cin.get();
	#endif
	#endif
				return false;
				}
			}
		}
	return false;
	}

//
//
// Methods for class Population
//
//

void ClearDebugLogs(){
	//most of the debug logs just append to the current log, so clear them out here
	#ifndef NDEBUG

//	ofstream pert("pertreport.log");
//	pert.close();

#ifdef MPI_VERSION
	ofstream subrec("subrec.log");
	subrec.close();
	ofstream optp("partscores.log");
	optp.close();
#endif
/*	ofstream brak("brakdebug.log");
	brak.close();
	ofstream brak2("brakmiss.log");
	brak2.close();
	ofstream opt("optimization.log");
	opt.close();
	ofstream optt("opttrees.tre");
	optt.close();
	ofstream opts("optscores.log");
	opts.close();

	ofstream optb("blendeb.log");
	optb.close();
*/
#endif
	}

Population::~Population()
{
	if(indiv != NULL){
		for (unsigned i = 0; i < total_size; ++i)	{
			for (unsigned j = 0; j < total_size; ++j)	{
				if (newindiv[i].treeStruct == indiv[j].treeStruct)	{
					newindiv[i].treeStruct = NULL;
					break;
					}
				}
			}
		}
	
	if( indiv!=NULL )
		MEM_DELETE_ARRAY(indiv); // indiv has length params.nindivs

	if( newindiv!=NULL )
		MEM_DELETE_ARRAY(newindiv); // newindiv has length params.nindivs

	ClearStoredTrees();

	if( cumfit!=NULL ) {
		for( unsigned i = 0; i < total_size; i++ )
			MEM_DELETE_ARRAY(cumfit[i]); // cumfit[i] has length 2
		MEM_DELETE_ARRAY(cumfit); // cumfit has length params.nindivs
	}

	if( treeString!=NULL)
		MEM_DELETE_ARRAY(treeString);

	for(vector<Tree*>::iterator vit=unusedTrees.begin();vit!=unusedTrees.end();vit++){
		delete *vit;
		}
	unusedTrees.clear();

	if(claMan!=NULL){
		delete claMan;
		}
	
	if(paraMan!=NULL){
		delete paraMan;
		}	
#ifdef INCLUDE_PERTURBATION
	if(pertMan!=NULL){
		delete pertMan;
		}
#endif

	if(Bipartition::str!=NULL) delete []Bipartition::str;
	
	for(vector<Tree*>::iterator delit=unusedTrees.begin();delit!=unusedTrees.end();delit++)
		delete *delit;
		
	if(adap!=NULL) delete adap;

	Tree::attemptedSwaps.ClearAttemptedSwaps();

//This shouldn't have been getting deleted.  It was created as a local in main and then just
//aliased in Population
//	if(rawPart != NULL)
//		rawPart->Delete();
		
}

void Population::ErrorMsg( char* msgstr, int len )
{
	switch( error )
	{
		case nomem:
			strncpy( msgstr, "not enough memory", len );
			break;
		case nofile:
			strncpy( msgstr, "parameter file not found", len );
			break;
		case baddimen:
			strncpy( msgstr, "bad dimensions specified", len );
			break;
		default:
			strncpy( msgstr, "undocumented error", len );
	}
}

void Population::CheckForIncompatibleConfigEntries(){
	//DEBUG - fill this in better

	//PARTITION - disallow a number of things that aren't implemented/tested with partitioned models
	if(dataPart->NumSubsets() > 1){
		if(conf->linkModels && modSpecSet.GetModSpec(0)->IsEmpiricalStateFrequencies())
			throw ErrorException("Sorry, empirical state frequencies can't be used with partitioned models when models are linked");
		}

	for(int ms = 0;ms < modSpecSet.NumSpecs();ms++){
		ModelSpecification *modSpec = modSpecSet.GetModSpec(ms);
		//if no model mutations will be performed, parameters cannot be estimated.
		if(conf->modWeight == ZERO_POINT_ZERO){
			if(modSpec->fixStateFreqs == false) 
				throw(ErrorException("if model mutation weight is set to zero,\nstatefrequencies cannot be set to estimate!"));
			if(modSpec->includeInvariantSites == true && modSpec->fixInvariantSites == false) 
				throw(ErrorException("if model mutation weight is set to zero,\ninvariantsites cannot be set to estimate!"));
			if(modSpec->IsAminoAcid() == false && modSpec->Nst() > 1 && modSpec->fixRelativeRates == false) 
				throw(ErrorException("if model mutation weight is set to zero, ratematrix\nmust be fixed or 1rate!"));
			if((modSpec->numRateCats > 1 && modSpec->IsFlexRateHet() == false && modSpec->fixAlpha == false && modSpec->IsCodon() == false) || (modSpec->IsCodon() && !modSpec->fixOmega)) 
				throw(ErrorException("if model mutation weight is set to zero,\nratehetmodel must be set to gammafixed, nonsynonymousfixed or none!"));
			}
		if((modSpec->IsNStateV() || modSpec->IsOrderedNStateV() || modSpec->IsBinaryNotAllZeros() || modSpec->IsOrientedGap()) && (_stricmp(conf->streefname.c_str(), "stepwise") == 0))
			throw ErrorException("Sorry, stepwise addition starting trees currently cannot be used when\n\ta conditioned model (datatype = standardvariable,\n\tstandardvariableordered, binarynotallzeros or indelmixturemodel)\n\tis used for any data.\n\tTry streefname = random, or provide your own starting tree.");
		if(conf->inferInternalStateProbs && ! (modSpec->IsNucleotide() || modSpec->IsAminoAcid() || modSpec->IsCodon()))
			throw ErrorException("Sorry, internal states can currently only be inferred for nucleotide, amino acid and codon models");
		}

	if(conf->inferInternalStateProbs && conf->bootstrapReps > 0) 
		throw(ErrorException("You cannont infer internal states during a bootstrap run!"));
	if(conf->outputSitelikelihoods > 0 && conf->bootstrapReps > 0) 
		throw(ErrorException("You cannont output site likelihoods during a bootstrap run!"));
	if(conf->startOptPrec < conf->minOptPrec)
		throw ErrorException("startoptprec must be equal to or greater than minoptprec");

	if(!conf->checkpoint && conf->workPhaseDivision)
		throw ErrorException("workphasedivision mode only makes sense if checkpoints are written (writecheckpoints = 1)");
	}

void Population::Setup(GeneralGamlConfig *c, DataPartition *d, DataPartition *rawD, int nprocs, int r){
	bool validateMode = false;
	if(r < 0){
		validateMode = true;
		r = 0;
		}
	stopwatch.Start();

	//most of the allocation occurs here or in children
	//set various things
	rank=r;
	conf=c;
	dataPart = d;
	rawPart = rawD;

	subtreeNode=0;

	CheckForIncompatibleConfigEntries();

	//put info that was read from the config file in its place

	if(rank == 0) total_size = conf->nindivs + nprocs-1;
	else total_size = conf->nindivs;
	swapTermThreshold = conf->swapTermThreshold;
	if(swapTermThreshold != 0)
		//this is a global that Tree needs access to
		swapBasedTerm = true;
	else 
		swapBasedTerm = false;

	//set two model statics
	Model::mutationShape = conf->gammaShapeModel;

	//check and warn if different codes have been selected for different subsets - this is experimental
	for(vector<ClaSpecifier>::iterator c = claSpecs.begin();c != claSpecs.end();c++){
		if(modSpecSet.GetModSpec((*c).modelIndex)->IsCodon()){
			for(vector<ClaSpecifier>::iterator c2 = c+1;c2 != claSpecs.end();c2++){
				if(modSpecSet.GetModSpec((*c2).modelIndex)->IsCodon()){
					if(modSpecSet.GetModSpec((*c).modelIndex)->geneticCode != modSpecSet.GetModSpec((*c2).modelIndex)->geneticCode){
						outman.UserMessage("\n################\nWARNING: Different genetic codes have been specified among partition subsets.");
						outman.UserMessage("This is experimental - check your results carefully!!!\n################\n");
						}
					}
				}
			}
		}

#ifdef INPUT_RECOMBINATION
	total_size = conf->nindivs + NUM_INPUT;
#endif

	adap=new Adaptation(conf);

#ifdef INCLUDE_PERTURBATION
	pertMan = new PerturbManager(conf);
#endif

	//instantiate the ParallelManager
	if(rank==0){
		MasterGamlConfig *mastConf = (MasterGamlConfig*) (conf);
		paraMan = new ParallelManager(dataPart->NTax(), nprocs, mastConf);
		}
	
	//use RTTI to check if the data subsets are nuclotide, and if so make ambig strings
	for(int ds = 0;ds < dataPart->NumSubsets();ds++){
		NucleotideData *nuc = dynamic_cast<NucleotideData *>(dataPart->GetSubset(ds));
		if(nuc != NULL)
			nuc->MakeAmbigStrings();
		}
		
	//allocate the treeString
	//remember that we also encode internal node numbers sometimes
	FLOAT_TYPE taxsize=log10((FLOAT_TYPE) ((FLOAT_TYPE)dataPart->NTax())*dataPart->NTax()*2);
	stringSize=(int)((dataPart->NTax()*2)*(10+DEF_PRECISION)+taxsize);
	treeString=new char[stringSize];
	treeString[stringSize - 1]='\0';

	//allocate the indiv array
	indiv = new Individual[total_size];
	newindiv = new Individual[total_size];

	for (unsigned i = conf->nindivs; i < total_size; i++)	{
		indiv[i].reproduced = indiv[i].willreproduce = 1;
		newindiv[i].reproduced = newindiv[i].willreproduce = 1;
		indiv[i].parent=i;
		newindiv[i].parent=i;	
		}

	cumfit = new FLOAT_TYPE*[total_size];
	for(unsigned i = 0; !error && i < total_size; i++ )
		cumfit[i] = new FLOAT_TYPE[2];

	//instantiate the clamanager and figure out how much memory to snatch
	FLOAT_TYPE memToUse;
	//this gives a bit of leeway in normal runs, when total mem usage may get significantly higher than the actual CLA usage
	//but not much is used when just scoring/optimizing one tree
	FLOAT_TYPE memUsageMult = ((conf->scoreOnly || conf->optimizeInputOnly) ? 1.05 : 1.25);
	outman.UserMessage("NOTE: Unlike many programs, the amount of system memory that Garli will\nuse can be controlled by the user.");
	if(conf->availableMemory > 0){
		outman.UserMessage("(This comes from the availablememory setting in the configuration file.");
		outman.UserMessage("Availablememory should NOT be set to more than the actual amount of ");
		outman.UserMessage("physical memory that your computer has installed)");
		memToUse=(FLOAT_TYPE)(1.0/memUsageMult)*conf->availableMemory;
		}
	else{
		outman.UserMessage("\nMemory to be used for conditional likelihood arrays specified as %.1f MB", conf->megsClaMemory);
		memToUse=conf->megsClaMemory;
		}
		
	const int KB = 1024;

	int numNodesPerIndiv = dataPart->NTax()-2;
	int idealClas =  3 * total_size * numNodesPerIndiv;
	int numClas;	

	int L0=(int) (numNodesPerIndiv * total_size * 2);//a downward and one upward set for each tree
	int L1=(int) (numNodesPerIndiv * total_size + 2*total_size + numNodesPerIndiv); //at least a downward set and a full root set for every tree, plus one other set
	int L2=(int) (numNodesPerIndiv * 2.0 + 2*total_size);//a downward set for the best, one other full set and enough for each root direction
	int L3;
	if(conf->scoreOnly || conf->optimizeInputOnly){
		L3=(int) (numNodesPerIndiv * 1.0 + 2);//one full set plus a few extra
		}
	else{
		L3=(int) (numNodesPerIndiv * 1.5 - 2 + 2*total_size);//one full set, enough to reserve at least all of the full internals of the 
															 //best indiv and enough for each root
		}
	//moved this instantiation to Population Constructor, so it can be called earlier
	//calcMan = new CalculationManager();
	outman.UserMessage("");
#ifdef USE_BEAGLE
	calcMan->OutputBeagleResources();

#ifndef BEAGLEPART
	if (modSpecSet.NumSpecs() > 1)
		throw ErrorException("still working on true partitioned beagle support. This version should work fine with a single subset");
#endif
#endif //USE_BEAGLE
	//For partitioned beagle each subset must have the same number of conditionals allocated, 
	//and the memory levels, etc will be the same.  In this mode the specified amount of memory 
	//is per instance (= per subset), so determine how many clas the largest subset (in terms of required size) 
	//could have within that restriction
	double claSizePerNodeKB = 0.0;
#ifdef USE_BEAGLE
	double subsetClaSizePerNodeKB = 0.0;
	vector<int> subsetClaSizes;
	for (vector<ClaSpecifier>::iterator subsetSpec = claSpecs.begin(); subsetSpec != claSpecs.end(); subsetSpec++) {
		subsetClaSizePerNodeKB = indiv[0].modPart.CalcRequiredSubsetCLAsizeKB(subsetSpec->claIndex, dataPart, true);
		subsetClaSizes.push_back(subsetClaSizePerNodeKB);
		claSizePerNodeKB = max(subsetClaSizePerNodeKB, claSizePerNodeKB);
	}
#else
	claSizePerNodeKB = indiv[0].modPart.CalcRequiredCLAsizeKB(dataPart);
#endif

	int maxClas = (int)((memToUse * KB) / claSizePerNodeKB);
	if(maxClas >= L0){
		numClas = min(maxClas, idealClas);
		memLevel = 0;		
		}
	else{
		numClas=maxClas;
	 	if(maxClas >= L1) memLevel = 1;
	 	else if(maxClas >= L2) memLevel = 2;
	 	else if(maxClas >= L3) memLevel = 3;
	 	else memLevel=-1;
		}

	outman.precision(4);
	outman.UserMessage("\nFor this dataset:");
	outman.UserMessage(" Mem level		availablememory setting");
	outman.UserMessage("  great			    >= %.0f MB", ceil(L0 * (claSizePerNodeKB/(FLOAT_TYPE)KB)) * memUsageMult);
	outman.UserMessage("  good			approx %.0f MB to %.0f MB", ceil(L0 * ((FLOAT_TYPE)claSizePerNodeKB/KB)) * memUsageMult - 1, ceil(L1 * ((FLOAT_TYPE)claSizePerNodeKB/KB)) * memUsageMult);
	outman.UserMessage("  low			approx %.0f MB to %.0f MB", ceil(L1 * ((FLOAT_TYPE)claSizePerNodeKB/KB)) * memUsageMult - 1, ceil(L2 * ((FLOAT_TYPE)claSizePerNodeKB/KB)) * memUsageMult);
	outman.UserMessage("  very low		approx %.0f MB to %.0f MB", ceil(L2 * ((FLOAT_TYPE)claSizePerNodeKB/KB)) * memUsageMult - 1, ceil(L3 * ((FLOAT_TYPE)claSizePerNodeKB/KB)) * memUsageMult);
	outman.UserMessage("the minimum required availablememory is %.0f MB", ceil(L3 * ((FLOAT_TYPE)claSizePerNodeKB/KB)) * memUsageMult );

	if(conf->scoreOnly || conf->optimizeInputOnly){
		outman.UserMessage("\nNOTE: Less memory is required when scoring or optimizing fixed trees.\n\tminimum of %.0f availablememory would be required to search\n", ceil(((int) (numNodesPerIndiv * 1.5 - 2 + 2*total_size)) * ((FLOAT_TYPE)claSizePerNodeKB/KB)) * 1.25);
		}

	outman.UserMessage("\nYou specified that Garli should use at most %.1f MB of memory.", conf->availableMemory);
	outman.UserMessage("\nwhen using the BEAGLE library, this is the maximum amount of memory per model partition (instance).");

	if( ! (conf->scoreOnly || conf->optimizeInputOnly)){
		if(memLevel == 0)
			outman.UserMessage("**Your memory level is: great (you don't need to change anything)**");
		else if(memLevel == 1)
			outman.UserMessage("**Your memory level is: good (you don't need to change anything)**");
		else if(memLevel == 2)
			outman.UserMessage("**Your memory level is: low\n\t(you may want to increase the availablememory setting)**");
		else if(memLevel == 3)
			outman.UserMessage("**Your memory level is: very low\n\t(if possible, you should increase the availablememory setting)**");
		else if(memLevel == -1)
			outman.UserMessage("**NOT ENOUGH MEMORY\n\t(you must increase the availablememory setting)**");
		}
	outman.UserMessage("\n#######################################################");
	if (memLevel == -1 && !validateMode)
		throw ErrorException("Not enough memory specified in config file (availablememory)!");

#ifndef USE_BEAGLE
	//increasing this more to allow for the possiblility of needing a set for all nodes for both the indiv and newindiv arrays
//if we do tons of recombination 
	idealClas *= 2;
	if (!validateMode)
		claMan = new ClaManager(dataPart->NTax() - 2, numClas, idealClas, &indiv[0].modPart, dataPart);
#else //BEAGLE
	//loop over subsets (for partitioned beagle)
	double totalMemAllocated = 0.0;
	for (vector<ClaSpecifier>::iterator subsetSpec = claSpecs.begin(); subsetSpec != claSpecs.end(); subsetSpec++) {
		SequenceData *subsetData = dataPart->GetSubset(subsetSpec->dataIndex);
		ModelSpecification *subsetModSpec = modSpecSet.GetModSpec(subsetSpec->modelIndex);
		if (subsetModSpec->IsMkTypeModel() || subsetModSpec->IsOrientedGap())
			throw ErrorException("Invalid datatype!\nSorry, only sequence datatypes (DNA, amino acid and codon) are supported by the BEAGLE-enabled version of GARLI");

		//outman.UserMessage("\nfor this partition subset, Garli will actually use approx. %.1f MB of memory", memUsageMult*(FLOAT_TYPE)numClas*(FLOAT_TYPE)claSizePerNodeKB/(FLOAT_TYPE)KB);
		double subsetMemUsage = memUsageMult * (FLOAT_TYPE)numClas * (FLOAT_TYPE)subsetClaSizes[subsetSpec->claIndex] / (FLOAT_TYPE)KB;
		outman.UserMessage("\nFor this partition subset, Garli will actually use approx. %.1f MB of memory", subsetMemUsage);
		totalMemAllocated += subsetMemUsage;
/*
	outman.precision(4);
	outman.UserMessage("allocating memory...\nusing %.1f MB for conditional likelihood arrays.  Memlevel= %d", (FLOAT_TYPE)numClas*(FLOAT_TYPE)claSizePerNode/(FLOAT_TYPE)MB, memLevel);
	outman.UserMessage("For this dataset:");
	outman.UserMessage("level 0: >= %.0f megs", ceil(L0 * (claSizePerNode/(FLOAT_TYPE)MB)));
	outman.UserMessage("level 1: %.0f megs to %.0f megs", ceil(L0 * ((FLOAT_TYPE)claSizePerNode/MB))-1, ceil(L1 * ((FLOAT_TYPE)claSizePerNode/MB)));
	outman.UserMessage("level 2: %.0f megs to %.0f megs", ceil(L1 * ((FLOAT_TYPE)claSizePerNode/MB))-1, ceil(L2 * ((FLOAT_TYPE)claSizePerNode/MB)));
	outman.UserMessage("level 3: %.0f megs to %.0f megs", ceil(L2 * ((FLOAT_TYPE)claSizePerNode/MB))-1, ceil(L3 * ((FLOAT_TYPE)claSizePerNode/MB)));
	outman.UserMessage("not enough mem: <= %.0f megs\n", ceil(L3 * ((FLOAT_TYPE)claSizePerNode/MB))-1);
*/

	//increasing this more to allow for the possiblility of needing a set for all nodes for both the indiv and newindiv arrays
	//if we do tons of recombination 
	idealClas *= 2;
		if (!validateMode && !claMan)
		claMan=new ClaManager(dataPart->NTax()-2, numClas, idealClas, &indiv[0].modPart, dataPart);

		CalculationManager::SetClaManager(claMan);
		//both the tips and internal branches need pmats, hence the x2
		int idealPmats = total_size * (2 * (dataPart->NTax() - 1));
		idealPmats *= 2;
		//for now using as many pmats as pmat holders
		//pmatMan = new PmatManager(numClas * 2, idealPmats, (modSpec.numRateCats + (modSpec.includeInvariantSites ? 1 : 0)), modSpec.nstates);
		if(!pmatMan)
			pmatMan = new PmatManager(idealPmats, idealPmats, (subsetModSpec->numRateCats + (subsetModSpec->includeInvariantSites ? 1 : 0)), subsetModSpec->nstates);
		CalculationManager::SetPmatManager(pmatMan);

#ifdef BEAGLEPART 
		CalculationManager::SetData(dataPart);
#else
		CalculationManager::SetData(subsetData);
#endif

		//invariable class needs to be treated as extra rate for beagle
		//calcMan->SetBeagleDetails(conf->gpuBeagle, conf->singlePrecBeagle, conf->rescaleBeagle, conf->ofprefix);
		calcMan->SetBeagleDetails(conf->preferredBeagleFlags, conf->requiredBeagleFlags, conf->deviceNumBeagle, conf->ofprefix);
#ifndef BEAGLEPART
		//calcMan->InitializeBeagle(data->NTax(), numClas, idealClas, subsetData->NStates(), subsetData->NChar(), (subsetModSpec.numRateCats + (subsetModSpec.includeInvariantSites ? 1 : 0)));
		calcMan->InitializeBeagleInstance(dataPart->NTax(), numClas, idealClas, subsetData->NStates(), subsetData->NChar(), (subsetModSpec->numRateCats + (subsetModSpec->includeInvariantSites ? 1 : 0)));
#else
		calcMan->AddSubsetInstance(numClas, idealClas, subsetData, subsetModSpec, subsetSpec->modelIndex);
#endif
	} // end loop over modelParts

	outman.UserMessage("\n\nThe total memory allocated for likelihood calculations is approx %.1f MB", totalMemAllocated);



	CalculationManager::SetClaManager(claMan);
	//CalculationManager::SetData(dataPart);
	NodeClaManager::SetClaManager(claMan);
	NodeClaManager::SetPmatManager(pmatMan);
#endif //!USE_BEAGLE

	//setup the bipartition statics
	Bipartition::SetBipartitionStatics(dataPart->NTax());

	//set the tree statics
#ifdef USE_BEAGLE
	Tree::SetTreeStatics(claMan, pmatMan, calcMan, dataPart, conf);
#else
	Tree::SetTreeStatics(claMan, dataPart, conf);
#endif

	//load any constraints
	GetConstraints();

	//try to get nexus starting tree/trees from file, which we don't want to do within the PerformSearch loop
	if((_stricmp(conf->streefname.c_str(), "random") != 0)  && (_stricmp(conf->streefname.c_str(), "stepwise") != 0))
		if(FileIsNexus(conf->streefname.c_str())){
			LoadNexusStartingConditions();
			}
	}

void Population::LoadNexusStartingConditions(){
	GarliReader & reader = GarliReader::GetInstance();
	NxsTaxaBlock *tax = NULL;
	NxsTreesBlock *treesblock = NULL;

	if(reader.GetNumTaxaBlocks() == 1)
		tax = reader.GetTaxaBlock(0);
	else //I think this check happens in NCL as well, but best to be safe
		throw ErrorException("multiple non-identical taxa blocks have been read");
	//this actually is also checked in SetTreeStatics
	if(Tree::rootWithDummy && !tax->IsAlreadyDefined("ROOT")){
		string n = "ROOT";
		tax->AppendNewLabel(n);
		}
	if(usedNCL && strcmp(conf->streefname.c_str(), conf->datafname.c_str()) == 0){
		//in this case we should have already read in the tree when getting the data, so check that we have either one
		//trees block for this taxa block or a garli block
		if(reader.GetNumTreesBlocks(tax) == 0 && reader.FoundModelString() == false)
			throw ErrorException("No nexus trees block or Garli block was found in file %s,\n     which was specified as source of starting tree and/or model", conf->streefname.c_str());
		else if(reader.GetNumTreesBlocks(tax) > 1)
			throw ErrorException("Expecting only one trees block in file %s (not sure which to use)", conf->streefname.c_str());
		else if(reader.GetNumTreesBlocks(tax) == 1)
			startingTreeInNCL = true;
		else startingTreeInNCL = false;
		}
	else{
		//use NCL to get trees from the specified file
		outman.UserMessage("Loading starting model and/or tree from file %s", conf->streefname.c_str());
		//it isn't easy to remove a previous trees block in factory mode, so we need to do this
		int initNumTreesBlocks = reader.GetNumTreesBlocks(tax);
		try{
		reader.ReadFilepath(conf->streefname.c_str(), MultiFormatReader::NEXUS_FORMAT);
			}
		catch (const NxsException & x){
			throw ErrorException("%s", x.msg.c_str());
			}
		int afterNumTreesBlocks = reader.GetNumTreesBlocks(tax);;
		if(afterNumTreesBlocks - initNumTreesBlocks > 1){//we added more than one trees block
			throw ErrorException("Expecting only one trees block in file %s (not sure which to use)", conf->streefname.c_str());
			}
		//otherwise we want the last one because others may have been read with the data
		else if(afterNumTreesBlocks == initNumTreesBlocks)//we didnt' add any tree blocks
			startingTreeInNCL = false;
		else //we found exactly one trees block.  WE NEED TO BE SURE THAT WE USE THE LATEST ONE LATER in SeeedPop
			startingTreeInNCL = true;

		//we read the file, but didn't find either
		if(startingTreeInNCL == false && reader.FoundModelString() == false)
			throw ErrorException("No nexus trees block or Garli block was found in file %s,\n     which was specified as the source of starting model and/or tree", conf->streefname.c_str());
		}
	if(reader.FoundModelString()) 
		startingModelInNCL = true;
	}

//return population more or less to what it looked like just after Setup()
void Population::Reset(){
	if(adap != NULL) delete adap;
	adap=new Adaptation(conf);
	lastTopoImprove = lastPrecisionReduction = gen = 0;
	//conf->restart indicates whether the current rep was restarted, so if we complete one rep and then
	//move on to another it should be false
	conf->restart = false;
	finishedRep = false;
	finishedGenerations = false;
	genTermination = false;
	bestFitness = prevBestFitness = -(FLT_MAX);
	initialRefinePass = finalRefinePass = 0;

	for(unsigned i=0;i<total_size;i++){
		if(indiv[i].treeStruct != NULL){
			indiv[i].treeStruct->RemoveTreeFromAllClas();
			for(unsigned j=0;j<total_size;j++)//because indivs and newindivs can share 
				//tree structures in some situations, this check is necessary to avoid double deletion
				if(newindiv[j].treeStruct == indiv[i].treeStruct) newindiv[j].treeStruct=NULL;
			delete indiv[i].treeStruct;
			indiv[i].treeStruct=NULL;
			}
		}
	for(unsigned i=0;i<total_size;i++){
		if(newindiv[i].treeStruct != NULL){
			newindiv[i].treeStruct->RemoveTreeFromAllClas();
			delete newindiv[i].treeStruct;
			newindiv[i].treeStruct=NULL;
			}
		}
	Tree::attemptedSwaps.ClearAttemptedSwaps();
	}

void Population::ApplyNSwaps(int numSwaps){

	Individual *ind0 = &newindiv[0];

	ind0->GetStartingConditionsFromFile(conf->streefname.c_str(), 0, dataPart->NTax());
	ind0->treeStruct->modPart = &ind0->modPart;
	//ind0->GetStartingConditionsFromNCL(	File(conf->streefname.c_str(), 0, data->NTax());
	
	Individual *repResult = new Individual(ind0);
	storedTrees.push_back(repResult);
	for(int s=0;s<numSwaps;s++){
		ind0->treeStruct->TopologyMutator(0.01, 10, 0);
		ind0->SetDirty();
		ind0->CalcFitness(0);
		Individual *repResult = new Individual(ind0);
		storedTrees.push_back(repResult);
		}

	WriteStoredTrees("swapped.tre");		
	}

void Population::SwapToCompletion(FLOAT_TYPE optPrecision){
	SeedPopulationWithStartingTree(currentSearchRep);
	InitializeOutputStreams();

	if(conf->runmode == 2)
		indiv[0].treeStruct->DeterministicSwapperByDist(&indiv[0], optPrecision, conf->limSPRrange, false);
	else if(conf->runmode == 3)
		indiv[0].treeStruct->DeterministicSwapperByCut(&indiv[0], optPrecision, conf->limSPRrange, false);
	else if(conf->runmode == 4)
		indiv[0].treeStruct->DeterministicSwapperRandom(&indiv[0], optPrecision, conf->limSPRrange);
	else if(conf->runmode == 5)
		indiv[0].treeStruct->DeterministicSwapperByDist(&indiv[0], optPrecision, conf->limSPRrange, true);
	else if(conf->runmode == 6)
		indiv[0].treeStruct->DeterministicSwapperByCut(&indiv[0], optPrecision, conf->limSPRrange, true);
	else if(conf->runmode == 13)
		indiv[0].treeStruct->GenerateTopologiesAtSprDistance(&indiv[0], optPrecision, conf->limSPRrange);

	bestIndiv = 0;
	FinalOptimization();
	WriteTreeFile(besttreefile.c_str(), -1);
/*	double imp = 999.9;
	do{
		imp = indiv[0].treeStruct->OptimizeAllBranches(optPrecision);
		optPrecision /= 1.5;
		}while(imp > 0.0);
*/	outman.UserMessage("final score: %f, %d sec", indiv[0].treeStruct->lnL, stopwatch.SplitTime());
	}

//this is mainly for debugging purposes, to ensure that we are able to make all trees or all trees 
//compatible with any constraints
void Population::GenerateTreesOnly(int nTrees){
	SeedPopulationWithStartingTree(1);
	InitializeOutputStreams();
	if((_stricmp(conf->streefname.c_str(), "random") == 0)){
		outman.UserMessageNoCR("Making random trees compatible with constraints... ");
		for(int i=0;i<nTrees;i++){
			indiv[0].MakeRandomTree(dataPart->NTax());
			AppendTreeToTreeLog(-1, 0);
			indiv[0].treeStruct->RemoveTreeFromAllClas();
			delete indiv[0].treeStruct;
			indiv[0].treeStruct=NULL;
			if(!(i % 100)) outman.UserMessageNoCR("%d ", i);
			}
		}
	else if((_stricmp(conf->streefname.c_str(), "stepwise") == 0)){
		outman.UserMessageNoCR("Making stepwise trees compatible with constraints... ");
		for(int i=0;i<nTrees;i++){
			indiv[0].MakeStepwiseTree(dataPart->NTax(), conf->attachmentsPerTaxon, adap->branchOptPrecision);
			AppendTreeToTreeLog(-1, 0);
			indiv[0].treeStruct->RemoveTreeFromAllClas();
			delete indiv[0].treeStruct;
			indiv[0].treeStruct=NULL;
			if(!(i % 100)) outman.UserMessageNoCR("%d ", i);
			}
		}
	FinalizeOutputStreams(0);
	}

void Population::RunTests(){
	//test a number of functions to ensure that any code changes haven't broken anything
	//it assumes that Setup has been called
	SeedPopulationWithStartingTree(1);
//	InitializeOutputStreams();

#ifdef NDEBUG
	outman.UserMessage("WARNING: You are running internal tests with NDEBUG defined!\nIt should not be defined for full error checking.");
#endif

	if(conf->bootstrapReps > 0){
		outman.UserMessage("\nBootstrap reweighting...");
		//if this is the first rep
		if(nextBootstrapSeed == 0)
			nextBootstrapSeed = rnd.seed();
		lastBootstrapSeed = nextBootstrapSeed;
		nextBootstrapSeed = dataPart->BootstrapReweight(lastBootstrapSeed, conf->resampleProportion);
		}


	//DEBUG
//	Individual *ind0 = &newindiv[0];
//	Individual *ind1 = &newindiv[1];
	Individual *ind0 = &indiv[0];
	Individual *ind1 = &indiv[1];
	Tree *tree0 = ind0->treeStruct;
	Tree *tree1 = ind1->treeStruct;

	//ind0->MakeRandomTree(data->NTax());
	//ind0->MakeStepwiseTree(dataPart->NTax(), conf->attachmentsPerTaxon, adap->branchOptPrecision);
	//ind0->treeStruct->modPart=&ind0->modPart;

#ifdef SINGLE_PRECISION_FLOATS
	int sigFigs = ceil(log10(-tree0->lnL));
	double tol = pow(10.0f, sigFigs-7) * 2.0;
#else
	double tol = 0.001;
#endif

	//check that the score was correct coming out of MakeStepwiseTree
	FLOAT_TYPE scr = tree0->lnL;
	tree0->MakeAllNodesDirty();
	ind0->SetDirty();
	ind0->CalcFitness(0);

	//this only really tests for major scoring problems in the optimization functions
	scr = tree0->lnL;
	tree0->OptimizeAllBranches(adap->branchOptPrecision);
	assert(tree0->lnL + tol > scr);
	assert(tree0->lnL * 2 < scr);

	//test rescaling
	scr = tree0->lnL;
	int r = Tree::rescaleEvery;
	Tree::rescaleEvery = 2;
	tree0->MakeAllNodesDirty();
	ind0->SetDirty();
	ind0->CalcFitness(0);

	if(FloatingPointEquals(ind0->Fitness(), scr, tol) == false){
		throw ErrorException("Failed rescaling test: freq %d=%f, freq 2=%f", r, scr, ind0->Fitness());
		}
	
	Tree::rescaleEvery = r;

	tree1=new Tree();
	ind1->CopySecByRearrangingNodesOfFirst(tree1, ind0);
	tree1->modPart=&ind1->modPart;

	ind0->SetDirty();
	ind0->CalcFitness(0);

	ind1->SetDirty();
	ind1->CalcFitness(0);

	assert(ind0->Fitness() == ind1->Fitness());
	tree0->MakeAllNodesDirty();
	tree1->MakeAllNodesDirty();

	for(int i=0;i<100;i++){

//BMERGE rerooting isn't working right with beagle currently
#ifndef USE_BEAGLE
		tree0->RerootHere(tree0->GetRandomInternalNode());
		tree1->RerootHere(tree1->GetRandomInternalNode());

		tree0->CalcBipartitions(true);
		tree1->CalcBipartitions(true);

		//check rerooting and bipartition comparisons
		assert(tree0->IdenticalTopologyAllowingRerooting(tree1->root));
#endif
		ind0->SetDirty();
		ind1->SetDirty();

		//check minimal recalculation scoring (proper readjustment of CLAs during rerooting)
		tree0->Score(tree0->GetRandomInternalNode());
		tree1->Score(tree1->GetRandomInternalNode());

		if(FloatingPointEquals(tree0->lnL, tree1->lnL, tol) == false){
			throw ErrorException("failed min recalc test: %f diff vs %f allowed",  tree0->lnL - tree1->lnL, tol);
			}

		//check full rescoring from arbitrary nodes in the trees		
		tree0->MakeAllNodesDirty();
		tree1->MakeAllNodesDirty();
		tree0->Score(tree0->GetRandomInternalNode());
		tree1->Score(tree1->GetRandomInternalNode());

		if(FloatingPointEquals(tree0->lnL, tree1->lnL, tol) == false){
			throw ErrorException("failed score at arbitrary node test: %f diff vs %f allowed",  tree0->lnL - tree1->lnL, tol);
			}

		//check that the derivative funcs are outputing the correct score
		TreeNode *nd = tree0->allNodes[tree0->GetRandomNonRootNode()];
		tree0->CalcDerivativesRateHet(nd->anc, nd);
		nd = tree1->allNodes[tree1->GetRandomNonRootNode()];
		tree1->CalcDerivativesRateHet(nd->anc, nd);

		if(FloatingPointEquals(tree0->lnL, tree1->lnL, tol) == false){
			throw ErrorException("failed derivative scoring test: %f diff vs %f allowed",  tree0->lnL - tree1->lnL, tol);
			}
		}
	}

void Population::ResetMemLevel(int numNodesPerIndiv, int numClas){
	assert(0);
	//Deprecated
/*
	const int KB = 1024;
	const int MB = KB*KB;
	
	int claSizePerNode = (4 * modSpec.numRateCats * data->NChar() * sizeof(FLOAT_TYPE)) + (data->NChar() * sizeof(int));
	int sizeOfIndiv = claSizePerNode * numNodesPerIndiv;
	int idealClas =  3 * total_size * numNodesPerIndiv;

	int L0=(int) (numNodesPerIndiv * total_size * 2);//a downward and one upward set for each tree
	int L1=(int) (numNodesPerIndiv * total_size + 2*total_size + numNodesPerIndiv); //at least a downward set and a full root set for every tree, plus one other set
	int L2=(int) (numNodesPerIndiv * 2.0 + 2*total_size);//a downward set for the best, one other full set and enough for each root direction
	int L3=(int) (numNodesPerIndiv * 1.5 - 2 + 2*total_size);//one full set, enough to reserve at least all of the full internals of the 
													 //best indiv and enough for each root	
	
	if(numClas >= L0) memLevel = 0;
	else if(numClas >= L1) memLevel = 1;
	else if(numClas >= L2) memLevel = 2;
	else if(numClas >= L3) memLevel = 3;
	else memLevel=-1;
	assert(memLevel >= 0);
*/
	}


void Population::GetConstraints(){
	//first see if there are any constraints
	if((strlen(conf->constraintfile.c_str()) != 0) && (_stricmp(conf->constraintfile.c_str(), "none") != 0)){
		if(FileIsNexus(conf->constraintfile.c_str())) throw ErrorException("Sorry, Garli doesn't allow constraint trees in Nexus format.\n     See the manual for proper constraint format.");
		ifstream con(conf->constraintfile.c_str());
		if(con.good() == false) throw ErrorException("Could not open constraint file %s!", conf->constraintfile.c_str());
		if(con.good()){
			outman.UserMessage("Loading constraints from file %s", conf->constraintfile.c_str());
			Tree::LoadConstraints(con, dataPart->NTax());
			}
		}
	}


//This is a stripped down version of SeedPopWithStartingTree that loads and validates
//starting conditions but doesn't score or require CLAs to have been allocated
void Population::ValidateInput(int rep){

	//create the first indiv, and then copy the tree and clas

	//this is really annoying and hacky - the maxPinv value is held by each model, and is data dependent (maxPinv can't be > obs pinv)
	//But, since a single model may apply to multiple data, need to be sure that the maxPinv is > the highest obs pinv of any of them
	//now always setting the model default for each data subset (which due to linkage might reset the model several times), but this 
	//shouldn't be problematic.  Note that the other data dependent model thing is empirical base freqs, but that will be disallowed
	//elsewhere when there is linkage.
	FLOAT_TYPE maxPinv = ZERO_POINT_ZERO;
	for(vector<ClaSpecifier>::iterator c = claSpecs.begin();c != claSpecs.end();c++){
		for(int m = 0;m < indiv[0].modPart.NumModels();m++){
			if((*c).modelIndex == m){
				indiv[0].modPart.GetModel(m)->SetDefaultModelParameters(dataPart->GetSubset((*c).dataIndex));
				if(indiv[0].modPart.GetModel(m)->MaxPinv() > maxPinv) maxPinv = indiv[0].modPart.GetModel(m)->MaxPinv();
				}
			}
		}
	//we should only need to do this crap if the models are linked, but not currently allowing linking of some models but not others
	if(conf->linkModels && modSpecSet.GetModSpec(0)->includeInvariantSites == true){
		assert(indiv[0].modPart.NumModels() == 1);
		if(maxPinv > ZERO_POINT_ZERO == false) throw ErrorException("invariantsites = estimate was specified, but no data subsets contained constant characters!");
		indiv[0].modPart.GetModel(0)->SetMaxPinv(maxPinv);
		indiv[0].modPart.GetModel(0)->SetPinv(maxPinv * 0.25, false);
		}

	//DEBUG - need to stick this in somewhere more natural so that it gets reset after a rep completes
	indiv[0].modPart.Reset();

	//This is getting very complicated.  Here are the allowable combinations.
	//streefname not specified (random or stepwise)
		//Case 1 - no gblock in datafile	
		//Case 2 - found gblock in datafile
	//streefname specified
		//specified file is same as datafile
			//Case 3 - Found trees block only
			//Case 4 - Found gblock only (create random tree)
			//Case 5 - Found both
		//specified file not same as datafile
			//NOTE that all of these are also possible with a gblock found in the datafile
			//3/25/08 Change - a second gblock is not allowed (it will throw an exception
			//upon reading the second in GarliReader::EnteringBlock), nor are both a garli block
			//with the data and model params in the old format in the streefname
			//specified streefname is Nexus
				//Case 6 - Found trees block only
				//Case 7 - Found gblock only (create random tree) (if a gblock was already read it will crap out)
				//Case 8 - Found both (if a gblock was already read it will crap out)
			//specified streefname is not Nexus
				//Case 9 - found a tree
				//Case 10 - found a model (create random tree) (if a gblock was already read it will crap out)
				//Case 11 - found both (if a gblock was already read it will crap out)

	GarliReader & reader = GarliReader::GetInstance();

#ifdef INPUT_RECOMBINATION
	if(0)
#else
	if((_stricmp(conf->streefname.c_str(), "random") != 0) && (_stricmp(conf->streefname.c_str(), "stepwise") != 0))
		//some starting file has been specified - Cases 3-11
#endif
	{
		//we already checked in Setup whether NCL has trees for us.  A starting model in Garli block will
		//be handled below, although both a garli block (in the data) and an old style model specification
		//are not allowed
		if(startingTreeInNCL){//cases 3, 5, 6 and 8
			//CAREFUL here - we may have more than one trees block because a tree could appear with the
			//dataset and in a different starting tree file.  The factory api allows this fine, so we
			//need to be sure to grab the last trees block.  Checking for whether the starting tree
			//file contained multiple trees blocks was already done in LoadNexusStartingConditions
			const NxsTreesBlock *treesblock = reader.GetTreesBlock(reader.GetTaxaBlock(0), reader.GetNumTreesBlocks(reader.GetTaxaBlock(0)) - 1);
			assert(treesblock != NULL);
			//this should verify some aspects of the tree description and change everything to taxon numbers
			treesblock->ProcessAllTrees();
			int numTrees = treesblock->GetNumTrees();
			if(numTrees > 0){
				int treeNum = (rank+rep-1) % numTrees;
				indiv[0].GetStartingTreeFromNCL(treesblock, treeNum, dataPart->NTax());
				outman.UserMessage("Obtained starting tree %d from Nexus", treeNum+1);
				}
			else throw ErrorException("Problem getting tree(s) from NCL!");
			}
		else if(strcmp(conf->streefname.c_str(), conf->datafname.c_str()) != 0 && !FileIsNexus(conf->streefname.c_str())){
			//cases 9-11 if the streef file is not the same as the datafile, and it isn't Nexus
			//use the old garli starting model/tree format
			outman.UserMessage("Obtaining starting conditions from file %s", conf->streefname.c_str());
			indiv[0].GetStartingConditionsFromFile(conf->streefname.c_str(), rank + rep - 1, dataPart->NTax());
			}
		indiv[0].SetDirty();
		}

	if(reader.FoundModelString()) 
		startingModelInNCL = true;

	if(startingModelInNCL || conf->parameterValueString.length() > 0){
		//crap out if we already got some parameters above in an old style starting conditions file
#ifndef SUBROUTINE_GARLI
		if(modSpecSet.GotAnyParametersFromFile() && (currentSearchRep == 1 && (conf->bootstrapReps == 0 || currentBootstrapRep == 1)))
			throw ErrorException("Found model parameters specified in a Nexus GARLI block with the dataset,\n\tand in the starting condition file (streefname).\n\tPlease use one or the other.");
#endif
		if(startingModelInNCL && conf->parameterValueString.length() > 0)
			throw ErrorException("Found model parameters specified in the configuration file and in the dataset or starting condition file (streefname).\n\tPlease use one or the other.");
		//model string from garli block, which could have come either in starting condition file
		//or in file with Nexus dataset.  Cases 2, 4, 5, 7 and 8 come through here.

		string modString;
		if(startingModelInNCL)
			modString = reader.GetModelString();
		else
			modString = conf->parameterValueString;

		if(modString.length() > 0)
			indiv[0].modPart.ReadGarliFormattedModelStrings(modString);

		if(startingModelInNCL)
			outman.UserMessage("Obtained starting or fixed model parameter values from Nexus:");
		else
			outman.UserMessage("Obtained starting or fixed model parameter values from configuration file:");
		}

	//The model params should be set to their initial values by now, so report them
	if(conf->bootstrapReps == 0 || (currentBootstrapRep == 1 && currentSearchRep == 1)){
		outman.UserMessage("MODEL REPORT - Parameters are at their INITIAL values (not yet optimized)");
		indiv[0].modPart.OutputHumanReadableModelReportWithParams();
		}

	outman.UserMessage("Starting with seed=%d\n", rnd.seed());

	//Here we'll error out if something was fixed but didn't appear
	for(int ms = 0;ms < modSpecSet.NumSpecs();ms++){
		const ModelSpecification *modSpec = modSpecSet.GetModSpec(ms);
		if((_stricmp(conf->streefname.c_str(), "random") == 0) || (_stricmp(conf->streefname.c_str(), "stepwise") == 0)){
			//if no streefname file was specified, the param values should be in a garli block with the dataset
			if(modSpec->IsNucleotide() && modSpec->IsUserSpecifiedStateFrequencies() && !modSpec->gotStateFreqsFromFile) 
				throw(ErrorException("state frequencies specified as fixed, but no\n\tGarli block found in %s!!" , conf->datafname.c_str()));
			else if(modSpec->fixAlpha && !modSpec->gotAlphaFromFile) 
				throw(ErrorException("alpha parameter specified as fixed, but no\n\tGarli block found in %s!!" , conf->datafname.c_str()));
			else if(modSpec->fixInvariantSites && !modSpec->gotPinvFromFile) 
				throw(ErrorException("proportion of invariant sites specified as fixed, but no\n\tGarli block found in %s!!" , conf->datafname.c_str()));
			else if(modSpec->IsUserSpecifiedRateMatrix() && !modSpec->gotRmatFromFile) 
				throw(ErrorException("relative rate matrix specified as fixed, but no\n\tGarli block found in %s!!" , conf->datafname.c_str()));
			else if(modSpec->IsCodon() && modSpec->fixOmega && !modSpec->gotOmegasFromFile) 
				throw(ErrorException("rate het model set to nonsynonymousfixed, but no\n\tGarli block found in %s!!" , conf->datafname.c_str()));
			}
		else{
			if((modSpec->IsNucleotide() || modSpec->IsAminoAcid()) && modSpec->IsUserSpecifiedStateFrequencies() && !modSpec->gotStateFreqsFromFile) 
				throw ErrorException("state frequencies specified as fixed, but no\n\tparameter values found in %s or %s!", conf->streefname.c_str(), conf->datafname.c_str());
			else if(modSpec->fixAlpha && !modSpec->gotAlphaFromFile) 
				throw ErrorException("alpha parameter specified as fixed, but no\n\tparameter values found in %s or %s!", conf->streefname.c_str(), conf->datafname.c_str());
			else if(modSpec->fixInvariantSites && !modSpec->gotPinvFromFile) 
				throw ErrorException("proportion of invariant sites specified as fixed, but no\n\tparameter values found in %s or %s!", conf->streefname.c_str(), conf->datafname.c_str());
			else if(modSpec->IsUserSpecifiedRateMatrix() && !modSpec->gotRmatFromFile) 
				throw ErrorException("relative rate matrix specified as fixed, but no\n\tparameter values found in %s or %s!", conf->streefname.c_str(), conf->datafname.c_str());
			else if(modSpec->IsCodon() && modSpec->fixOmega && !modSpec->gotOmegasFromFile) 
				throw ErrorException("rate het model set to nonsynonymousfixed, but no\n\tparameter values found in %s or %s!", conf->streefname.c_str(), conf->datafname.c_str());
			}
		}

	//the treestruct could be null if there was a start file that contained no tree
	if((_stricmp(conf->streefname.c_str(), "random") != 0) && (_stricmp(conf->streefname.c_str(), "stepwise") != 0) && (indiv[0].treeStruct != NULL)){
		bool foundPolytomies = indiv[0].treeStruct->ArbitrarilyBifurcate();
		if(foundPolytomies) outman.UserMessage("WARNING: Polytomies found in start tree.  These were arbitrarily resolved.");
	
		indiv[0].treeStruct->root->CheckTreeFormation();
		indiv[0].treeStruct->root->CheckforPolytomies();
		}
	
	//if there are not mutable params in the model, remove any weight assigned to the model
	if(indiv[0].modPart.NumMutableParams() == 0) {
		if((conf->bootstrapReps == 0 && currentSearchRep == 1) || (currentBootstrapRep == 1 && currentSearchRep == 1))
			outman.UserMessage("NOTE: Model contains no mutable parameters!\nSetting model mutation weight to zero.\n");
		adap->modelMutateProb=ZERO_POINT_ZERO;
		adap->UpdateProbs();
		}
	}

void Population::SeedPopulationWithStartingTree(int rep){
	for(unsigned i=0;i<total_size;i++){
		if(indiv[i].treeStruct != NULL) indiv[i].treeStruct->RemoveTreeFromAllClas();
		if(newindiv[i].treeStruct != NULL) newindiv[i].treeStruct->RemoveTreeFromAllClas();
		}

	//create the first indiv, and then copy the tree and clas

	//this is really annoying and hacky - the maxPinv value is held by each model, and is data dependent (maxPinv can't be > obs pinv)
	//But, since a single model may apply to multiple data, need to be sure that the maxPinv is > the highest obs pinv of any of them
	//now always setting the model default for each data subset (which due to linkage might reset the model several times), but this 
	//shouldn't be problematic.  Note that the other data dependent model thing is empirical base freqs, but that will be disallowed
	//elsewhere when there is linkage.
	FLOAT_TYPE maxPinv = ZERO_POINT_ZERO;
	for(vector<ClaSpecifier>::iterator c = claSpecs.begin();c != claSpecs.end();c++){
		for(int m = 0;m < indiv[0].modPart.NumModels();m++){
			if((*c).modelIndex == m){
				indiv[0].modPart.GetModel(m)->SetDefaultModelParameters(dataPart->GetSubset((*c).dataIndex));
				if(indiv[0].modPart.GetModel(m)->MaxPinv() > maxPinv) maxPinv = indiv[0].modPart.GetModel(m)->MaxPinv();
				}
			}
		}
	//we should only need to do this crap if the models are linked, but not currently allowing linking of some models but not others
	if(conf->linkModels && modSpecSet.GetModSpec(0)->includeInvariantSites == true){
		assert(indiv[0].modPart.NumModels() == 1);
		if(maxPinv > ZERO_POINT_ZERO == false) throw ErrorException("invariantsites = estimate was specified, but no data subsets contained constant characters!");
		indiv[0].modPart.GetModel(0)->SetMaxPinv(maxPinv);
		indiv[0].modPart.GetModel(0)->SetPinv(maxPinv * 0.25, false);
		}

	//DEBUG - need to stick this in somewhere more natural so that it gets reset after a rep completes
	indiv[0].modPart.Reset();

	//This is getting very complicated.  Here are the allowable combinations.
	//streefname not specified (random or stepwise)
		//Case 1 - no gblock in datafile	
		//Case 2 - found gblock in datafile
	//streefname specified
		//specified file is same as datafile
			//Case 3 - Found trees block only
			//Case 4 - Found gblock only (create random tree)
			//Case 5 - Found both
		//specified file not same as datafile
			//NOTE that all of these are also possible with a gblock found in the datafile
			//3/25/08 Change - a second gblock is not allowed (it will throw an exception
			//upon reading the second in GarliReader::EnteringBlock), nor are both a garli block
			//with the data and model params in the old format in the streefname
			//specified streefname is Nexus
				//Case 6 - Found trees block only
				//Case 7 - Found gblock only (create random tree) (if a gblock was already read it will crap out)
				//Case 8 - Found both (if a gblock was already read it will crap out)
			//specified streefname is not Nexus
				//Case 9 - found a tree
				//Case 10 - found a model (create random tree) (if a gblock was already read it will crap out)
				//Case 11 - found both (if a gblock was already read it will crap out)

	GarliReader & reader = GarliReader::GetInstance();

#ifdef INPUT_RECOMBINATION
	if(0)
#else
	if((_stricmp(conf->streefname.c_str(), "random") != 0) && (_stricmp(conf->streefname.c_str(), "stepwise") != 0))
		//some starting file has been specified - Cases 3-11
#endif
	{
		//we already checked in Setup whether NCL has trees for us.  A starting model in Garli block will
		//be handled below, although both a garli block (in the data) and an old style model specification
		//are not allowed
		if(startingTreeInNCL){//cases 3, 5, 6 and 8
			//CAREFUL here - we may have more than one trees block because a tree could appear with the
			//dataset and in a different starting tree file.  The factory api allows this fine, so we
			//need to be sure to grab the last trees block.  Checking for whether the starting tree
			//file contained multiple trees blocks was already done in LoadNexusStartingConditions
			const NxsTreesBlock *treesblock = reader.GetTreesBlock(reader.GetTaxaBlock(0), reader.GetNumTreesBlocks(reader.GetTaxaBlock(0)) - 1);
			assert(treesblock != NULL);
			//this should verify some aspects of the tree description and change everything to taxon numbers
			treesblock->ProcessAllTrees();
			int numTrees = treesblock->GetNumTrees();
			if(numTrees > 0){
				int treeNum = (rank+rep-1) % numTrees;
				indiv[0].GetStartingTreeFromNCL(treesblock, treeNum, dataPart->NTax());
				outman.UserMessage("Obtained starting tree %d from Nexus", treeNum+1);
				}
			else throw ErrorException("Problem getting tree(s) from NCL!");
			}
		else if(strcmp(conf->streefname.c_str(), conf->datafname.c_str()) != 0 && !FileIsNexus(conf->streefname.c_str())){
			//cases 9-11 if the streef file is not the same as the datafile, and it isn't Nexus
			//use the old garli starting model/tree format
			outman.UserMessage("Obtaining starting conditions from file %s", conf->streefname.c_str());
			indiv[0].GetStartingConditionsFromFile(conf->streefname.c_str(), rank + rep - 1, dataPart->NTax());
			}
		indiv[0].SetDirty();
		}

	if(reader.FoundModelString()) 
		startingModelInNCL = true;

	if(startingModelInNCL || conf->parameterValueString.length() > 0){
		//crap out if we already got some parameters above in an old style starting conditions file
#ifndef SUBROUTINE_GARLI
		if(modSpecSet.GotAnyParametersFromFile() && (currentSearchRep == 1 && (conf->bootstrapReps == 0 || currentBootstrapRep == 1)))
			throw ErrorException("Found model parameters specified in a Nexus GARLI block with the dataset,\n\tand in the starting condition file (streefname).\n\tPlease use one or the other.");
#endif
		if(startingModelInNCL && conf->parameterValueString.length() > 0)
			throw ErrorException("Found model parameters specified in the configuration file and in the dataset or starting condition file (streefname).\n\tPlease use one or the other.");
		//model string from garli block, which could have come either in starting condition file
		//or in file with Nexus dataset.  Cases 2, 4, 5, 7 and 8 come through here.

		string modString;
		if(startingModelInNCL)
			modString = reader.GetModelString();
		else
			modString = conf->parameterValueString;

		if(modString.length() > 0)
			indiv[0].modPart.ReadGarliFormattedModelStrings(modString);

		if(startingModelInNCL)
			outman.UserMessage("Obtained starting or fixed model parameter values from Nexus:");
		else
			outman.UserMessage("Obtained starting or fixed model parameter values from configuration file:");
		}

	//Here we'll error out if something was fixed but didn't appear
	for(int ms = 0;ms < modSpecSet.NumSpecs();ms++){
		const ModelSpecification *modSpec = modSpecSet.GetModSpec(ms);
		if((_stricmp(conf->streefname.c_str(), "random") == 0) || (_stricmp(conf->streefname.c_str(), "stepwise") == 0)){
			//if no streefname file was specified, the param values should be in a garli block with the dataset
			if(modSpec->IsNucleotide() && modSpec->IsUserSpecifiedStateFrequencies() && !modSpec->gotStateFreqsFromFile) 
				throw(ErrorException("state frequencies specified as fixed, but no\n\tGarli block found in %s!!" , conf->datafname.c_str()));
			else if(modSpec->fixAlpha && !modSpec->gotAlphaFromFile) 
				throw(ErrorException("alpha parameter specified as fixed, but no\n\tGarli block found in %s!!" , conf->datafname.c_str()));
			else if(modSpec->fixInvariantSites && !modSpec->gotPinvFromFile) 
				throw(ErrorException("proportion of invariant sites specified as fixed, but no\n\tGarli block found in %s!!" , conf->datafname.c_str()));
			else if(modSpec->IsUserSpecifiedRateMatrix() && !modSpec->gotRmatFromFile) 
				throw(ErrorException("relative rate matrix specified as fixed, but no\n\tGarli block found in %s!!" , conf->datafname.c_str()));
			else if(modSpec->IsCodon() && modSpec->fixOmega && !modSpec->gotOmegasFromFile) 
				throw(ErrorException("rate het model set to nonsynonymousfixed, but no\n\tGarli block found in %s!!" , conf->datafname.c_str()));
			}
		else{
			if((modSpec->IsNucleotide() || modSpec->IsAminoAcid()) && modSpec->IsUserSpecifiedStateFrequencies() && !modSpec->gotStateFreqsFromFile) 
				throw ErrorException("state frequencies specified as fixed, but no\n\tparameter values found in %s or %s!", conf->streefname.c_str(), conf->datafname.c_str());
			else if(modSpec->fixAlpha && !modSpec->gotAlphaFromFile) 
				throw ErrorException("alpha parameter specified as fixed, but no\n\tparameter values found in %s or %s!", conf->streefname.c_str(), conf->datafname.c_str());
			else if(modSpec->fixInvariantSites && !modSpec->gotPinvFromFile) 
				throw ErrorException("proportion of invariant sites specified as fixed, but no\n\tparameter values found in %s or %s!", conf->streefname.c_str(), conf->datafname.c_str());
			else if(modSpec->IsUserSpecifiedRateMatrix() && !modSpec->gotRmatFromFile) 
				throw ErrorException("relative rate matrix specified as fixed, but no\n\tparameter values found in %s or %s!", conf->streefname.c_str(), conf->datafname.c_str());
			else if(modSpec->IsCodon() && modSpec->fixOmega && !modSpec->gotOmegasFromFile) 
				throw ErrorException("rate het model set to nonsynonymousfixed, but no\n\tparameter values found in %s or %s!", conf->streefname.c_str(), conf->datafname.c_str());
			}
		if(conf->modWeight == ZERO_POINT_ZERO)
			if(modSpec->IsCodon() && modSpec->gotOmegasFromFile == false) 
				throw(ErrorException("sorry, to turn off model mutations you must provide omega values in a codon model.\nSet modweight to > 0.0 or provide omega values."));
		}
		
	//The model params should be set to their initial values by now, so report them
	if(conf->bootstrapReps == 0 || (currentBootstrapRep == 1 && currentSearchRep == 1)){
		outman.UserMessage("MODEL REPORT - Parameters are at their INITIAL values (not yet optimized)");
		indiv[0].modPart.OutputHumanReadableModelReportWithParams();
		}

	outman.UserMessage("Starting with seed=%d\n", rnd.seed());

	//A random tree specified, or a starting file was specified but contained no tree
	if(_stricmp(conf->streefname.c_str(), "stepwise") == 0){
		if(Tree::constraints.empty()) outman.UserMessage("creating likelihood stepwise addition starting tree...");
		else outman.UserMessage("creating likelihood stepwise addition starting tree (compatible with constraints)...");
		//5/20/08 If we're making a stepwise tree, we depend on the extern globalBest being zero to keep the optimization
		//during the stepwise creation to be localized to just the three branches (the radius optimization only happens if
		//the lnL of created tree is within a threshold of the global best).  Having global best = zero effectively turns
		//off all radius opt.  There was a bug here because it wasn't getting reset before starting search reps after the 
		//first.  This caused the stepwise to be slow, and to not be reproducible when the seed from a rep > 1 was specified
		//as the initial seed for a new run
		globalBest = ZERO_POINT_ZERO;
		//DEBUG - haven't worked this out with the rooted tree yet, since the fake root needs to be in the tree before it is scored
		assert(!indiv[0].treeStruct->rootWithDummy);
		indiv[0].MakeStepwiseTree(dataPart->NTax(), conf->attachmentsPerTaxon, adap->branchOptPrecision);
		}
	else if(_stricmp(conf->streefname.c_str(), "random") == 0 || indiv[0].treeStruct == NULL){
		if(Tree::constraints.empty()) outman.UserMessage("creating random starting tree...");
		else outman.UserMessage("creating random starting tree (compatible with constraints)...");
		indiv[0].MakeRandomTree(dataPart->NTax());
		indiv[0].SetDirty();
		}

	assert(indiv[0].treeStruct != NULL);
	bool foundPolytomies = indiv[0].treeStruct->ArbitrarilyBifurcate();
	if(foundPolytomies) outman.UserMessage("WARNING: Polytomies found in start tree.  These were arbitrarily resolved.");
	
	indiv[0].treeStruct->root->CheckTreeFormation();
	indiv[0].treeStruct->root->CheckforPolytomies();

//BMERGE - rebalancing totally hoses things with beagle.  avoiding for now
#ifndef USE_BEAGLE
	if(!indiv[0].treeStruct->rootWithDummy)
		indiv[0].treeStruct->CheckBalance();
#endif
	indiv[0].treeStruct->modPart=&indiv[0].modPart;
	
	try{
		indiv[0].CalcFitness(0);
		}
	catch(UnscoreableException &ex){
		throw ErrorException("Initial individual unscorable, perhaps due to poorness of starting tree.\n\tTry providing a tree if you previously tried random.");
		}

	//check the current likelihood now to know how accurate we can expect them to be later
#ifdef SINGLE_PRECISION_FLOATS
	Tree::expectedPrecision = pow(10.0, - (double) ((int) FLT_DIG - ceil(log10(-indiv[0].Fitness()))));
#else
	Tree::expectedPrecision = pow(10.0, - (double) ((int) DBL_DIG - ceil(log10(-indiv[0].Fitness()))));
#endif
//	outman.UserMessage("expected likelihood precision = %.4e", Tree::expectedPrecision);

	//if there are not mutable params in the model, remove any weight assigned to the model
	if(indiv[0].modPart.NumMutableParams() == 0) {
		if((conf->bootstrapReps == 0 && currentSearchRep == 1) || (currentBootstrapRep == 1 && currentSearchRep == 1))
			outman.UserMessage("NOTE: Model contains no mutable parameters!\nSetting model mutation weight to zero.\n");
		adap->modelMutateProb=ZERO_POINT_ZERO;
		adap->UpdateProbs();
		}

	outman.precision(10);
	outman.UserMessage("Initial ln Likelihood: %.4f", indiv[0].Fitness());
#ifdef MAC_FRONTEND
	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	[[MFEInterfaceClient sharedClient] didBeginInitializingSearch];
	[pool release];
#endif		

	if(conf->refineStart==true && !conf->scoreOnly){
		//12/26/07 now only passing the first argument here ("optModel") as false if no model muts are used
		//if single parameters are fixed that will be checked in the Refine function itself
		//5/15/14 Moved the initial refinement phase to the Population level, which makes more sense and
		//mirrors the final optimization behavior
		//indiv[0].RefineStartingConditions(adap->modWeight != ZERO_POINT_ZERO, adap->branchOptPrecision);
		InitialOptimization(&indiv[0], adap->modWeight != ZERO_POINT_ZERO, adap->branchOptPrecision);
		indiv[0].CalcFitness(0);
		outman.UserMessage("lnL after optimization: %.4f", indiv[0].Fitness());
		}	

	globalBest=bestFitness=prevBestFitness=indiv[0].Fitness();

	//don't bother allocating any further indivs if we will only use one
	if(conf->optimizeInputOnly || conf->scoreOnly)
		return;

#ifndef INPUT_RECOMBINATION
	for(unsigned i=1;i<total_size;i++){
		if(indiv[i].treeStruct == NULL) 
			indiv[i].treeStruct = new Tree();
		indiv[i].CopySecByRearrangingNodesOfFirst(indiv[i].treeStruct, &indiv[0]);
		indiv[i].treeStruct->modPart=&indiv[i].modPart;
		}
#else
	for(unsigned i=1;i<conf->nindivs;i++){
		if(indiv[i].treeStruct == NULL) 
			indiv[i].treeStruct = new Tree();
		indiv[i].CopySecByRearrangingNodesOfFirst(indiv[i].treeStruct, &indiv[0]);
		indiv[i].treeStruct->modPart=&indiv[i].modPart;
		}
	
	//string inputs="sphinx.input6000.goodmod.tre";
	for(unsigned i=conf->nindivs;i<total_size;i++){
		indiv[i].GetStartingConditionsFromFile(conf->streefname.c_str(), i-conf->nindivs, dataPart->NTax());
		indiv[i].treeStruct->modPart=&indiv[i].modPart;
		indiv[i].SetDirty();
		//indiv[i].RefineStartingConditions((adap->modWeight == ZERO_POINT_ZERO || modSpec->fixAlpha == true) == false, adap->branchOptPrecision);
		indiv[i].CalcFitness(0);
		}
#endif

	CalcAverageFitness();
	}

//Copied almost exactly from Individual::RefineStartingConditions.  For various reasons it is easier 
//to have it at the population level.
void Population::InitialOptimization(Individual *ind, bool optModel, FLOAT_TYPE branchPrec){
	bool optOmega, optAlpha, optFlex, optPinv, optFreqs, optRelRates, optSubsetRates;
	optOmega = optAlpha = optFlex = optPinv = optFreqs = optRelRates = optSubsetRates = false;

	bool optInsDel = false;

	ModelPartition &modPart = ind->modPart;
	Tree *treeStruct = ind->treeStruct;

	if(optModel){
		for(int modnum = 0;modnum < modPart.NumModels();modnum++){
			Model *mod = modPart.GetModel(modnum);
			const ModelSpecification *modSpec = mod->GetCorrespondingSpec();
			if(modSpec->numRateCats > 1 && modSpec->IsNonsynonymousRateHet() == false && modSpec->IsFlexRateHet() == false && modSpec->fixAlpha == false) 
				optAlpha = true;
			if(modSpec->IsFlexRateHet()) 
				optFlex = true;
			if(modSpec->includeInvariantSites && modSpec->fixInvariantSites == false) 
				optPinv = true;
			if(modSpec->IsCodon() && !modSpec->fixOmega) 
				optOmega = true;
			if(modSpec->IsOrientedGap()) 
				optInsDel = true;

			if(modSpec->IsCodon() == false && modSpec->fixStateFreqs == false && modSpec->IsEqualStateFrequencies() == false && modSpec->IsEmpiricalStateFrequencies() == false)
				optFreqs = true;
			//this is the case of forced freq optimization with codon models.  For everything to work they must be set as both not fixed but empirical
			if(modSpec->IsCodon() && modSpec->fixStateFreqs == false && modSpec->IsEqualStateFrequencies() == false && modSpec->IsEmpiricalStateFrequencies() == true)
				optFreqs = true;
			if(modSpec->fixRelativeRates == false && (modSpec->Nst() > 1 || modSpec->IsEstimateAAMatrix() || modSpec->IsTwoSerineRateMatrix()))
				optRelRates = true;
			}
		//oops, bug fixed 10/2/12 - subset rates weren't getting opt in linked models
		//modSpecSet.inferSubsetRates is already getting set only if conf.inferSubsetRates
		//is true and there are multiple matrices, but not necessarily multiple models
		//if(modSpecSet.InferSubsetRates() && modSpecSet.NumSpecs() > 1)
		if(modSpecSet.InferSubsetRates())
			optSubsetRates = true;
		}

	outman.UserMessageNoCR("optimizing: starting branch lengths");
	if(optAlpha) outman.UserMessageNoCR(", alpha shape");
	if(optPinv) outman.UserMessageNoCR(", prop. invar");
	if(optRelRates) outman.UserMessageNoCR(", rel rates");
	if(optFreqs) outman.UserMessageNoCR(", eq freqs");
	if(optOmega) outman.UserMessageNoCR(", dN/dS (aka omega) parameters");
	if(optInsDel){
		outman.UserMessageNoCR(", ins rate");
		outman.UserMessageNoCR(", del rate");
		}
	if(optSubsetRates) outman.UserMessageNoCR(", subset rates");
	outman.UserMessage("...");
	FLOAT_TYPE improve=(FLOAT_TYPE)999.9;
	ind->CalcFitness(0);

	assert(initialRefinePass < 1);
	for(initialRefinePass = 1;improve > branchPrec;initialRefinePass++){
		FLOAT_TYPE alphaOptImprove=0.0, pinvOptImprove = 0.0, omegaOptImprove = 0.0, flexOptImprove = 0.0, optImprove=0.0, scaleOptImprove=0.0, subsetRateImprove=0.0, rateOptImprove=0.0;
		FLOAT_TYPE freqOptImprove=0.0, insDelImprove = 0.0;
		
		ind->CalcFitness(0);
		FLOAT_TYPE passStart = ind->Fitness();
		
		optImprove=treeStruct->OptimizeAllBranches(branchPrec);
		ind->CalcFitness(0);

		FLOAT_TYPE trueImprove = ind->Fitness() - passStart;
		assert(trueImprove >= -1.0);
		if(trueImprove < ZERO_POINT_ZERO) trueImprove = ZERO_POINT_ZERO;

		vector<FLOAT_TYPE> blens;
		treeStruct->StoreBranchlengths(blens);
		scaleOptImprove=treeStruct->OptimizeTreeScale(branchPrec);
		ind->CalcFitness(0);
		//if some of the branch lengths were at the minimum or maximum boundaries the scale optimization
		//can actually worsen the score.  If so, return them to their original lengths.
		if(scaleOptImprove < ZERO_POINT_ZERO){
			treeStruct->RestoreBranchlengths(blens);
			ind->CalcFitness(0);
			scaleOptImprove = ZERO_POINT_ZERO;
			}

		ind->CalcFitness(0);
		if(optModel){
			for(int modnum = 0;modnum < modPart.NumModels();modnum++){
				Model *mod = modPart.GetModel(modnum);
				const ModelSpecification *modSpec = mod->GetCorrespondingSpec();
				if(modSpec->IsCodon()){
					if(!modSpec->fixOmega)
						omegaOptImprove += treeStruct->OptimizeOmegaParameters(branchPrec, modnum);
					}
				else if(mod->NRateCats() > 1){
					if(modSpec->IsFlexRateHet()){//Flex rates
						//no longer doing alpha first, it was too hard to know if the flex rates had been partially optimized
						//already during making of a stepwise tree
						//if(i == 1) rateOptImprove = treeStruct->OptimizeAlpha(branchPrec);
						//if(i == 1 && modSpec.gotFlexFromFile==false) rateOptImprove = treeStruct->OptimizeBoundedParameter(branchPrec, mod->Alpha(), 0, 1.0e-8, 999.9, &Model::SetAlpha);
						flexOptImprove += treeStruct->OptimizeFlexRates(branchPrec, modnum);
						}
					else if(modSpec->fixAlpha == false){//normal gamma
						//rateOptImprove = treeStruct->OptimizeAlpha(branchPrec);
						//do NOT let alpha go too low here - on bad or random starting trees the branch lengths get crazy long
						//rateOptImprove = treeStruct->OptimizeBoundedParameter(branchPrec, mod->Alpha(), 0, 1.0e-8, 999.9, &Model::SetAlpha);
						//alphaOptImprove += treeStruct->OptimizeBoundedParameter(branchPrec, mod->Alpha(), 0, 0.05, 999.9, modnum, &Model::SetAlpha);
						alphaOptImprove += treeStruct->OptimizeBoundedParameter(modnum, branchPrec, mod->Alpha(), 0, 0.05, 999.9, &Model::SetAlpha);
						}
					}
				if(modSpec->includeInvariantSites && !modSpec->fixInvariantSites)
					pinvOptImprove += treeStruct->OptimizeBoundedParameter(modnum, branchPrec, mod->PropInvar(), 0, 1.0e-8, mod->maxPropInvar, &Model::SetPinv);
				if(modSpec->IsOrientedGap()){
					insDelImprove += treeStruct->OptimizeInsertDeleteRates(branchPrec, modnum);
					}
				if(modSpec->IsCodon() == false && modSpec->fixStateFreqs == false && modSpec->IsEqualStateFrequencies() == false && modSpec->IsEmpiricalStateFrequencies() == false)
					freqOptImprove += treeStruct->OptimizeEquilibriumFreqs(branchPrec, modnum);
				if(modSpec->fixRelativeRates == false && (modSpec->Nst() > 1 || modSpec->IsEstimateAAMatrix() || modSpec->IsTwoSerineRateMatrix()))
					rateOptImprove += treeStruct->OptimizeRelativeNucRates(branchPrec, modnum);
				}
			if(optSubsetRates){
				subsetRateImprove += treeStruct->OptimizeSubsetRates(branchPrec);
				}
			}
		improve=scaleOptImprove + trueImprove + alphaOptImprove + pinvOptImprove + flexOptImprove + omegaOptImprove + rateOptImprove + freqOptImprove + subsetRateImprove + insDelImprove;
		outman.precision(8);
		outman.UserMessageNoCR("pass%2d:+%9.3f (branch=%7.2f scale=%6.2f", initialRefinePass, improve, trueImprove, scaleOptImprove);
		if(optOmega) outman.UserMessageNoCR(" omega=%6.2f", omegaOptImprove);
		if(optAlpha) outman.UserMessageNoCR(" alpha=%6.2f", alphaOptImprove);

		if(optFreqs) outman.UserMessageNoCR(" freqs=%6.2f", freqOptImprove);
		if(optRelRates) outman.UserMessageNoCR(" rel rates=%6.2f", rateOptImprove);

		if(optFlex) outman.UserMessageNoCR(" flex=%6.2f", flexOptImprove);
		if(optPinv) outman.UserMessageNoCR(" pinv=%6.2f", pinvOptImprove);
		if(optInsDel){
			outman.UserMessageNoCR(" ins/del=%6.2f", insDelImprove);
			}
		if(optSubsetRates) outman.UserMessageNoCR(" subset rates=%6.2f", subsetRateImprove);
		outman.UserMessageNoCR(")");
		
		UpdateFractionDone(1);
		if(conf->reportRunProgress)
			outman.UserMessageNoCR(" %14.2f %14.2f", 0.01 * (int) ceil(rep_fraction_done * 100), 0.01 * (int) ceil(tot_fraction_done * 100));
		outman.UserMessage("");
		}
	initialRefinePass = -1;
	treeStruct->nodeOptVector.clear();
	}

/* This is deprecated in favor of model based function
void Population::OutputModelReport(){
	//Report on the model setup
	outman.UserMessage("MODEL REPORT:");
	if(modSpec->IsCodon()){
		if(modSpec->IsVertMitoCode()) outman.UserMessage("  Number of states = 60 (codon data, vertebrate mitochondrial code)");
		else if(modSpec->IsInvertMitoCode()) outman.UserMessage("  Number of states = 62 (codon data, invertebrate mitochondrial code)");
		else outman.UserMessage("  Number of states = 61 (codon data, standard code)");
		}
	else if(modSpec->IsAminoAcid())
		outman.UserMessage("  Number of states = 20 (amino acid data)");
	else 
		outman.UserMessage("  Number of states = 4 (nucleotide data)");
	
	if(modSpec->IsAminoAcid() == false){
		if(modSpec->IsCodon() && modSpec->numRateCats == 1) outman.UserMessageNoCR("  One estimated dN/dS ratio (aka omega)\n");
		if(modSpec->IsCodon()) outman.UserMessageNoCR("  Nucleotide Relative Rate Matrix Assumed by Codon Model:\n     ");
		else outman.UserMessageNoCR("  Nucleotide Relative Rate Matrix: ");
		if(modSpec->Nst() == 6){
			if(modSpec->IsArbitraryRateMatrix()) outman.UserMessage("User specified matrix type: %s", modSpec->arbitraryRateMatrixString.c_str());
			else outman.UserMessage("6 rates");
			if(modSpec->fixRelativeRates == true) outman.UserMessage("values specified by user (fixed)");
			}
		else if(modSpec->Nst() == 2) outman.UserMessage("2 rates (transition and transversion)");
		else outman.UserMessage("1 rate");
		}
	else{
		outman.UserMessageNoCR("  Amino Acid Rate Matrix: ");
		if(modSpec->IsJonesAAMatrix()) outman.UserMessage("Jones");
		else if(modSpec->IsDayhoffAAMatrix()) outman.UserMessage("Dayhoff");
		else if(modSpec->IsPoissonAAMatrix()) outman.UserMessage("Poisson");
		else if(modSpec->IsWAGAAMatrix()) outman.UserMessage("WAG");
		else if(modSpec->IsMtMamAAMatrix()) outman.UserMessage("MtMam");
		else if(modSpec->IsMtRevAAMatrix()) outman.UserMessage("MtRev");
		}

	outman.UserMessageNoCR("  Equilibrium State Frequencies: ");
	if(modSpec->IsEqualStateFrequencies()){
		if(modSpec->IsCodon()){
			if(modSpec->IsVertMitoCode()) outman.UserMessage("equal (1/60 = 0.01667, fixed)");
			else if(modSpec->IsInvertMitoCode()) outman.UserMessage("equal (1/62 = 0.01613, fixed)");
			else outman.UserMessage("equal (1/61 = 0.01639, fixed)");
			}
		else if(modSpec->IsAminoAcid())
			outman.UserMessage("equal (0.05, fixed)");
		else 
			outman.UserMessage("equal (0.25, fixed)");
		}
	else if(modSpec->IsF3x4StateFrequencies()) outman.UserMessage("empirical values calculated by F3x4 method (fixed)");
	else if(modSpec->IsF1x4StateFrequencies()) outman.UserMessage("empirical values calculated by F1x4 method (fixed)");
	else if(modSpec->IsEmpiricalStateFrequencies()) outman.UserMessage("empirical values (fixed)");
	else if(modSpec->IsJonesAAFreqs()) outman.UserMessage("Jones");
	else if(modSpec->IsWAGAAFreqs()) outman.UserMessage("WAG");
	else if(modSpec->IsMtMamAAFreqs()) outman.UserMessage("MtMam");
	else if(modSpec->IsMtRevAAFreqs()) outman.UserMessage("MtRev");
	else if(modSpec->IsDayhoffAAFreqs()) outman.UserMessage("Dayhoff");
	else if(modSpec->IsUserSpecifiedStateFrequencies()) outman.UserMessage("specified by user (fixed)");
	else outman.UserMessage("estimated");

	outman.UserMessage("  Rate Heterogeneity Model:");
	if(modSpec->numRateCats == 1){
		if(modSpec->includeInvariantSites == false) outman.UserMessage("    no rate heterogeneity");
		else{
			if(modSpec->fixInvariantSites == true) outman.UserMessage("    only an invariant (invariable) site category,\n    proportion specified by user (fixed)");
			else outman.UserMessage("    only an invariant (invariable) site category,\n    proportion estimated");
			}
		}
	else{
		outman.UserMessageNoCR("    %d ", modSpec->numRateCats);
		if(modSpec->IsNonsynonymousRateHet()){
			outman.UserMessage("nonsynonymous rate categories, rate and proportion of each estimated\n     (this is effectively the M3 model of PAML)");
			}
		else if(modSpec->IsFlexRateHet() == false){
			if(modSpec->fixAlpha == true) outman.UserMessage("discrete gamma distributed rate cats,\n    alpha param specified by user (fixed)");
			else outman.UserMessage("discrete gamma distributed rate cats, alpha param estimated");
			if(modSpec->includeInvariantSites == true){
				if(modSpec->fixInvariantSites == true) outman.UserMessage("    with an invariant (invariable) site category,\n    proportion specified by user (fixed)");				
				else outman.UserMessage("    with an invariant (invariable) site category, proportion estimated");
				}
			}
		else{
			outman.UserMessage("FLEX rate categories, rate and proportion of each estimated");
			if(modSpec->includeInvariantSites == true){
				if(modSpec->fixInvariantSites == true) outman.UserMessage("    with an invariant (invariable) site category,\n    proportion specified by user (fixed)");				
				else outman.UserMessage("    with an invariant (invariable) site category, proportion estimated");
				}
			}
		}
	outman.UserMessage("");
	}
*/
	/*
void Population::WriteStateFiles(){
	char str[100];

	//write the adaptation info checkpoint in binary format
	sprintf(str, "%s.adap.check", conf->ofprefix.c_str());
	ofstream out(str, ios::binary | ios::out);
	adap->WriteToCheckpoint(out);
	out.close();

	//write the state of the population, including the seed, generation, elapsed time,
	//lastTopoImprove and specifications of the current individuals
	sprintf(str, "%s.pop.check", conf->ofprefix.c_str());
	ofstream pout(str);
	pout.precision(10);
	WritePopulationCheckpoint(pout);
	pout.close();

	//if we are keeping track of swaps, write a checkpoint for that
	if(conf->uniqueSwapBias != ONE_POINT_ZERO){
		sprintf(str, "%s.swaps.check", conf->ofprefix.c_str());
		ofstream sout(str);
		Tree::attemptedSwaps.WriteSwapCheckpoint(sout);
		sout.close();
		}	
	}
*/

//#define OLD_CHECK

#ifdef OLD_CHECK
void Population::WriteStateFiles(){
	char name[100];

	//write the adaptation info checkpoint in binary format
	sprintf(name, "%s.adap.check", conf->ofprefix.c_str());
#ifdef BOINC
	char physical_name[256];
	boinc_resolve_filename(name, physical_name, sizeof(physical_name));
	MFILE out;
	out.open(physical_name, "wb");
#else
	ofstream out(name, ios::out | ios::binary);
#endif
	adap->WriteToCheckpoint(out);
	out.close();

	//write the state of the population, including the seed, generation, elapsed time,
	//lastTopoImprove and specifications of the current individuals
	sprintf(name, "%s.pop.check", conf->ofprefix.c_str());
#ifdef BOINC
	MFILE pout;
	boinc_resolve_filename(name, physical_name, sizeof(physical_name));
	pout.open(physical_name, "wb");
#else
	ofstream pout(name, ios::out | ios::binary);
#endif
	WritePopulationCheckpoint(pout);
	pout.close();

	//if we are keeping track of swaps, write a checkpoint for that
	if(conf->uniqueSwapBias != ONE_POINT_ZERO){
		sprintf(name, "%s.swaps.check", conf->ofprefix.c_str());
#ifdef BOINC
		MFILE sout;
		boinc_resolve_filename(name, physical_name, sizeof(physical_name));
		sout.open(physical_name, "wb");
#else
		ofstream sout(name, ios::out | ios::binary);
#endif
		Tree::attemptedSwaps.WriteSwapCheckpoint(sout);
		sout.close();
		}
	}

#else
void Population::WriteStateFiles(){
	char aname[128];
	char pname[128];
	char sname[128];

	sprintf(aname, "%s.adap.check", conf->ofprefix.c_str());
	sprintf(pname, "%s.pop.check", conf->ofprefix.c_str());
	sprintf(sname, "%s.swaps.check", conf->ofprefix.c_str());

	//1. write the adaptation info checkpoint in binary format
	//2. write the state of the population, including the seed, generation, elapsed time,
	//	 lastTopoImprove and specifications of the current individuals
	//3. if we are keeping track of swaps, write a checkpoint for that
#ifdef BOINC
	//The BOINC provided MFILE class handily allows writing to it before it is actually
	//open and attached to any file.  It buffers the information, and flushes it upon closing
	//of the file.  This is good, because all of the checkpoint information can be gathered
	//and then written all at once, making it less likely that the checkpoint will be corrupted 
	//by the program being terminated mid-checkpoint 
	MFILE aout;
	MFILE pout;
	MFILE sout;

	adap->WriteToCheckpoint(aout);
	WritePopulationCheckpoint(pout);

	if(conf->uniqueSwapBias != ONE_POINT_ZERO)
		Tree::attemptedSwaps.WriteSwapCheckpoint(sout);
	
	char aphysical_name[512];
	char pphysical_name[512];
	char sphysical_name[512];
	boinc_resolve_filename(aname, aphysical_name, sizeof(aphysical_name));
	boinc_resolve_filename(pname, pphysical_name, sizeof(pphysical_name));
	boinc_resolve_filename(sname, sphysical_name, sizeof(sphysical_name));

	aout.open(aphysical_name, "wb");
	aout.close();

	pout.open(pphysical_name, "wb");
	pout.close();

	if(conf->uniqueSwapBias != ONE_POINT_ZERO){
		sout.open(sphysical_name, "wb");
		sout.close();
		}

	boinc_checkpoint_completed();
#else
	//it would be nice to be able to do something like what is done for BOINC above (buffering output and
	//then writing all at once), maybe with a stringstream, but I couldn't get that to work
	ofstream aout(aname, ios::out | ios::binary);
	adap->WriteToCheckpoint(aout);
	aout.close();
	
	ofstream pout(pname, ios::out | ios::binary);
	WritePopulationCheckpoint(pout);
	pout.close();

	if(conf->uniqueSwapBias != ONE_POINT_ZERO){
		ofstream sout(sname, ios::out | ios::binary);
		Tree::attemptedSwaps.WriteSwapCheckpoint(sout);
		sout.close();
		}
#endif
	}
#endif

//Returns whether or not checkpoints were actually found and read
bool Population::ReadStateFiles(){
	char name[100];

	//read the adaptation binary checkpoint
	sprintf(name, "%s.adap.check", conf->ofprefix.c_str());
	FILE *in;
#ifdef BOINC
	char physical_name[100];
	boinc_resolve_filename(name, physical_name, sizeof(physical_name));
	in = boinc_fopen(physical_name, "rb");

#else
	if(FileExists(name) == false){
	#if defined(SUBROUTINE_GARLI) || defined(OLD_SUBROUTINE_GARLI)
		//for the MPI version we don't care if checkpoint files weren't found
		return false;
	#else
		throw(ErrorException("Could not find checkpoint file %s!\nEither the previous run was not writing checkpoints (checkpoint = 0),\nthe checkpoint files were moved/deleted or the ofprefix setting\nin the config file was changed.", name));
	#endif
		}
	in = fopen(name, "rb");
#endif
	adap->ReadFromCheckpoint(in);
	fclose(in);

	//Read the population checkpoint
	ReadPopulationCheckpoint();

	//need to reset these here, although really only because asserts check that the values never decrease
	rep_fraction_done = tot_fraction_done = 0.0;

#ifdef BOINC
	boinc_fraction_done(tot_fraction_done);
#endif

	//Read the swap checkpoint, if necessary
	if(conf->uniqueSwapBias != ONE_POINT_ZERO){
		sprintf(name, "%s.swaps.check", conf->ofprefix.c_str());
		FILE *sin;
#ifdef BOINC
		boinc_resolve_filename(name, physical_name, sizeof(physical_name));
		sin = boinc_fopen(physical_name, "rb");
#else
		if(FileExists(name) == false) 
			throw(ErrorException("Could not find checkpoint file %s!\nEither the previous run was not writing checkpoints (checkpoint = 0),\nthe file was moved/deleted or the ofprefix setting\nin the config file was changed.", name));
		sin = fopen(name, "rb");
#endif
		Tree::attemptedSwaps.ReadBinarySwapCheckpoint(sin);
		fclose(sin);
		}
	return true;
	}
/*
void Population::WritePopulationCheckpoint(ofstream &out) {
	long currentSeed = rnd.seed();
	out.write((char*) &currentSeed, sizeof(currentSeed));
	long currentTime = stopwatch.SplitTime();
	out.write((char*) &currentTime, sizeof(currentTime));

	//7/13/07 changing this to calculate the actual size of the chunk of scalars
	//(the number of bytes between the start of the object and the first nonscalar
	//data member) rather than counting the number of each type and adding it up 
	//manually.  This should make it work irrespective of things like memory padding
	//for data member alignment, which could vary between platforms and compilers
	intptr_t scalarSize = (intptr_t) &fraction_done - (intptr_t) this  + sizeof(fraction_done);
	out.write((char*) this, (streamsize) scalarSize);
		
	for(unsigned i=0;i<total_size;i++){
		assert(out.good());
		indiv[i].mod->OutputBinaryFormattedModel(out);
		indiv[i].treeStruct->OutputBinaryFormattedTree(out);
		}
	}
*/

void Population::WritePopulationCheckpoint(OUTPUT_CLASS &out) {
	assert(!timeTermination && !userTermination);
	long currentSeed = rnd.seed();
	out.WRITE_TO_FILE(&currentSeed, sizeof(currentSeed), 1);
	int currentTime = stopwatch.SplitTime();
	out.WRITE_TO_FILE(&currentTime, sizeof(currentTime), 1);

	//7/13/07 changing this to calculate the actual size of the chunk of scalars
	//(the number of bytes between the start of the object and the first nonscalar
	//data member) rather than counting the number of each type and adding it up 
	//manually.  This should make it work irrespective of things like memory padding
	//for data member alignment, which could vary between platforms and compilers
	intptr_t scalarSize = (intptr_t) &rep_fraction_done - (intptr_t) this + sizeof(rep_fraction_done);
	out.WRITE_TO_FILE(this, (streamsize) scalarSize, 1);

	//save the current members of the population
	for(unsigned i=0;i<total_size;i++){
		indiv[i].modPart.WriteModelPartitionCheckpoint(out);
		indiv[i].treeStruct->OutputBinaryFormattedTree(out);
		}

	//write any individuals that we may have stored from previous search reps
	for(vector<Individual*>::iterator it = storedTrees.begin(); it != storedTrees.end() ; it++){
		(*it)->modPart.WriteModelPartitionCheckpoint(out);
		(*it)->treeStruct->OutputBinaryFormattedTree(out);
		}
	}


void Population::ReadPopulationCheckpoint(){
	char str[100];
	sprintf(str, "%s.pop.check", conf->ofprefix.c_str());
	if(FileExists(str) == false) throw(ErrorException("Could not find checkpoint file %s!\nEither the previous run was not writing checkpoints (checkpoint = 0),\nthe file was moved/deleted or the ofprefix setting\nin the config file was changed.", str));

#ifdef BOINC
	char physical_name[100];
	boinc_resolve_filename(str, physical_name, sizeof(physical_name));
	FILE *pin = boinc_fopen(physical_name, "rb");

#else
	FILE *pin = fopen(str, "rb");
#endif

	long seed;
	fread((char *) &seed, sizeof(seed), 1, pin);
	if(ferror(pin) || feof(pin)){//this mainly checks for a zero-byte file
		throw ErrorException("Error reading checkpoint file %s.\n\tA problem may have occured writing the file to disk, or the file may have been overwritten or truncated.\n\tUnfortunately you'll need to start the run again from scratch.", str);
		}
	rnd.set_seed(seed);

	int t;
	fread((char *) &t, sizeof(t), 1, pin);
	stopwatch.AddPreviousTime(t);

	//7/13/07 changing this to calculate the actual size of the chunk of scalars
	//(the number of bytes between the start of the object and the first nonscalar
	//data member) rather than counting the number of each type and adding it up 
	//manually.  This should make it work irrespective of things like memory padding
	//for data member alignment, which could vary between platforms and compilers
	intptr_t scalarSize = (intptr_t) &rep_fraction_done - (intptr_t) this + sizeof(rep_fraction_done);
	fread(this, scalarSize, 1, pin);

	//if were restarting a bootstrap run we need to change to the bootstrapped data
	//now, so that scoring below is correct
	if(conf->bootstrapReps > 0){
		int s = dataPart->BootstrapReweight(lastBootstrapSeed, conf->resampleProportion);
		//this should be the case because what was written to the checkpoint for nextBootstrapSeed
		//should have come out of BootstrapReweight when it was originally called with lastBootstrapSeed
		assert(s == nextBootstrapSeed);
		}

	for(unsigned i=0;i<total_size;i++){
		assert(modSpecSet.NumSpecs() == indiv[i].modPart.NumModelSets());
		for(int m = 0;m < modSpecSet.NumSpecs();m++){
			//it would make more sense to have this happen at a lower level, but the data are needed
			indiv[i].modPart.GetModelSet(m)->SetDefaultModelSetParameters(dataPart->GetSubset(m));
			}
		indiv[i].modPart.ReadModelPartitionCheckpoint(pin);

		indiv[i].treeStruct = new Tree();
		indiv[i].treeStruct->ReadBinaryFormattedTree(pin);
		indiv[i].treeStruct->AssignCLAsFromMaster();
		indiv[i].treeStruct->modPart=&indiv[i].modPart;
		indiv[i].SetDirty();
		indiv[i].treeStruct->root->CheckTreeFormation();
		indiv[i].CalcFitness(0);
		}

	//if we are doing multiple reps, there should have been one tree per completed rep written to file
	//remember that currentSearchRep starts at 1
	for(int i=1;i<(finishedRep == false ? currentSearchRep : currentSearchRep+1);i++){
		Individual *ind = new Individual;
		for(int m = 0;m < modSpecSet.NumSpecs();m++){
			//it would make more sense to have this happen at a lower level, but the data are needed
			ind->modPart.GetModelSet(m)->SetDefaultModelSetParameters(dataPart->GetSubset(m));
			}
		ind->modPart.ReadModelPartitionCheckpoint(pin);

		ind->treeStruct = new Tree();
		ind->treeStruct->ReadBinaryFormattedTree(pin);
		ind->treeStruct->AssignCLAsFromMaster();
		ind->treeStruct->modPart=&ind->modPart;
		ind->SetDirty();
		ind->treeStruct->root->CheckTreeFormation();
		ind->CalcFitness(0);
		ind->treeStruct->RemoveTreeFromAllClas();
		storedTrees.push_back(ind);
		}

	//as far as the TopologyList is concerned, each individual will be considered different
	ntopos = total_size;
	if(fabs(bestFitness - indiv[bestIndiv].Fitness()) > 0.01)
		throw ErrorException("Problem reading checkpoint files.  Scores of stored trees don't match calculated scores.");
	CalcAverageFitness();
	globalBest = bestFitness;
	}

//Depending on the generation, output to various files during the GA search
void Population::WriteGenerationOutput(){
	if(conf->outputMostlyUselessFiles) 
		OutputFate();
	if(conf->logevery > 0 && !(gen % conf->logevery)) 
		OutputLog();
	if(conf->saveevery > 0 && !(gen % conf->saveevery)){
		if(best_output & WRITE_CONTINUOUS){
			string outname = besttreefile;
			outname += ".current";
			WriteTreeFile( outname.c_str(), -1);
			}

		outman.UserMessageNoCR("%-8d %-14.4f   %-9.3f  %6d ", gen, BestFitness(), adap->branchOptPrecision, lastTopoImprove);

		if(swapBasedTerm)
			outman.UserMessageNoCR("%14d ", indiv[bestIndiv].treeStruct->attemptedSwaps.GetUnique());

		if(conf->reportRunProgress)
			outman.UserMessageNoCR("%14.2f %14.2f", 0.01 * (int) ceil(rep_fraction_done * 100), 0.01 * (int) ceil(tot_fraction_done * 100));

		outman.UserMessage("");
		
		if(conf->outputMostlyUselessFiles){
			swapLog << gen << "\t";
			indiv[bestIndiv].treeStruct->attemptedSwaps.SwapReport(swapLog);
			}
		}
	}

void Population::Run(){
//	calcCount=0;
	optCalcs=0;

#ifdef VARIABLE_OPTIMIZATION
//	var << "type\tdist\tinitlnL\tnoBail.01\tnoBail.5\t3B.01\t3B.5\tdef.01\tdefdef\n";
#endif
/*	if(conf->restart == false){
		if(conf->bootstrapReps == 0) outman.UserMessage("Running Genetic Algorithm with initial seed=%d", rnd.init_seed());
		}
	else{
		outman.UserMessage("Restarting Genetic Algorithm from checkpoint");
		outman.UserMessage("generation %d, seed %d, best lnL %.3f", gen, rnd.init_seed(), BestFitness());
		}
*/
#ifdef MAC_FRONTEND
	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	[[MFEInterfaceClient sharedClient] didBeginRun];
	[pool release];
#endif	
	
	CalcAverageFitness();

	outman.precision(6);

	outman.UserMessageNoCR("%-8s %-14s %-8s  %-14s ", "gen", "current_lnL", "precision", "last_tree_imp");

	if(swapBasedTerm)
		outman.UserMessageNoCR("%-14s ", "swaps_on_cur");

	if(conf->reportRunProgress)
		outman.UserMessageNoCR("%-14s %-14s", "rep_prop_done", "tot_prop_done");

	outman.UserMessage("");

	outman.UserMessageNoCR("%-8d %-14.4f   %-9.3f  %6d ", gen, BestFitness(), adap->branchOptPrecision, lastTopoImprove);

	if(swapBasedTerm)
		outman.UserMessageNoCR("%14d ", indiv[bestIndiv].treeStruct->attemptedSwaps.GetUnique());

	if(finishedGenerations == false)
		UpdateFractionDone(2);
	if(conf->reportRunProgress)
		outman.UserMessageNoCR("%14.2f %14.2f", 0.01 * (int) ceil(rep_fraction_done * 100), 0.01 * (int) ceil(tot_fraction_done * 100));

	outman.UserMessage("");

	OutputLog();
	if(conf->outputMostlyUselessFiles) 
		OutputFate();	

	gen++;
	for (; gen < conf->stopgen+1; ++gen){
		//this is set true if the generation loop was exited normally but final optimization was not done
		if(finishedGenerations)
			break;

		NextGeneration();
		UpdateFractionDone(2);
		if(swapBasedTerm){
			if(uniqueSwapTried){
				lastUniqueSwap = gen;
				uniqueSwapTried = false;
				}
			}
		keepTrack();
		
        WriteGenerationOutput();
			
#ifndef BOINC
		userTermination = CheckForUserSignal();
		if(userTermination){
			outman.UserMessage("NOTE: ****Run terminated by user interuption ...");
			break;
			}
#endif

#ifdef PERIODIC_SCORE_DEBUG
		if(gen % 500 == 0 ||gen==1)
			OutputFilesForScoreDebugging(&indiv[bestIndiv], tempGlobal++);
#endif

#ifdef NNI_SPECTRUM
		if(gen % 1000 == 0 || gen==1)
			NNISpectrum(bestIndiv);
#endif

		if(!(gen%adap->intervalLength)){
			outman.precision(10);
			bool reduced=false;
			if(gen-(max(lastTopoImprove, lastPrecisionReduction)) >= adap->intervalsToStore*adap->intervalLength){
				//this allows the program to bail if numPrecReductions < 0, which can be handy to get to this point
				//with checkpointing in and then restart from the same checkpoint with various values of numPrecReductions
				if(adap->numPrecReductions < 0) return;
				reduced=adap->ReducePrecision();
				}
			//optimize params if we just reduced prec or if we are at the min prec and we've run for a while since the last reduction
			if(reduced || ((gen - lastPrecisionReduction >= (adap->intervalLength * 50)) && (gen % (adap->intervalLength * 50) == 0) && (FloatingPointEquals(adap->branchOptPrecision, conf->minOptPrec, 1.0e-8)))){
				if(reduced){
					lastPrecisionReduction=gen;
					outman.UserMessage("Optimization precision reduced ");
					}
				//Added in this optimization of rate het params at prec reduction,
				//mainly to help with optimization in partitioned models
				FLOAT_TYPE improve = 0.0;
				Tree *bestTree = indiv[bestIndiv].treeStruct;
				for(int modnum = 0;modnum < indiv[bestIndiv].modPart.NumModels();modnum++){
					Model *mod = indiv[bestIndiv].modPart.GetModel(modnum);
					const ModelSpecification *modSpec = mod->GetCorrespondingSpec();
					if(modSpec->IsCodon())//optimize omega even if there is only 1
						improve += bestTree->OptimizeOmegaParameters(adap->branchOptPrecision, modnum);
					else if(mod->NRateCats() > 1){
						if(modSpec->IsFlexRateHet()){//Flex rates
							improve += bestTree->OptimizeFlexRates(adap->branchOptPrecision, modnum);
							}
						else if(modSpec->fixAlpha == false){//normal gamma
							improve += bestTree->OptimizeBoundedParameter(modnum, adap->branchOptPrecision, mod->Alpha(), 0, min(mod->Alpha(), 0.05), 999.9, &Model::SetAlpha);
							}
						}
					if(modSpec->includeInvariantSites && !modSpec->fixInvariantSites)
						improve += bestTree->OptimizeBoundedParameter(modnum, adap->branchOptPrecision, mod->PropInvar(), 0, min(mod->PropInvar(), 1.0e-8), mod->maxPropInvar, &Model::SetPinv);
					if(modSpec->IsCodon() == false && modSpec->fixStateFreqs == false && modSpec->IsEqualStateFrequencies() == false && modSpec->IsEmpiricalStateFrequencies() == false){
						FLOAT_TYPE paramOpt = bestTree->OptimizeEquilibriumFreqs(adap->branchOptPrecision, modnum);
						if(paramOpt < ZERO_POINT_ZERO && paramOpt > -1e-8)//avoid printing very slightly negative values
							paramOpt = ZERO_POINT_ZERO;
						improve += paramOpt;
						outman.DebugMessage("eq freq opt = %.4f", paramOpt);
						}
					if(modSpec->fixRelativeRates == false && modSpec->Nst() > 1 && modSpec->IsAminoAcid() == false){
						FLOAT_TYPE paramOpt = bestTree->OptimizeRelativeNucRates(adap->branchOptPrecision, modnum);
						if(paramOpt < ZERO_POINT_ZERO && paramOpt > -1e-8)//avoid printing very slightly negative values
							paramOpt = ZERO_POINT_ZERO;
						improve += paramOpt;
						outman.DebugMessage("rel rates opt = %.4f", paramOpt);
						}
					if(modSpec->IsEstimateAAMatrix() || (modSpec->IsTwoSerineRateMatrix() && !modSpec->fixRelativeRates)){
						FLOAT_TYPE paramOpt = bestTree->OptimizeRelativeNucRates(adap->branchOptPrecision, modnum);
						if(paramOpt < ZERO_POINT_ZERO && paramOpt > -1e-8)//avoid printing very slightly negative values
							paramOpt = ZERO_POINT_ZERO;
						improve += paramOpt;
						outman.DebugMessage("rel rates opt = %.4f", paramOpt);
						}
					if(modSpec->IsOrientedGap()){
						FLOAT_TYPE paramOpt = bestTree->OptimizeInsertDeleteRates(adap->branchOptPrecision, modnum);
						if(paramOpt < ZERO_POINT_ZERO && paramOpt > -1e-8)//avoid printing very slightly negative values
							paramOpt = ZERO_POINT_ZERO;
						improve += paramOpt;
						outman.DebugMessage("ins/del opt = %.4f", paramOpt);
						}
					}
				if(modSpecSet.InferSubsetRates()){
					improve += bestTree->OptimizeSubsetRates(adap->branchOptPrecision);
					}
				if(!(FloatingPointEquals(adap->modWeight, 0.0, 1e-8)))
					outman.UserMessage("   Optimizing parameters...    improved %8.3f lnL", improve);
				/////

				FLOAT_TYPE before=bestFitness;
				//under some conditions (very steep lopsided likelihood curve for branch lengths)
				//the blen opt can actually make the score worse
				bestTree->OptimizeAllBranches(adap->branchOptPrecision);
				indiv[bestIndiv].SetDirty();
				CalcAverageFitness();
				FLOAT_TYPE bImp = bestFitness - before - improve;
				if(bImp < ZERO_POINT_ZERO && bImp > -1e-4)//avoid printing very slightly negative values
					bImp = ZERO_POINT_ZERO;
				outman.UserMessage("   Optimizing branchlengths... improved %8.3f lnL", bImp);
				//This is important so that new better topos can be properly identified in the next gen!
				adap->lastgenscore = BestFitness();
				}

			//automatic termination conditions
			if(conf->enforceTermConditions == true){
				if(swapBasedTerm && !FloatingPointEquals(adap->topoMutateProb, ZERO_POINT_ZERO, max(1.0e-8, GARLI_FP_EPS * 2.0))){
					assert(swapTermThreshold != 0);
					if(swapTermThreshold < 0 && (gen - lastUniqueSwap > abs(swapTermThreshold))){
						break;
						}
					else {
						if(swapTermThreshold > 0 && (gen - lastUniqueSwap > swapTermThreshold)
						&& (gen-max(lastTopoImprove, lastPrecisionReduction) > conf->lastTopoImproveThresh || FloatingPointEquals(adap->topoMutateProb, ZERO_POINT_ZERO, max(1.0e-8, GARLI_FP_EPS * 2.0)))
						&& (gen > adap->intervalsToStore * adap->intervalLength)
						&& adap->improveOverStoredIntervals < conf->improveOverStoredIntervalsThresh
						&& (FloatingPointEquals(adap->branchOptPrecision, adap->minOptPrecision, max(1.0e-8, GARLI_FP_EPS * 2.0)) || adap->numPrecReductions==0)){
							if(adap->topoMutateProb > ZERO_POINT_ZERO) 
								outman.UserMessage("Reached termination condition!\nlast topological improvement at gen %d", lastTopoImprove);
							else 
								outman.UserMessage("Reached termination condition!\n");
							outman.UserMessage("Improvement over last %d gen = %.5f", adap->intervalsToStore*adap->intervalLength, adap->improveOverStoredIntervals);
							outman.UserMessage("Last new topology swap at gen %d", lastUniqueSwap);
							break;
							}
						}
					}
				else{
					if((gen-max(lastTopoImprove, lastPrecisionReduction) > conf->lastTopoImproveThresh || FloatingPointEquals(adap->topoMutateProb, ZERO_POINT_ZERO, max(1.0e-8, GARLI_FP_EPS * 2.0)))
					&& (gen > adap->intervalsToStore * adap->intervalLength)
					&& (adap->improveOverStoredIntervals < conf->improveOverStoredIntervalsThresh)
					&& (FloatingPointEquals(adap->branchOptPrecision, adap->minOptPrecision, max(1.0e-8, GARLI_FP_EPS * 2.0)) || adap->numPrecReductions==0)){
						if(adap->topoMutateProb > ZERO_POINT_ZERO) 
							outman.UserMessage("Reached termination condition!\nlast topological improvement at gen %d", lastTopoImprove);
						else 
							outman.UserMessage("Reached termination condition!\n");
						outman.UserMessage("Improvement over last %d gen = %.5f", adap->intervalsToStore*adap->intervalLength, adap->improveOverStoredIntervals);
						break;
						}
					}
				}

#ifdef INCLUDE_PERTURBATION
			CheckPerturbSerial();
#endif
			}

		if(ShouldCheckpoint(true) == true)
			WriteStateFiles();

		if(stopwatch.ThisExecutionSplitTime() > conf->stoptime){
			outman.UserMessage("NOTE: ****Specified time limit (%d seconds) reached...", conf->stoptime);
			//Time termination can be used a sort of "pause" along with checkpointing.  Checkpoints may be
			//written very infrequently though (large saveevery), so spit one out now.
			//Always do this in BOINC mode.	
#ifndef BOINC
			if(conf->checkpoint)
#endif
				WriteStateFiles();	
			timeTermination = true;
			break;
			}
		if(gen == conf->stopgen){
			//stopgen is essentially a stopping condition, treated just like the flexible automated criterion.
			outman.UserMessage("NOTE: ****Specified generation limit (%d) reached...", conf->stopgen);
			genTermination = true;
			}
#ifdef INCLUDE_PERTURBATION
		if(pertMan->pertAbandoned==true && pertMan->restartAfterAbandon==true && (gen - pertMan->lastPertGeneration > pertMan->gensBeforeRestart)){
			params->starting_tree="";
			pertMan->lastPertGeneration=gen;
			pertMan->pertAbandoned=false;
			pertMan->numPertsNoImprove=0;
			bestSinceRestart.SetFitness(-1e100);
			if(BestFitness() > allTimeBest.Fitness()) StoreAllTimeBest();
			SeedPopulationWithStartingTree();
			outman.UserMessage("restarting ....");
			}
#endif
		}
	//Allow killing during FinalOpt
	TurnOffSignalCatching();

	//checkpoint immediately before final opt, if finishedGenerations isn't set, which will
	//indicate that we've already written and restarted from a checkpoint written here
	if(!(finishedGenerations || timeTermination || userTermination)){
		//this will indicate that we've finished the loop over generations through automatic means
		//and that we are just before finalOpt (i.e., finishedGenerations is true but finishedRep is
		//false).  This will be critical for restarting from a checkpoint written just before finalOpt
		if(!(timeTermination || userTermination))
			finishedGenerations = true;

#ifndef BOINC
		//non-BOINC checkpointing
		if(ShouldCheckpoint(false))
#endif
			WriteStateFiles();
		if(conf->workPhaseDivision){
			//second workphasedivision exit
			WriteStateFiles();
			outman.UserMessage("\nNOTE: Terminating run before final optimization and writing checkpoint");
			outman.UserMessage("because workphasedivision configuration entry was set.");
			workPhaseTermination = true;
			return;
			}
		}

	//don't optimize if checkpointing is happening and the run was prematurely killed
	if(conf->refineEnd && !(conf->checkpoint && (timeTermination || userTermination))){
		UpdateFractionDone(3);
		BetterFinalOptimization();
		finishedRep = true;
		//finishedRep must be true for the following two functions to output correctly
		OutputLog();
		if(conf->outputTreelog && treeLog.is_open())
			AppendTreeToTreeLog(-1);
		}

	//outman.UserMessage("Maximum # clas used = %d out of %d", claMan->MaxUsedClas(), claMan->NumClas());
	//outman.UserMessage("%d conditional likelihood calculations\n%d branch optimization passes", calcCount, optCalcs);
	UpdateFractionDone(4);
	}


FLOAT_TYPE Population::GenerationFractionDone(){
	//This just pulls out some complicated (and ad hoc) code out of UpdateFractionDone that assigns a proportion done 
	//to a point during the search (generation) phase of a run.
	bool willReduce = (FloatingPointEquals(adap->startOptPrecision, adap->minOptPrecision, 1e-6) == false) 
		&& (adap->numPrecReductions > 0);
	//the 0.45 here is for rounding purposes.  Don't want to round down if precReductionFactor ends up being 
	//slightly more than the diff due to floating point rep
	int reduction_number = willReduce ? (int) (0.45 + ((adap->startOptPrecision - adap->branchOptPrecision) / adap->precReductionFactor)) : 0;
	int remaining_reductions = willReduce ? adap->numPrecReductions - reduction_number : 0;
	double evalInterval = (adap->intervalLength * adap->intervalsToStore);

	FLOAT_TYPE genProportionDone = 0.0;

	//these should add up to one
	const FLOAT_TYPE preReduceSplit = 0.2;
	const FLOAT_TYPE duringReduceSplit = 0.5;
	const FLOAT_TYPE postReduceSplit = 0.3;

	if(willReduce && remaining_reductions == adap->numPrecReductions){
		//We've done a decent number of gen, but haven't yet reduced the prec.  This will be linear until the minimum
		//possible number of generations before a prec reduction could happen have passed, thereafter it will be asymptotic
		FLOAT_TYPE linearProportion = 0.5;
		if(gen <= evalInterval)
			genProportionDone = preReduceSplit * linearProportion * (gen / evalInterval);
		else
			genProportionDone = preReduceSplit * (linearProportion + (1.0 - linearProportion) * (1.0 - (evalInterval / (FLOAT_TYPE) gen)));
		}

	else if(willReduce && remaining_reductions > 0){
		//during reduction phase - divide evenly among the precision reductions
		FLOAT_TYPE perReduction = duringReduceSplit / ((FLOAT_TYPE) adap->numPrecReductions - 1.0);
		//between one reduction and the next, linear, then asymptotic
		FLOAT_TYPE sinceLastReduction = gen - lastPrecisionReduction;
		FLOAT_TYPE linearProportion = 0.5;
		genProportionDone = preReduceSplit + (perReduction * (reduction_number - 1));
		if(sinceLastReduction <= evalInterval)
			//genProportionDone +=((reduction_number - 1) * perReduction) + (duringReduceSplit * perReduction * ((FLOAT_TYPE) sinceLastReduction / evalInterval));
			genProportionDone += perReduction * linearProportion * ((FLOAT_TYPE) sinceLastReduction / evalInterval);
		else
			//genProportionDone += ((reduction_number - 1) * perReduction) + (duringReduceSplit * perReduction) + (1.0 - duringReduceSplit) * perReduction * (1.0 - (evalInterval / (FLOAT_TYPE) sinceLastReduction));
			genProportionDone += perReduction * (linearProportion + (1.0 - linearProportion) * (1.0 - (evalInterval / (FLOAT_TYPE) sinceLastReduction)));
		}

	else if(remaining_reductions == 0){
		//this is linear until we get past the absolute minimum point that the run could have finished
		FLOAT_TYPE sinceLastReduction = gen - lastPrecisionReduction;
		//the chance of going over the minimum # gen in the last phase is small with lower # of taxa, which makes for a big jump
		//in proportion because the asymptotic phase isn't entered at all.  Scale the proportion where the asymp phase starts
		//with the # of taxa
		double linearProportion = max(0.5, 0.9 - 0.10 * (dataPart->NTax() / 50.0));
		genProportionDone = preReduceSplit + duringReduceSplit;
		if(sinceLastReduction <= conf->lastTopoImproveThresh){
			genProportionDone += postReduceSplit * linearProportion * (sinceLastReduction / (FLOAT_TYPE) conf->lastTopoImproveThresh);
			}
		//thereafter it is conservatively asymptotic
		else{
			assert( (1.0 - (conf->lastTopoImproveThresh / (FLOAT_TYPE) sinceLastReduction)) >= 0.0);
			genProportionDone += postReduceSplit * (linearProportion + (1.0 - linearProportion) * min(1.0, (1.0 - (conf->lastTopoImproveThresh / (FLOAT_TYPE) sinceLastReduction))));
			}
		}

	assert(genProportionDone >= 0.0 && genProportionDone <= 1.0);
	return genProportionDone;
	}

void Population::UpdateFractionDone(int phase){
#ifndef BOINC
	if(!conf->reportRunProgress)
		return;
#endif

	//Update the proportion done.  This is mainly for BOINC, but might be used elsewhere.
	//
	//-This uses the "search phase" and other info to automatically determine the progress using a fairly arbitrary algorithm.
	//-The fractions done are members of Population (rep_fraction_done and tot_fraction_done), so values will be maintained
	//	until updated, and stored in checkpoints.
	//-The fractions should never reduce during the course of a run, or after restarting from checkpoint.
	//-It is generally not clear when various phases will end, so the intent is to have a portion of each phase increase the
	//	fraction linearly, and then switch to an assymptotic approach to some value.
	//-The searchphasedivision setting complicates things a bit.  In this case the run terms immediately after initial opt,
	//	and immediately before final opt.  During the generation cycle it will generally be time limited, but doesn't need 
	//	to be.  During each of those phases it will be assumed that the whole fraction done refers to just that phase, and 
	//	that upon restarting from checkpoint the fraction done will reset, as opposed to continuing from the previous value.

	//phases:
	//0 - data has been read and everything allocated OR a search has just started a new replicate
	//1 - during initial optimization, based on optimization pass
	//2 - during generation cycle, this can be difficult, because no real a priori way to know if run will (or is intended to)
	//		stop due to stoptime, stopgen or auto termination.  So, will use the max of those.
	//		The calculations for fraction done assuming autotermination also define these phases:
	//			pre-reduction - before first precision reduction
	//			reduction - while reductions are happening
	//			terminal - remaining gens after min prec reached
	//3 - final opt, based on optimization pass
	//4 - a replicate is entirely done

	//Called from:
	//0 - if not restart at start of loop over reps in PerformSearch (reset rep_fraction_done)
	//1


	FLOAT_TYPE afterSetup = 0.05;
	FLOAT_TYPE beforeTerm = 0.95;
	FLOAT_TYPE workDuration = beforeTerm - afterSetup;

	//under normal conditions, these should add up to 1.0, and divide workDuration into sub-intervals
	FLOAT_TYPE initialOptProportionOfWork = -1.0;
	FLOAT_TYPE finalOptProportionOfWork = -1.0;
	FLOAT_TYPE genProportionOfWork = -1.0;

	if(conf->workPhaseDivision){
		//because in this mode it stops right after initial opt, and starts right before final
		initialOptProportionOfWork = finalOptProportionOfWork = genProportionOfWork = 1.0;
		}
	else{
		initialOptProportionOfWork = 0.15;
		finalOptProportionOfWork = 0.15;
		genProportionOfWork = 1.0 - initialOptProportionOfWork - finalOptProportionOfWork;
		}

	FLOAT_TYPE newFract = -1.0;
	FLOAT_TYPE newRepFract = -1.0;

	if(phase == 0){
		//This will reset for the start of a new rep.  The fraction of completion contibuted by previous finished reps
		//will be included below.
		//newFract = afterSetup;
		rep_fraction_done = newRepFract = 0.0;
		}
	else if(phase == 1){
		//during initial optimization
		FLOAT_TYPE linearUntilPass = 10;
		FLOAT_TYPE linearProportion = 0.5;
		if(initialRefinePass < 0) //this means that initial opt has completed
			newRepFract = initialOptProportionOfWork;
		else if(initialRefinePass <= linearUntilPass)
			newRepFract = initialOptProportionOfWork * linearProportion * min(1.0, (initialRefinePass / linearUntilPass));
		else
			newRepFract = initialOptProportionOfWork * (linearProportion + (1.0 - linearProportion) * max(0.0, 1.0 - linearUntilPass / (FLOAT_TYPE) initialRefinePass));
		}
	else if(phase == 2){
		//during generations
		newRepFract = conf->workPhaseDivision == true ? 0.0 : initialOptProportionOfWork;
		//No real way of knowing what will cause termination of this run.  Choose the max.
		//In the case of time term it will "reset" upon a restart from checkpoint, which might
		//cause odd changes in the fraction done in some cases, but whatever. 
		FLOAT_TYPE timeFract = min(1.0, (FLOAT_TYPE) stopwatch.ThisExecutionSplitTime() / conf->stoptime);
		FLOAT_TYPE genLimitFract = min(1.0, (FLOAT_TYPE) gen / conf->stopgen);
		FLOAT_TYPE genFract = conf->enforceTermConditions ? GenerationFractionDone() : 0.0;
		
		//An odd case can happen here, where the generation phase is nearly done (giving high fraction done) 
		//but a run is likely to be time limited, and after a restart the generation fraction calculation 
		//will still be very high.  That would make the fraction done immediately upon restart to be very high
		//and to just hang there for a long time.  So, to be conservative downweight the generation fraction if the 
		//stoptime is < 5 hr.
		if(conf->stoptime < 5 * 60 * 60)
			newRepFract += genProportionOfWork * max(max(timeFract, genLimitFract), genFract * 0.5);
		else
			newRepFract += genProportionOfWork * max(max(timeFract, genLimitFract), genFract);
		}
	else if(phase == 3){
		//during final optimization
		FLOAT_TYPE linearUntilPass = 20;
		FLOAT_TYPE linearProportion = 0.5;
		newRepFract = conf->workPhaseDivision == true ? 0.0 : initialOptProportionOfWork + genProportionOfWork;
		if(finalRefinePass <= linearUntilPass)
			newRepFract += finalOptProportionOfWork * linearProportion * min(1.0, finalRefinePass / linearUntilPass);
		else
			newRepFract += finalOptProportionOfWork * (linearProportion + (1.0 - linearProportion) * max(0.0, 1.0 - linearUntilPass / finalRefinePass));
		}
	else if(phase == 4){
		newRepFract = 1.0;
		}

	//outman.DebugMessage("newrep, oldrep: %f, %f, %d", newRepFract, rep_fraction_done, phase);
	assert(newRepFract >= 0.0 && newRepFract <= 1.0);
	assert(newRepFract + 0.00001 >= rep_fraction_done);

	//Now, prorate the newRepFract, which corresponds to the work portion of a rep and ranges from 0.0 to 1.0
	if(newFract < 0.0){
		//In searchphasereduction, don't prorate, since each segment makes up its own 0-1 fraction
		if(conf->workPhaseDivision)
			newFract = afterSetup + workDuration * newRepFract;
		else{
			int totSearches = conf->searchReps * (conf->bootstrapReps > 0 ? conf->bootstrapReps : 1);
			int curSearch = currentSearchRep + (currentBootstrapRep > 0 ? currentBootstrapRep - 1 : 0) * conf->searchReps;
			FLOAT_TYPE perRepProportion = 1.0 / totSearches;
			FLOAT_TYPE completedRepProportion = (curSearch - 1.0) * perRepProportion;
			newFract = afterSetup + workDuration * (completedRepProportion + (newRepFract * perRepProportion));
			}
		}

	//outman.DebugMessage("newtot, oldtot: %f, %f, %d", newFract, tot_fraction_done, phase);
	assert(newFract >= 0.0 && newFract <= 1.0);
	assert(newFract + 0.00001 >= tot_fraction_done);
	rep_fraction_done = newRepFract;
	tot_fraction_done = newFract;

#ifdef BOINC
	boinc_fraction_done(tot_fraction_done);
#endif
	//outman.DebugMessage("newrep, newtot: %f, %f, %d", rep_fraction_done, tot_fraction_done, phase);
	}

//this is a final opt adapted from final opt of trunk version 1.0
void Population::BetterFinalOptimization(){
	outman.setf(ios::fixed);
	outman.precision(5);
	outman.UserMessage("Current score = %.4f", BestFitness());

#ifdef INCLUDE_PERTURBATION
	if(pertMan->ratcheted) TurnOffRatchet();

	if(allTimeBest != NULL){
		if(BestFitness() < allTimeBest->Fitness()){
			RestoreAllTimeBest();
			}
		}
#endif

	//This was a little dangerous since any subsequent scoring of any of the trees would cause problems
	//probably not that important anyway.
/*	for(unsigned i=0;i<total_size;i++){
		if(i != bestIndiv) indiv[i].treeStruct->RemoveTreeFromAllClas();
		}
*/
	outman.UserMessage("Performing final optimizations...");
#ifdef MAC_FRONTEND
	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	[[MFEInterfaceClient sharedClient] didBeginBranchOptimization];
	[pool release];
#endif
	FLOAT_TYPE incr;

	assert(finalRefinePass < 1);
	finalRefinePass = 1;

	double freqOptImprove, nucRateOptImprove, pinvOptImprove, alphaOptImprove, omegaOptImprove, flexOptImprove, subRateOpt, insDelOptImprove;
	double paramOpt, blenOptImprove;
	paramOpt = blenOptImprove = freqOptImprove = nucRateOptImprove = pinvOptImprove = alphaOptImprove = omegaOptImprove = flexOptImprove = subRateOpt = insDelOptImprove = ZERO_POINT_ZERO;

	FLOAT_TYPE precThisPass = max(adap->branchOptPrecision * pow(ZERO_POINT_FIVE, finalRefinePass), (FLOAT_TYPE)1e-10);
	FLOAT_TYPE paramPrecThisPass = max(adap->branchOptPrecision*0.1, 0.01);
	bool optAnyModel = FloatingPointEquals(conf->modWeight, ZERO_POINT_ZERO, 1e-8) == false;
	bool goingToExit;

	Individual *optInd = &indiv[bestIndiv];
	Tree *optTree = optInd->treeStruct;

	adap->branchOptPrecision = min(adap->branchOptPrecision, 0.01);

	string blenS;
	do{
		//during each pass we'll keep track of a few things
		//  incr = total improvement this pass. this controls termination of opt
		//  summed improvement for each param/blens since last output.  If not outputting
		//    every pass it won't be zeroed and so will accumulate.  The output string
		//    will be constructed each pass, but only output sometimes
		precThisPass = max(adap->branchOptPrecision * pow(ZERO_POINT_FIVE, finalRefinePass), (FLOAT_TYPE)1e-10);
		paramPrecThisPass = max(precThisPass, 1e-5);

		string outString;
		char temp[50];
		FLOAT_TYPE passStart=optInd->Fitness();

		//back up the current branch lengths in case something goes wrong in the blen optimization
		vector<FLOAT_TYPE> blens;
		optTree->StoreBranchlengths(blens);

		//remember that what is returned from OptAllBranches isn't the true increase in score, just an estimate
		incr=optTree->OptimizeAllBranches(precThisPass);
		optInd->CalcFitness(0);

		FLOAT_TYPE trueImprove= optInd->Fitness() - passStart;

		//In very rare cases the score can come out very slightly worse (or apparently worse due to numerical instability issues) after
		//optimizing all of the branches.  In general this is taken care of at a lower level, but if it percolates up to here we'll ignore
		//the last set of changes and pretend they never happened.
		if(trueImprove < ZERO_POINT_ZERO){
			outman.DebugMessage("OptimizeAllBranches worsened score by %f.  Restoring previous branch lengths...", trueImprove);
			optTree->RestoreBranchlengths(blens);
			trueImprove = ZERO_POINT_ZERO;
			optInd->SetDirty();
			optInd->CalcFitness(0);
			}

		blenOptImprove += trueImprove;
		incr = trueImprove;

		sprintf(temp, "(branch= %4.4f", blenOptImprove);
		blenS = temp;

		optInd->CalcFitness(0);

		//these strings will be overwritten each time one of the parameter types are optimized
		//always with the sum total of improvement due to that param, be it over models, passes, etc.
		//this means that each will only appear once, even in partitioned models
		string omegaS, alphaS, flexS, pinvS, freqsS, relRatesS, insDelS, subsetS;
		omegaS = alphaS = flexS = pinvS = freqsS = relRatesS = insDelS = subsetS = "";
		for(int m = 0;m < indiv[bestIndiv].modPart.NumModels();m++){
			Model *mod = indiv[bestIndiv].modPart.GetModel(m);
			const ModelSpecification *modSpec = mod->GetCorrespondingSpec();

			bool optOmega, optAlpha, optFlex, optPinv, optFreqs, optRelRates, optInsDel;
			optOmega = optAlpha = optFlex = optPinv = optFreqs = optRelRates = optInsDel = false;

			if(modSpec->IsCodon() && ! modSpec->fixOmega)
				optOmega = true;
			else if(modSpec->numRateCats > 1 && !modSpec->IsCodon()){
				if(modSpec->IsFlexRateHet())
					optFlex = true;
				else if(modSpec->fixAlpha == false)
					optAlpha = true;
				}
			if(modSpec->includeInvariantSites && !modSpec->fixInvariantSites)
				optPinv = true;

			if(modSpec->IsCodon() == false && modSpec->fixStateFreqs == false && modSpec->IsEqualStateFrequencies() == false && modSpec->IsEmpiricalStateFrequencies() == false)
				optFreqs = true;
			//this is the case of forced freq optimization with codon models.  For everything to work they must be set as both not fixed but empirical
			if(modSpec->IsCodon() && modSpec->fixStateFreqs == false && modSpec->IsEqualStateFrequencies() == false && modSpec->IsEmpiricalStateFrequencies() == true)
				optFreqs = true;
			if((modSpec->fixRelativeRates == false) && ((modSpec->Nst() > 1 && modSpec->IsAminoAcid() == false) || modSpec->IsEstimateAAMatrix() || modSpec->IsTwoSerineRateMatrix()))
				optRelRates = true;
			if(modSpec->IsOrientedGap())
				optInsDel = true;

			//this is taken from the improved version in the trunk, and is a bit redundant in this context.  
			//the output strings will be generated every time that any of the params are optimized, and will
			//then be updated the next time the same parameter type is optimized in a different model.  The
			//last model to optimize a given param will make the correct string that will eventually get output
			if(optOmega) {
				paramOpt = optTree->OptimizeOmegaParameters(paramPrecThisPass, m);
				if(paramOpt < ZERO_POINT_ZERO && paramOpt > -1e-8)//avoid printing very slightly negative values
					paramOpt = ZERO_POINT_ZERO;
				omegaOptImprove += paramOpt;
				sprintf(temp, "  omega= %4.4f", omegaOptImprove);
				omegaS = temp;
				incr += paramOpt;
				}
			if(optAlpha){
				paramOpt = optTree->OptimizeBoundedParameter(m, paramPrecThisPass, mod->Alpha(), 0, min(0.05, mod->Alpha()), max(999.9, mod->Alpha()), &Model::SetAlpha);
				if(paramOpt < ZERO_POINT_ZERO && paramOpt > -1e-8)//avoid printing very slightly negative values
					paramOpt = ZERO_POINT_ZERO;
				alphaOptImprove += paramOpt;
				sprintf(temp, "  alpha= %4.4f", alphaOptImprove);
				alphaS = temp;
				incr += paramOpt;
				}
			if(optFlex){
				//Flex opt is tough, give it more passes if they are helping
				FLOAT_TYPE p = 0.0;
				paramOpt = 0.0;
				int innerPass = 0;
				do{
					p = optTree->OptimizeFlexRates(paramPrecThisPass, m);
					paramOpt += p;
					}while(p > trueImprove && innerPass++ < 5);
				if(paramOpt < ZERO_POINT_ZERO && paramOpt > -1e-8)//avoid printing very slightly negative values
					paramOpt = ZERO_POINT_ZERO;
				flexOptImprove += paramOpt;
				sprintf(temp, "  flex rates= %4.4f", flexOptImprove);
				flexS = temp;
				incr += paramOpt;
				}
			if(optPinv){
				paramOpt = optTree->OptimizeBoundedParameter(m, paramPrecThisPass, mod->PropInvar(), 0, min(1.0e-8,mod->PropInvar()), mod->maxPropInvar, &Model::SetPinv);
				if(paramOpt < ZERO_POINT_ZERO && paramOpt > -1e-8)//avoid printing very slightly negative values
					paramOpt = ZERO_POINT_ZERO;
				pinvOptImprove += paramOpt;
				sprintf(temp, "  pinv= %4.4f", pinvOptImprove);
				pinvS = temp;
				incr += paramOpt;
				}
			if(optFreqs){
				paramOpt = optTree->OptimizeEquilibriumFreqs(paramPrecThisPass, m);
				if(paramOpt < ZERO_POINT_ZERO && paramOpt > -1e-8)//avoid printing very slightly negative values
					paramOpt = ZERO_POINT_ZERO;
				freqOptImprove += paramOpt;
				sprintf(temp, "  eq freqs= %4.4f", freqOptImprove);
				freqsS = temp;
				incr += paramOpt;
				}
			if(optRelRates){
				paramOpt = optTree->OptimizeRelativeNucRates(paramPrecThisPass, m);
				if(paramOpt < ZERO_POINT_ZERO && paramOpt > -1e-8)//avoid printing very slightly negative values
					paramOpt = ZERO_POINT_ZERO;
				nucRateOptImprove += paramOpt;
				sprintf(temp, "  rel rates= %4.4f", nucRateOptImprove);
				relRatesS = temp;
				incr += paramOpt;
				}
			if(optInsDel){
				paramOpt = optTree->OptimizeInsertDeleteRates(paramPrecThisPass, m);
				if(paramOpt < ZERO_POINT_ZERO && paramOpt > -1e-8)//avoid printing very slightly negative values
					paramOpt = ZERO_POINT_ZERO;
				insDelOptImprove += paramOpt;
				sprintf(temp, "  ins/del rates= %4.4f", insDelOptImprove);
				insDelS = temp;
				incr += paramOpt;
				}
			
			optInd->CalcFitness(0);
			}
		
		if(modSpecSet.InferSubsetRates()){
			paramOpt = indiv[bestIndiv].treeStruct->OptimizeSubsetRates(max(adap->branchOptPrecision*0.1, 0.001));
			if(paramOpt < ZERO_POINT_ZERO && paramOpt > -1e-8)//avoid printing very slightly negative values
				paramOpt = ZERO_POINT_ZERO;
			subRateOpt += paramOpt;
			sprintf(temp, "  subset rates= %4.4f", subRateOpt);
			subsetS = temp;
			paramOpt += subRateOpt;
			}
		optInd->CalcFitness(0);
		
		outString = blenS + omegaS + alphaS + flexS + pinvS + freqsS + relRatesS + insDelS + subsetS;
		goingToExit = !(incr > 1.0e-5 || precThisPass > 1.0e-4 || finalRefinePass < 10);

		UpdateFractionDone(3);
		if(finalRefinePass < 20 || (finalRefinePass % 10 == 0) || goingToExit){
			if(finalRefinePass > 20 && (goingToExit || (finalRefinePass % 10 == 0)))
				outman.UserMessage(" optimization up to ...");
			outman.UserMessageNoCR("pass %-2d: %.4f   %s)", finalRefinePass, optInd->Fitness(), outString.c_str());
			if(conf->reportRunProgress)
				outman.UserMessageNoCR(" %14.2f %14.2f", 0.01 * (int) ceil(rep_fraction_done * 100), 0.01 * (int) ceil(tot_fraction_done * 100));
			outman.UserMessage("");
			paramOpt = blenOptImprove = freqOptImprove = nucRateOptImprove = pinvOptImprove = alphaOptImprove = omegaOptImprove = flexOptImprove = subRateOpt = insDelOptImprove = ZERO_POINT_ZERO;
			}

		finalRefinePass++;
		}while(!goingToExit);
	finalRefinePass = -1;
#ifdef PUSH_TO_MIN_BLEN
	double init = indiv[bestIndiv].treeStruct->lnL;
	int num = indiv[bestIndiv].treeStruct->PushBranchlengthsToMin();
	indiv[bestIndiv].treeStruct->Score();
	double aft = indiv[bestIndiv].treeStruct->lnL;
	double imp=indiv[bestIndiv].treeStruct->OptimizeAllBranches(precThisPass);
	indiv[bestIndiv].treeStruct->Score();
	double fin = indiv[bestIndiv].treeStruct->lnL;
	outman.UserMessage("Looking for minimum length branches...");
	indiv[bestIndiv].CalcFitness(0);
	outman.DebugMessage("%d branches pushed to min.\nScore after opt: %.9f\nScore after push: %.9f\nScore after reopt: %.9f", num, init, aft, fin);
#endif

	outman.UserMessage("Final score = %.4f", indiv[bestIndiv].Fitness());
	unsigned totalSecs = stopwatch.SplitTime();
	unsigned secs = totalSecs % 60;
	totalSecs -= secs;
	unsigned min = (totalSecs % 3600)/60;
	totalSecs -= min * 60;
	unsigned hours = totalSecs / 3600;
	if(conf->searchReps == currentSearchRep && (conf->bootstrapReps == 0 || conf->bootstrapReps == currentBootstrapRep ))
		outman.UserMessage("Time used = %d hours, %d minutes and %d seconds", hours, min, secs);
	else
		outman.UserMessage("Time used so far = %d hours, %d minutes and %d seconds", hours, min, secs);

	log << "Score after final optimization: " << indiv[bestIndiv].Fitness() << endl;
	//not sure how this would be done partitioned
/*
	if(modSpec.IsCodon()){
		vector<FLOAT_TYPE> sProps;
		indiv[bestIndiv].treeStruct->mod->CalcSynonymousBranchlengthProportions(sProps);
		outman.UserMessage("Proportion of branchlengths that are Synonymous: %.5f", sProps[sProps.size()-1]); 
		}
*/
#ifdef MAC_FRONTEND
	pool = [[NSAutoreleasePool alloc] init];
	[[MFEInterfaceClient sharedClient] reportFinalScore:BestFitness()];
	[pool release];
#endif

	outman.unsetf(ios::fixed);

	if(conf->outputTreelog && treeLog.is_open())
		AppendTreeToTreeLog(-1);

#ifdef ENABLE_CUSTOM_PROFILER
	char fname[100];
	sprintf(fname, "%s.profileresults.log", conf->ofprefix.c_str());
#ifdef BOINC
	char physical_name[100];
	boinc_resolve_filename(fname, physical_name, sizeof(physical_name));
	ofstream prof(physical_name);
	//MFILE prof;
	//prof.open(physical_name, "w");
#else
	ofstream prof(fname);
#endif

	/*
//FROM WRITETREEFILE
#ifdef BOINC
	char physical_name[100];
	boinc_resolve_filename(fname, physical_name, sizeof(physical_name));
	MFILE outf;
	outf.open(physical_name, "w");
#else
	ofstream outf;
	outf.open( filename.c_str() );
	outf.precision(8);
#endif
//...
#ifdef BOINC
	const char *s = trans.c_str();
	outf.write(s, sizeof(char), trans.length());
	s = str.c_str();
	outf.write(s, sizeof(char), str.length());
	theInd->treeStruct->root->MakeNewick(treeString, false, true);
	size_t len = strlen(treeString);
	outf.write(treeString, sizeof(char), len);
	str = ";\nend;\n";
	s = str.c_str();
	outf.write(s, sizeof(char), str.length());
#else
	outf << trans;
	outf << str;
	outf.setf( ios::floatfield, ios::fixed );
	outf.setf( ios::showpoint );
	theInd->treeStruct->root->MakeNewick(treeString, false, true);
	outf << treeString << ";\n";
	outf << "end;\n";
#endif	
*/
	char str[256];
	sprintf(str, "dataset: %s\tstart:%s\n", conf->datafname.c_str(), conf->streefname.c_str());
	
	prof << "dataset: " << conf->datafname << "\t" << "start: " << conf->streefname << endl;
	prof << "seed: " << conf->randseed << "\t" << "refine: " << (conf->refineStart == true) << endl;
	prof << "start prec: " << conf->startOptPrec << "\t" << "final prec: " << adap->branchOptPrecision << endl;

#ifdef SINGLE_PRECISION_FLOATS
	prof << "Single precision\n";
#else
	prof << "Double precision\n";
#endif
	unsigned s = stopwatch.SplitTime();
	prof << "Total Runtime: " << s << "\tnumgen: " << gen << "\tFinalScore: " << indiv[bestIndiv].Fitness() << "\n";
	outman.SetOutputStream(prof);
	indiv[bestIndiv].modPart.OutputHumanReadableModelReportWithParams();

	prof << "Function\t\tcalls\ttime\tTperC\t%runtime" << endl;
	ProfIntInt.Report(prof, s);
	ProfIntTerm.Report(prof, s);
	ProfTermTerm.Report(prof, s);
	ProfRescale.Report(prof, s);
	ProfScoreInt.Report(prof, s);
	ProfScoreTerm.Report(prof, s);
	ProfIntDeriv.Report(prof, s);
	ProfTermDeriv.Report(prof, s);
	ProfCalcPmat.Report(prof, s);
	ProfCalcEigen.Report(prof, s);
	ProfModDeriv.Report(prof, s);
	ProfNewton.Report(prof, s);
	ProfEQVectors.Report(prof, s);
	prof.close();
	outman.SetOutputStream(cout);
#endif
	/*	cout << "intterm calls " << inttermcalls << " time " << inttermtime/(double)(ticspersec.QuadPart) << endl;
	cout << "termterm calls " << termtermcalls << " time " << termtermtime/(double)(ticspersec.QuadPart) << endl;
	cout << "rescale calls " << rescalecalls << " time " << rescaletime/(double)(ticspersec.QuadPart) << " numrescales " << numactualrescales << endl;
	cout << "totalopt calls " << totaloptcalls << " time " << totalopttime/(double)(ticspersec.QuadPart) << endl;
	cout << "calcderiv calls " << calcderivcalls << " time " << calcderivtime/(double)(ticspersec.QuadPart) << endl;
	cout << "derivgetclas calls " << derivgetclascalls << " time " << derivgetclastime/(double)(ticspersec.QuadPart) << endl;
	cout << "derivint calls " << derivintcalls << " time " << derivinttime/(double)(ticspersec.QuadPart) << endl;
	cout << "derivterm calls " << derivtermcalls << " time " << derivtermtime/(double)(ticspersec.QuadPart) << endl;
	cout << "modderiv calls " << modderivcalls << " time " << modderivtime/(double)(ticspersec.QuadPart) << endl;
	cout << "pmat calls " << pmatcalls << " time " << pmattime/(double)(ticspersec.QuadPart) << endl;
*/	}

//this is the original partitioned final opt
void Population::FinalOptimization(){
	//DEPRECATED in favor of BetterFinalOptimization
	assert(0);

	outman.setf(ios::fixed);
	outman.precision(5);
	outman.UserMessage("Current score = %.4f", BestFitness());
	
#ifdef INCLUDE_PERTURBATION
	if(pertMan->ratcheted) TurnOffRatchet();
	
	if(allTimeBest != NULL){
		if(BestFitness() < allTimeBest->Fitness()){
			RestoreAllTimeBest();
			}
		}
#endif

	for(unsigned i=0;i<total_size;i++){
		if(i != bestIndiv) indiv[i].treeStruct->RemoveTreeFromAllClas();
		}
	
	outman.UserMessage("Performing final optimization...");
#ifdef MAC_FRONTEND
	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	[[MFEInterfaceClient sharedClient] didBeginBranchOptimization];
	[pool release];
#endif	
	int pass=1;
	FLOAT_TYPE incr;

	double paramOpt, subRateOpt, paramTot;
	paramTot = ZERO_POINT_ZERO;
	do{
		paramOpt = ZERO_POINT_ZERO;
		for(int m = 0;m < indiv[bestIndiv].modPart.NumModels();m++){
			const ModelSpecification *modSpec = indiv[bestIndiv].modPart.GetModel(m)->GetCorrespondingSpec();
			if(modSpec->IsFlexRateHet()) paramOpt = indiv[bestIndiv].treeStruct->OptimizeFlexRates(max(adap->branchOptPrecision*0.1, 0.001), m);
			else if(modSpec->IsCodon()) paramOpt = indiv[bestIndiv].treeStruct->OptimizeOmegaParameters(max(adap->branchOptPrecision*0.1, 0.001), m);
			paramTot += paramOpt;
			if(modSpec->IsFlexRateHet()){
				outman.UserMessage("Flex optimization: %f", paramTot);
				}
			else if(modSpec->IsCodon()){
				outman.UserMessage("Omega optimization: %f", paramTot);
				}
			}
		if(modSpecSet.InferSubsetRates()){
			subRateOpt = indiv[bestIndiv].treeStruct->OptimizeSubsetRates(max(adap->branchOptPrecision*0.1, 0.001));
			outman.UserMessage("Subset rate optimization: %f", subRateOpt);
			paramOpt += subRateOpt;
			}
		}while(paramOpt > ZERO_POINT_ZERO);

	do{
		incr=indiv[bestIndiv].treeStruct->OptimizeAllBranches(max(adap->branchOptPrecision * pow(ZERO_POINT_FIVE, pass), (FLOAT_TYPE)1e-10));

		indiv[bestIndiv].CalcFitness(0);
		outman.UserMessage("\tpass %d %.4f", pass++, indiv[bestIndiv].Fitness());
		}while(incr > .00001 || pass < 10);
	outman.UserMessage("Final score = %.4f", indiv[bestIndiv].Fitness());
	unsigned totalSecs = stopwatch.SplitTime();
	unsigned secs = totalSecs % 60;
	totalSecs -= secs;
	unsigned min = (totalSecs % 3600)/60;
	totalSecs -= min * 60;
	unsigned hours = totalSecs / 3600;
	if(conf->searchReps == currentSearchRep && (conf->bootstrapReps == 0 || conf->bootstrapReps == currentBootstrapRep ))
		outman.UserMessage("Time used = %d hours, %d minutes and %d seconds", hours, min, secs);
	else 
		outman.UserMessage("Time used so far = %d hours, %d minutes and %d seconds", hours, min, secs);
		
	log << "Score after final optimization: " << indiv[bestIndiv].Fitness() << endl;
#ifdef MAC_FRONTEND
	pool = [[NSAutoreleasePool alloc] init];
	[[MFEInterfaceClient sharedClient] reportFinalScore:BestFitness()];
	[pool release];
#endif	

	outman.unsetf(ios::fixed);
	
#ifdef ENABLE_CUSTOM_PROFILER
	char fname[100];
	sprintf(fname, "%s.profileresults.log", conf->ofprefix.c_str());
	ofstream prof(fname);
	prof << "dataset: " << conf->datafname << "\t" << "start: " << conf->streefname << endl;
	prof << "seed: " << conf->randseed << "\t" << "refine: " << (conf->refineStart == true) << endl;
	prof << "start prec: " << conf->startOptPrec << "\t" << "final prec: " << adap->branchOptPrecision << endl;

#ifdef SINGLE_PRECISION_FLOATS
	prof << "Single precision\n";
#else
	prof << "Double precision\n";
#endif
	unsigned s = stopwatch.SplitTime();
	prof << "Total Runtime: " << s << "\tnumgen: " << gen << "\tFinalScore: " << indiv[bestIndiv].Fitness() << "\n";
	outman.SetOutputStream(prof);
	indiv[bestIndiv].modPart.OutputHumanReadableModelReportWithParams();

	prof << "Function\t\tcalls\ttime\tTperC\t%runtime" << endl;
	ProfIntInt.Report(prof, s);
	ProfIntTerm.Report(prof, s);
	ProfTermTerm.Report(prof, s);
	ProfRescale.Report(prof, s);
	ProfScoreInt.Report(prof, s);
	ProfScoreTerm.Report(prof, s);
	ProfIntDeriv.Report(prof, s);
	ProfTermDeriv.Report(prof, s);
	ProfCalcPmat.Report(prof, s);
	ProfCalcEigen.Report(prof, s);
	ProfModDeriv.Report(prof, s);
	ProfNewton.Report(prof, s);
	ProfEQVectors.Report(prof, s);
	prof.close();
	outman.SetOutputStream(cout);
#endif
	/*	cout << "intterm calls " << inttermcalls << " time " << inttermtime/(double)(ticspersec.QuadPart) << endl;
	cout << "termterm calls " << termtermcalls << " time " << termtermtime/(double)(ticspersec.QuadPart) << endl;
	cout << "rescale calls " << rescalecalls << " time " << rescaletime/(double)(ticspersec.QuadPart) << " numrescales " << numactualrescales << endl;
	cout << "totalopt calls " << totaloptcalls << " time " << totalopttime/(double)(ticspersec.QuadPart) << endl;
	cout << "calcderiv calls " << calcderivcalls << " time " << calcderivtime/(double)(ticspersec.QuadPart) << endl;
	cout << "derivgetclas calls " << derivgetclascalls << " time " << derivgetclastime/(double)(ticspersec.QuadPart) << endl;
	cout << "derivint calls " << derivintcalls << " time " << derivinttime/(double)(ticspersec.QuadPart) << endl;
	cout << "derivterm calls " << derivtermcalls << " time " << derivtermtime/(double)(ticspersec.QuadPart) << endl;
	cout << "modderiv calls " << modderivcalls << " time " << modderivtime/(double)(ticspersec.QuadPart) << endl;
	cout << "pmat calls " << pmatcalls << " time " << pmattime/(double)(ticspersec.QuadPart) << endl;
*/	}

//figures out the best individual that has been stored and returns index, optionally summarizes the final trees/models that have been stored
int Population::EvaluateStoredTrees(bool report){
	double bestL=-FLT_MAX;
	int bestRep;
	if(report){
		outman.UserMessage("\n#######################################################\n\nCompleted %d replicate search(es) (of %d).", storedTrees.size(), conf->searchReps);
		if(conf->searchReps > 1 && (storedTrees.size() > 1))
			outman.UserMessage("\nNOTE: Unless the following output indicates that search replicates found the\n\tsame topology, you should assume that they found different topologies.");
		outman.UserMessage("Results:");
		}
	for(unsigned r=0;r<storedTrees.size();r++){
		storedTrees[r]->treeStruct->CalcBipartitions(true);	
		if(storedTrees[r]->Fitness() > bestL){
			bestL = storedTrees[r]->Fitness();
			bestRep = r;
			}
		}
	if(report){
		for(unsigned r=0;r<storedTrees.size();r++){
			unsigned r2;
			for(r2=0;r2<r;r2++){
				if(conf->collapseBranches){
/*					if(storedTrees[r]->treeStruct->IdenticalTopologyAllowingRerooting(storedTrees[r2]->treeStruct->root)
						&& storedTrees[r2]->treeStruct->IdenticalTopologyAllowingRerooting(storedTrees[r]->treeStruct->root))
						break;
*/		
					//This is where only collapsing branches upon output gets annoying.  We really want to check
					//whether the collapsed trees are the same, but we're no longer storing them.  So, generate the collapsed
					//trees and check.  A set of collapsed trees could be generated in adavance, so doing this every time is 
					//a bit of extra work
					Individual tempInd, tempInd2;
					tempInd.DuplicateIndivWithoutCLAs(storedTrees[r]);
					tempInd2.DuplicateIndivWithoutCLAs(storedTrees[r2]);
					int num = 0;
					tempInd.treeStruct->root->CollapseMinLengthBranches(num);
					tempInd2.treeStruct->root->CollapseMinLengthBranches(num);
					tempInd.treeStruct->CalcBipartitions(true);
					tempInd2.treeStruct->CalcBipartitions(true);
					if(tempInd.treeStruct->IdenticalTopologyAllowingRerooting(tempInd2.treeStruct->root)
						&& tempInd2.treeStruct->IdenticalTopologyAllowingRerooting(tempInd.treeStruct->root))
						break;
					}
				else
					if(storedTrees[r]->treeStruct->IdenticalTopologyAllowingRerooting(storedTrees[r2]->treeStruct->root)) break;
				}
			if(r == bestRep && conf->searchReps > 1) 
				outman.UserMessageNoCR("Replicate %d : %.4f (best)", r+1, storedTrees[r]->Fitness());
			else 
				outman.UserMessageNoCR("Replicate %d : %.4f       ", r+1, storedTrees[r]->Fitness());
			if(r2 < r) 
				outman.UserMessageNoCR(" (same topology as %d)", r2+1);
			if((userTermination || timeTermination) && r == storedTrees.size() - 1)
				outman.UserMessageNoCR(" (TERMINATED PREMATURELY) ", r2+1);
			outman.UserMessage("");
			}

		if(conf->searchReps > 1)
			outman.UserMessage("\nParameter estimates across search replicates:");
		else 
			outman.UserMessage("\nParameter estimates:");
		for(int part = 0;part < storedTrees[0]->modPart.NumModels();part++){
			if(storedTrees[0]->modPart.NumModels() > 1)
				outman.UserMessage("\nPartition model subset %d:", part + 1);
			Model *tree0mod = storedTrees[0]->modPart.GetModel(part);
			if(tree0mod->GetMutableParameters()->size() > 0){
				string s;
				tree0mod->FillModelOrHeaderStringForTable(s, false);
				outman.UserMessage("       %s", s.c_str());
				for(unsigned i=0;i<storedTrees.size();i++){
					storedTrees[i]->modPart.GetModel(part)->FillModelOrHeaderStringForTable(s, true);
					outman.UserMessage("rep%2d: %s", i+1, s.c_str());
					}
				if(storedTrees[0]->modPart.GetModel(part)->GetModSpec()->IsOrientedGap()){
					outman.UserMessage("\t **ins = proportion of columns that experienced an insertion");
					outman.UserMessage("\t **del = deletion rate relative to nucleotide substitution rate");
					}
				}
			else{
				outman.UserMessage("\t Model contains no estimated parameters");
				}
			}
		outman.UserMessageNoCR("\nTreelengths");
		if(modSpecSet.InferSubsetRates())
			outman.UserMessageNoCR(" and subset rate multipliers");
		outman.UserMessage(":");
		string line;
		char cStr[100];
		sprintf(cStr, " %4s ", "TL");
		line = cStr;
		if(modSpecSet.InferSubsetRates()){
			for(int d = 0;d < dataPart->NumSubsets();d++){
				char oStr[10];
				sprintf(oStr, "R(%d)", d + 1);
				sprintf(cStr, " %5s", oStr);
				line += cStr;
				}
			}
		outman.UserMessage("       %s", line.c_str());
		for(unsigned i=0;i<storedTrees.size();i++){
			sprintf(cStr, " %5.3f", storedTrees[i]->treeStruct->Treelength());
			line = cStr;
			if(modSpecSet.InferSubsetRates()){
				for(int d = 0;d < dataPart->NumSubsets();d++){
					sprintf(cStr, " %5.3f", storedTrees[i]->modPart.SubsetRate(d));
					line += cStr;
					}
				}
			outman.UserMessage("rep%2d: %s", i+1, line.c_str());
			}
		}

	bool firstEstAArmat = true;
	for(int part = 0;part < storedTrees[0]->modPart.NumModels();part++){
		const Model *tree0mod = storedTrees[0]->modPart.GetModel(part);
		const ModelSpecification *modSpec = tree0mod->GetCorrespondingSpec();

		for(unsigned i=0;i<storedTrees.size();i++){
			Model *treeImod = storedTrees[i]->modPart.GetModel(part);
			if((modSpec->IsEstimateAAMatrix() || modSpec->IsTwoSerineRateMatrix()) && conf->bootstrapReps == 0){
				string n = conf->ofprefix.c_str();
				n += ".AArmatrix.dat";
				ofstream mat;
				if(firstEstAArmat){
					mat.open(n.c_str());
					treeImod->OutputAminoAcidRMatrixMessage(mat);
					firstEstAArmat = false;
					outman.UserMessage("Estimated amino acid rate matrix/matrices saved to %s.AArmatrix.dat", conf->ofprefix.c_str());
					}
				else
					mat.open(n.c_str(), ios::app);
				treeImod->OutputAminoAcidRMatrixArray(mat, part, i);
				mat << endl;
				mat.close();
				}
			}
		}

	return bestRep;
	}

void Population::ClearStoredTrees(){
	for(vector<Individual*>::iterator it=storedTrees.begin();it!=storedTrees.end();it++){
		delete (*it)->treeStruct;
		(*it)->treeStruct=NULL;
		delete (*it);
		}
	storedTrees.clear();
	}

void Population::Bootstrap(){

	//if we're not restarting
	if(conf->restart == false) 
		currentBootstrapRep = 1;

	for( ;currentBootstrapRep<=conf->bootstrapReps;currentBootstrapRep++){
#ifdef MAC_FRONTEND
		NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
		[[MFEInterfaceClient sharedClient] didBeginBootstrapReplicate:rep];
		[pool release];
#endif
		if(conf->restart == false){
			outman.UserMessage("\nBootstrap reweighting...");
			//if this is the first rep use the bootstrapseed if one was specified,
			//or the current seed (which could have come from a specified randseed or could have been generated randomly)
			if(nextBootstrapSeed == 0){
				assert(currentBootstrapRep == 1);
				if(conf->bootstrapSeed > 0)
					nextBootstrapSeed = conf->bootstrapSeed;
				else
					nextBootstrapSeed = rnd.seed();
				}
			lastBootstrapSeed = nextBootstrapSeed;
			nextBootstrapSeed = dataPart->BootstrapReweight(lastBootstrapSeed, conf->resampleProportion);
			}
		
		PerformSearch();
		//In workphasedivision mode we could have gotten here because PerformSearch returned early after initial
		//optimization, before final optimzation or after final optimization.  
		//if(workPhaseTermination && !(currentSearchRep == conf->searchReps && (conf->bootstrapReps == 0 || currentBootstrapRep == conf->bootstrapReps)))
		if(workPhaseTermination)
			break;
		Reset();
		
		if(!userTermination && !timeTermination){

#ifdef MAC_FRONTEND
			pool = [[NSAutoreleasePool alloc] init];
			[[MFEInterfaceClient sharedClient] didCompleteBoostrapReplicate:rep];
			[pool release];
#endif		
			}
		else {
			if(userTermination)
				outman.UserMessage("abandoning bootstrap rep %d.... terminating\n", currentBootstrapRep);
			break;
			}
		}
	}

/* OLD VERSION
void Population::Bootstrap(){
	
	data->ReserveOriginalCounts();

	stopwatch.Start();
	CatchInterrupt();

	for(int rep=1;rep <= (int) conf->bootstrapReps;rep++){
		lastTopoImprove = lastPrecisionReduction = gen = 0;
		outman.UserMessage("bootstrap replicate %d (seed %d)", rep, rnd.seed());
#ifdef MAC_FRONTEND
		NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
		[[MFEInterfaceClient sharedClient] didBeginBootstrapReplicate:rep];
		[pool release];
#endif				
		data->BootstrapReweight();
		
		SeedPopulationWithStartingTree();
		Run();
		
		if(userTermination == false){
			adap->branchOptPrecision = adap->startOptPrecision;
			FinishBootstrapRep(rep);
			outman.UserMessage("finished with bootstrap rep %d\n", rep);
#ifdef MAC_FRONTEND
			pool = [[NSAutoreleasePool alloc] init];
			[[MFEInterfaceClient sharedClient] didCompleteBoostrapReplicate:rep];
			[pool release];
#endif		
			}
		else {
			outman.UserMessage("abandoning bootstrap rep %d ....terminating", rep);
			break;
			}
		}
	FinalizeOutputStreams();
	}
*/

//this function manages multiple search replicates, setting up the population
//and then calling Run().  It can be called either directly from main(), or 
//from Bootstrap()
void Population::PerformSearch(){
	if(conf->restart == false)
		currentSearchRep = 1;
	else{
		outman.UserMessage("\nRestarting from checkpoint...");
		if(finishedRep == true){
			//if we've restarted but the last checkpoint written apparently represents 
			//the state of the population immediately after the completion of a replicate
			currentSearchRep++;
			if(currentSearchRep > conf->searchReps && (conf->bootstrapReps == 0 || currentBootstrapRep == conf->bootstrapReps)){
				outman.UserMessage("The checkpoint loaded indicates that this run already completed.\nTo start a new run set restart to 0 and change the output\nfile prefix (ofprefix).");
				restartedAfterTermination = true;
				return;
				}
			else{//we need to initialize the output here, while the population still knows that this was a restart (before calling Reset)
				InitializeOutputStreams();
				Reset();
				}
			}
		else
			InitializeOutputStreams();
		}

	for(;currentSearchRep<=conf->searchReps;currentSearchRep++){
		string s;
		if(conf->restart == false){
			//this will reset the rep_fraction_done at the start of each rep
			UpdateFractionDone(0);
			if(currentSearchRep > 1){
				Reset();
				//this just changes what the rng has stored as the init seed ix0, for output purposes
				rnd.set_seed(rnd.seed());
				}
			}
		
		//ensure that the user can ctrl-c kill the program during creation of each stepwise addition tree
		TurnOffSignalCatching();

		GetRepNums(s);

		if(conf->restart == false){
			//the fraction done is set to 1% here, indicating the this rep is ready to go
			//if we restarted, the fraction should already have been set when reading the state files
			if(s.length() > 0) 
				outman.UserMessage("\n>>>%s<<<", s.c_str());
			SeedPopulationWithStartingTree(currentSearchRep);
			//can't initialize output until after the pop is seeded, unless restarting (which happens above)
			InitializeOutputStreams();
			//write a checkpoint, since the refinement (and maybe making a stepwise tree) could have taken a good while
#ifndef BOINC
			//non-BOINC checkpointing
			if(ShouldCheckpoint(false) == true)
#endif
				WriteStateFiles();
			if(conf->workPhaseDivision){
				WriteStateFiles();
				//first workphasedivision exit point
				outman.UserMessage("\nNOTE: Terminating run after initial optimization and writing");
				outman.UserMessage("checkpoint because workphasedivision configuration entry was set.");
				workPhaseTermination = true;
				UpdateFractionDone(4);
				break;
				}
			UpdateFractionDone(1);
			}
		else{
			//in this case the progress should have been read from checkpoint and will maintain that value going forward
			adap->SetChangeableVariablesFromConfAfterReadingCheckpoint(conf);
			if(currentSearchRep > conf->searchReps) 
				throw ErrorException("rep number in checkpoint (%d) is larger than total rep specified in config (%d)", currentSearchRep, conf->searchReps);
			outman.UserMessage("%s generation %d, seed %d, best lnL %.3f", s.c_str(), gen, rnd.init_seed(), indiv[bestIndiv].Fitness());
			}

#ifndef BOINC
		//Start catching Ctrl-C's
		TurnOnSignalCatching();
#endif				
		if(!conf->scoreOnly)
			Run();

		//for most purposes, these two types of termination are premature and treated identically
		//gen termination is treated as normal termination besides some warnings
		//bool prematureTermination = (userTermination || timeTermination);
		bool prematureTermination = conf->checkpoint ? genTermination : (genTermination | timeTermination || userTermination);

		//if we're checkpointing and terminated prematurely just bail without doing anything else
		if(conf->checkpoint && (userTermination || timeTermination || workPhaseTermination)){
#ifndef BOINC
			if(!workPhaseTermination){
				outman.UserMessage("\nNOTE: A CHECKPOINTED RUN (writecheckpoints = 1) WAS PREMATURELY");
				outman.UserMessage("TERMINATED.  OUTPUT FILES (tree files, etc.) WILL NOT BE" ); 
				outman.UserMessage("FINALIZED SO THAT THE RUN CAN BE RESTARTED WHERE IT LEFT OFF");
				outman.UserMessage("(set restart = 1 in the config file).  IF YOU WANT TO USE THE");
				outman.UserMessage("PARTIAL OUTPUT FILES WITHOUT RESTARTING YOU WILL NEED TO MANUALLY");
				outman.UserMessage("ADD \"end;\" TO THE TREE FILES.\n");
				}
#endif
			UpdateFractionDone(4);
			break;
			}

		outman.UserMessage("");
		if(userTermination)
			outman.UserMessage("MODEL REPORT - SEARCH TERMINATED BY USER");
		else if(timeTermination)
			outman.UserMessage("MODEL REPORT - SEARCH TERMINATED AFTER REACHING TIME LIMIT");
		else if(genTermination)
			outman.UserMessage("MODEL REPORT - SEARCH TERMINATED AFTER REACHING GENERATION LIMIT");
		else
			outman.UserMessage("MODEL REPORT - Parameter values are FINAL");
		indiv[bestIndiv].modPart.OutputHumanReadableModelReportWithParams();

		//this rep is over
		//11/28/09 We will now always store the final individual in the stored trees array, 
		//even if prematureTerm
		if(Tree::outgroup != NULL) 
			OutgroupRoot(&indiv[bestIndiv], bestIndiv);
		//this individual will be stored in the storedTrees array until population deletes it much later,
		Individual *repResult = new Individual(&indiv[bestIndiv]);
		//Note that the collapsed individual is intentionally not stored here.  It will be re-collapsed on
		//output to file, and the collapsing here is just for this message
		if(conf->collapseBranches){
			Individual repResultColl(&indiv[bestIndiv]);
			int numCollapsed = 0;
			repResultColl.treeStruct->root->CollapseMinLengthBranches(numCollapsed);
			outman.UserMessage("NOTE: Collapsing of minimum length branches was requested (collapsebranches = 1)");\
			if(numCollapsed == 0)
				outman.UserMessage("    No branches were short enough to be collapsed.\n");
			else 
				outman.UserMessage("    %d branches were collapsed.\n", numCollapsed);
			if(repResult->treeStruct->constraints.empty() == false){
				for(vector<Constraint>::iterator con=repResult->treeStruct->constraints.begin();con!=repResult->treeStruct->constraints.end();con++){
					if(con->IsPositive()){
						outman.UserMessage("\nNOTE: If collapsing of minimum length branches is requested (collapsebranches = 1) in a run with\n\ta positive constraint, it is possible for a constrained branch itself to be collapsed.\n\tIf you care, be careful to check whether this has happened or turn off branch collapsing.\n");
						break;
						}
					}
				}
			}
		storedTrees.push_back(repResult);

		//output site likelihoods if requested
		if(conf->outputSitelikelihoods > 0){
			outman.UserMessage("Saving site likelihoods to file %s.sitelikes.log ...", conf->ofprefix.c_str());
			for( int set = 0;set < dataPart->NumSubsets();set++){ 
				if( dataPart->GetSubset(set)->DidUseDefaultWeightsets() == true)
					outman.UserMessage("WARNING: Site likelihoods are being output when wtset %s is in effect.\n Sites with weight > 1 will only be output once!", dataPart->GetSubset(set)->WeightsetName().c_str());
				}
			
			if( (userTermination || timeTermination || genTermination) ){
				outman.UserMessage("WARNING: Site likelihoods being output on prematurely terminated search ...");
				}

			//This has to work somewhat differently for partitioning.  As far as the tree functions know we will always
			//be in append mode for the sitelike output (negative sitelike value).  The pop will have to nuke any 
			//existing file here the first time through and put in the header
			indiv[bestIndiv].treeStruct->sitelikeLevel = -(int) conf->outputSitelikelihoods;

			ofstream ordered;
			indiv[bestIndiv].treeStruct->ofprefix = conf->ofprefix;
			string oname = indiv[bestIndiv].treeStruct->ofprefix + ".sitelikes.log";
			if(currentSearchRep == 1){
				ordered.open(oname.c_str());
				ordered << "Tree\t-lnL\tSite\t-lnL";
				if(conf->outputSitelikelihoods > 1) 
					ordered << "\tunder1\tunder2";
				ordered << "\n";
				ordered.close();
				}
	
			indiv[bestIndiv].treeStruct->Score();
			ordered.open(oname.c_str(), ios::app);
			ordered.precision(12);
			ordered << currentSearchRep << "\t" << -indiv[bestIndiv].treeStruct->lnL << "\n";
			ordered.close();
			}

		//warn if the normal auto-term conditions weren't used
		if(userTermination || timeTermination || genTermination){
			if(s.length() > 0 && (userTermination || timeTermination))
				outman.UserMessage(">>>Terminated %s<<<", s.c_str());
			outman.UserMessage("%s", TerminationWarningMessage().c_str());
			}
		else{
			//I think that this should only be reported here if there is > 1 search rep per boot rep, since it should really be noting that
			//a given rep has finished, and the overall boot rep doesn't really finish until after the summary across search reps
			if(s.length() > 0 && (conf->bootstrapReps == 0 || (conf->bootstrapReps > 0 && conf->searchReps > 1)))
				outman.UserMessage(">>>Completed %s<<<", s.c_str());
			}

		int best=0;
		//If this is the last search of a run, bootstrap rep, or it was killed prematurely
		//Note that EvaluateStoredTrees does some summary output for the model estimates from those trees
		if((currentSearchRep == conf->searchReps) || prematureTermination){
			if(storedTrees.size() > 0){
				best=EvaluateStoredTrees(true);
				//recombine final trees
	/*			if(total_size > 2){
					for(int i=0;i<storedTrees.size();i++)
						if(storedTrees[i]->treeStruct->root->claIndexDown == -1)
							storedTrees[i]->treeStruct->AssignCLAsFromMaster();
					for(int i=0;i<this->total_size;i++){
						//only the best indiv has clas assigned at this point
						if(i != bestIndiv) indiv[i].CopySecByRearrangingNodesOfFirst(indiv[i].treeStruct, storedTrees[best], false);
						else indiv[i].CopySecByRearrangingNodesOfFirst(indiv[i].treeStruct, storedTrees[best], true);
						}
					int holdover = 2;
					for(int rounds=0;rounds<10;rounds++){
						for(int i=holdover;i<total_size;i++){
							int partner = rnd.random_int(storedTrees.size());
							indiv[i].CrossOverWith(*storedTrees[partner], 0.01);
							}
						this->CalcAverageFitness();
						for(int i=holdover;i<total_size;i++){
							if(indiv[i].Fitness() > indiv[0].Fitness()){
								indiv[1].CopySecByRearrangingNodesOfFirst(indiv[1].treeStruct, &indiv[0], true);
								indiv[0].CopySecByRearrangingNodesOfFirst(indiv[0].treeStruct, &indiv[i], true);
								}
							else if(indiv[i].Fitness() > indiv[1].Fitness()){
								indiv[1].CopySecByRearrangingNodesOfFirst(indiv[1].treeStruct, &indiv[i], true);
								}
							}
						for(int i=holdover;i<total_size;i++){
							int parent = rnd.random_int(holdover);
							indiv[i].CopySecByRearrangingNodesOfFirst(indiv[i].treeStruct, &indiv[parent], true);
							}
						this->CalcAverageFitness();
						}
					}
				string s = "recom.";
				s += besttreefile;
				this->WriteTreeFile(s.c_str());
				Individual *repResult = new Individual(&indiv[0]);
				storedTrees.push_back(repResult);
				outman.UserMessage("Best topology created by recombination: %f", indiv[0].Fitness());
				
	*/			}
			}

		//write the best trees from all completed reps:
		//at the end of each rep
		//at the end of all reps
		//if termination was premature and we're told to write in that case (the premature tree will be included)
		if( ( (! prematureTermination) && (all_best_output & WRITE_REP_TERM)) ||
			( (! prematureTermination) && (currentSearchRep == conf->searchReps) && (all_best_output & WRITE_REPSET_TERM)) ||
			( (prematureTermination) && (storedTrees.size() > 1) && (all_best_output & WRITE_PREMATURE))){
			if(storedTrees.size() > 0){
				if(prematureTermination || currentSearchRep == conf->searchReps)//message only if last
					outman.UserMessage("\nSaving final trees from all search reps to %s.all.tre", besttreefile.c_str());
				WriteStoredTrees(besttreefile.c_str());
				}
			}

		//write the best overall tree:
		//at the end of each rep
		//at the end of all reps
		////if termination was premature and we're told to write in that case (the premature tree will be written if it is best)
		if( ( (! prematureTermination) && (best_output & WRITE_REP_TERM)) ||
			( (! prematureTermination) && (currentSearchRep == conf->searchReps) && (best_output & WRITE_REPSET_TERM)) ||
			( (prematureTermination) && (best_output & WRITE_PREMATURE))){
			//the first two options here write trees from the storedTrees array, the last writes the best from the current population
			outman.UserMessage("\nSaving final tree from best search rep (#%d) to %s.tre", best + 1, besttreefile.c_str());
			if(conf->searchReps > 1 && storedTrees.size() > 0){
				WriteTreeFile(besttreefile.c_str(), best, conf->collapseBranches);
				}
			else if(storedTrees.size() == 1)
				WriteTreeFile(besttreefile.c_str(), 0, conf->collapseBranches);
			else
				WriteTreeFile(besttreefile.c_str(), -1, conf->collapseBranches);
			}
		
		if(conf->bootstrapReps > 0){
			//write best boot tree if:
			//end of single search rep of many, and we're supposed to (not normal)
			//end of search rep set or single search
			//premature termination and we're told to
			//premature termination and we've already stored a tree (due to change, even a single termed run will now be in here)
			if( ( (! prematureTermination) && (bootlog_output & WRITE_REP_TERM)) ||
				( (! prematureTermination) && (currentSearchRep == conf->searchReps) && (bootlog_output & WRITE_REPSET_TERM)) ||
				( (prematureTermination) && (bootlog_output & WRITE_PREMATURE)) ||
				( (prematureTermination) && storedTrees.size() > 0)){
				if(conf->searchReps > 1 && storedTrees.size() > 0){
					//we're doing multiple searches per boot rep, and have successfully completed at least one replicate
					//(although the present replicate could have been prematurely terminated)
					char temp_buf[100];
					char suffix[100];
					sprintf(suffix, "boot.tre");
					DetermineFilename(bootlog_output, temp_buf, suffix);
					outman.UserMessage("\nSaving tree from best search rep (#%d) to bootstrap file %s\n", best+1, temp_buf);
					if(prematureTermination && best == storedTrees.size() - 1)
						outman.UserMessage("WARNING: Tree from prematurely terminated search saved to bootstrap file");
					FinishBootstrapRep(storedTrees[best], currentBootstrapRep);
					}

				else if(storedTrees.size() == 1){
					//We just successfully completed a one-search-rep bootstrap replicate 
					FinishBootstrapRep(storedTrees[0], currentBootstrapRep);	
					char temp_buf[100];
					char suffix[100];
					sprintf(suffix, "boot.tre");
					DetermineFilename(bootlog_output, temp_buf, suffix);
					outman.UserMessage("\nSaving best tree to bootstrap file %s\n", temp_buf);
					if(prematureTermination)
						outman.UserMessage("WARNING: Tree from prematurely terminated search saved to bootstrap file");
					}
				else //This rep was prematurely killed, but we're supposed to write it
					FinishBootstrapRep(&indiv[bestIndiv], currentBootstrapRep);
				if(!prematureTermination){
					outman.UserMessage(">>>Completed Bootstrap rep %d<<<", currentBootstrapRep);
					}
				}
			else{ 
				if(prematureTermination && !(bootlog_output & WRITE_PREMATURE)) 
					outman.UserMessage("Not saving search rep to bootstrap file due to early termination");
				}
			}

		if(conf->inferInternalStateProbs == true){
			//don't infer internals states unless at least one rep successfully completed
			if((prematureTermination == false && currentSearchRep == conf->searchReps) || (prematureTermination && storedTrees.size() > 0)){
				//this is important to ensure that there are enough free clas for a temp set to be allocated,
				//since recycling won't happen in this usage (it could if implemented, but I don't see the benefit)
				for(int i = 0;i < total_size;i++){
					if(indiv[i].treeStruct != NULL)
						indiv[i].treeStruct->MakeAllNodesDirty();
					if(newindiv[i].treeStruct != NULL)
						newindiv[i].treeStruct->MakeAllNodesDirty();
					}
				if(storedTrees.size() > 0){//careful here, the trees in the storedTrees array don't have clas assigned
					outman.UserMessage("Inferring internal state probabilities on best tree... saving to file %s.internalstates.log\n", conf->ofprefix.c_str());
					Individual *theInd;
					Individual tempInd;
					if(Tree::outgroup != NULL){
						tempInd.DuplicateIndivWithoutCLAs(storedTrees[best]);
						OutgroupRoot(&tempInd, -1);
						theInd = &tempInd;
						}
					else
						theInd = storedTrees[best];
					//InferAllInternalStateProbs will deal with assigning clas, since neither the tree in storedTrees nor the potentially temp tree have them
#ifndef USE_BEAGLE
					theInd->treeStruct->InferAllInternalStateProbs(conf->ofprefix.c_str());
					if(prematureTermination && best == storedTrees.size() - 1)
						outman.UserMessage("WARNING: Internal states inferred on tree from prematurely terminated search\n");
#else
					outman.UserMessage("**Internal state inference must be updated for beagle");
#endif
					}
				}
			else if(prematureTermination){
				outman.UserMessage(">>>Internal state probabilities not inferred due to premature termination<<<\n");
				}
			}
		//finalize anything that needs it at rep end
		FinalizeOutputStreams(0);
		//finalize anything that needs it at the end of the repset
		if(currentSearchRep == conf->searchReps || prematureTermination) {
			FinalizeOutputStreams(1);
			outman.UserMessage("#######################################################");
			}

		if(userTermination == true || timeTermination == true) 
			break;

#ifndef BOINC
		if(ShouldCheckpoint(false) == true || conf->workPhaseDivision )
#endif
			//write a checkpoint that will indicate that the rep is done and results have been written to file
			//the gen will be UINT_MAX, as it is after a rep has terminated, which will tell the function that reads
			//the checkpoint to set finishedrep = true.  This automatically happens in the BOINC case
			WriteStateFiles();

		if(conf->workPhaseDivision && !(currentSearchRep == conf->searchReps && (conf->bootstrapReps == 0 || currentBootstrapRep == conf->bootstrapReps))){
			//third workphasedivision exit point - if this is the end of the whole run, don't do this.
			outman.UserMessage("\nNOTE: Terminating run after final optimization and writing checkpoint");
			outman.UserMessage("because workphasedivision configuration entry was set.");
			workPhaseTermination = true;
			UpdateFractionDone(4);
			break;
			}
		//this needs to be set here so that the population is reset at the top of this loop before the next rep
		conf->restart = false;
		}

	ClearStoredTrees();
	}

void Population::OptimizeInputAndWriteSitelikelihoods(){
	log_output = fate_output = swaplog_output = treelog_output = problog_output = Population::DONT_OUTPUT;
	InitializeOutputStreams();

	//find out how many trees we have
	GarliReader & reader = GarliReader::GetInstance();
	const NxsTreesBlock *treesblock = reader.GetTreesBlock(reader.GetTaxaBlock(0), reader.GetNumTreesBlocks(reader.GetTaxaBlock(0)) - 1);
	if(treesblock == NULL || !strcmp(conf->streefname.c_str(), "random") || !strcmp(conf->streefname.c_str(), "stepwise"))
		throw ErrorException("You must specify a nexus treefile to use this runmode.");
	int numTrees = treesblock->GetNumTrees();

	string oname = conf->ofprefix + ".sitelikes.log";
	ofstream ordered;
	ordered.open(oname.c_str());
	ordered << "Tree\t-lnL\tSite\t-lnL\n";
	ordered.close();

	bestIndiv = 0;
	conf->searchReps = numTrees;
	//loop over the trees
	for(int t = 1;t <= numTrees;t++){
		currentSearchRep = t;
		if(!conf->scoreOnly){
			outman.UserMessage("Optimizing tree %d ...", t);
			SeedPopulationWithStartingTree(t);
			bestIndiv = 0;
			BetterFinalOptimization();
			}
		else
			outman.UserMessage("Scoring tree %d ...", t);

		outman.UserMessage("Writing site likelihoods for tree %d ...", t);
		indiv[0].treeStruct->sitelikeLevel = - (max((int)conf->outputSitelikelihoods, 1));
		indiv[0].treeStruct->ofprefix = conf->ofprefix;
		indiv[0].treeStruct->Score();
		
		ordered.open(oname.c_str(), ios::app);
		ordered.precision(10);
		ordered << t << "\t" << -indiv[0].treeStruct->lnL << "\n";
		ordered.close();

		Individual *repResult = new Individual(&indiv[0]);
		storedTrees.push_back(repResult);
		Reset();
		}
	bool coll = conf->collapseBranches;
	conf->collapseBranches = false;
	EvaluateStoredTrees(true);
	if(coll)
		outman.UserMessage("\nNOTE: collapsebranches setting ignored when writing and comparing optimized trees...");
	outman.UserMessage("\nWriting optimized trees and models to %s.all.tre", besttreefile.c_str());
	WriteStoredTrees(besttreefile.c_str());
	FinalizeOutputStreams(0);
	FinalizeOutputStreams(1);
	FinalizeOutputStreams(2);
	}

void Population::OptimizeInputAndWriteSitelikelihoodsAndTryRootings(){
	log_output = fate_output = swaplog_output = treelog_output = problog_output = Population::DONT_OUTPUT;
	//log_output = fate_output = swaplog_output = problog_output = Population::DONT_OUTPUT;
	InitializeOutputStreams();	
	//assert(Tree::someOrientedGap);
	//find out how many trees we have
	GarliReader & reader = GarliReader::GetInstance();
	const NxsTreesBlock *treesblock = reader.GetTreesBlock(reader.GetTaxaBlock(0), reader.GetNumTreesBlocks(reader.GetTaxaBlock(0)) - 1);
	if(treesblock == NULL || treesblock->GetNumTrees() > 1)
		throw ErrorException("You must specify a treefile with exactly one tree to use this runmode.");
//	int numTrees = treesblock->GetNumTrees();
//	assert(numTrees == 1);

	bestIndiv = 0;

	//start the sitelike file
	string oname = conf->ofprefix + ".sitelikes.log";
	ofstream ordered;
	ordered.open(oname.c_str());
	ordered << "Tree\t-lnL\tSite\t-lnL\n";
	ordered.close();

	Tree::useOptBoundedForBlen = true;

	currentSearchRep = 1;
	outman.UserMessage("Optimizing tree %d ...", 1);

	conf->refineStart = false;
	SeedPopulationWithStartingTree(currentSearchRep);
	bestIndiv = 0;

	assert(indiv[0].treeStruct->dummyRoot);
	//the number of branches on which the root could be attached
	conf->searchReps = (dataPart->NTax() - 1) * 2 - 3;

	double initIns = 0.05, initDel = 0.1;

	for(int m = 0;m < indiv[0].modPart.NumModels();m++){
		if(indiv[0].modPart.GetModel(m)->IsOrientedGap()){
			initIns = indiv[0].modPart.GetModel(m)->InsertRate();
			initDel = indiv[0].modPart.GetModel(m)->DeleteRate();
			}
		}
#ifdef OPT_BOUNDED_LOG
	char name[50];
	sprintf(name, "%s.optbounded.log", conf->ofprefix.c_str());
	ofstream log(name);
	log.close();
	indiv[0].treeStruct->ofprefix = conf->ofprefix;
	indiv[1].treeStruct->ofprefix = conf->ofprefix;
#endif
	BetterFinalOptimization();
	//this will stick the current tree into the treelog
	conf->outputTreelog = true;
	InitializeOutputStreams();

	outman.UserMessage("Writing site likelihoods for tree %d ...", 1);
	indiv[0].treeStruct->sitelikeLevel = -1;
	indiv[0].treeStruct->ofprefix = conf->ofprefix;
	indiv[0].treeStruct->Score();
	
	//put the score of the initial indiv in the file
	ordered.open(oname.c_str(), ios::app);
	ordered.precision(10);
	ordered << "0\t" << -indiv[0].treeStruct->lnL << "\t";
	indiv[0].treeStruct->root->MakeNewick(treeString, false, true, false);
	ordered << treeString << "\n";
	ordered.close();

	//store the indiv 
	Individual *repResult = new Individual(&indiv[0]);
	storedTrees.push_back(repResult);

	//return the rates to initial vals to have them reopt each time
		bool resetInsDel = false;
	if(resetInsDel){
		for(int m = 0;m < indiv[0].modPart.NumModels();m++){
			if(indiv[0].modPart.GetModel(m)->IsOrientedGap()){
				indiv[0].modPart.GetModel(m)->SetInsertRate(0, initIns);
				indiv[0].modPart.GetModel(m)->SetDeleteRate(0, initDel);
				indiv[0].treeStruct->MakeAllNodesDirty();
				}
			}
		}
	else
		indiv[0].treeStruct->MakeAllNodesDirty();
	indiv[0].CalcFitness(0);

	//copy the tree and model into indiv[1]
	outman.UserMessage("Rooting at nodes across tree...");
	Tree *indiv1Tree = indiv[1].treeStruct;
	indiv[1].CopySecByRearrangingNodesOfFirst(indiv1Tree, &indiv[0]);

	//get ready to swap on indiv[1]
	indiv1Tree->GatherValidReconnectionNodes(99999, indiv1Tree->dummyRoot, NULL);
	indiv1Tree->sprRang.SortByDist();

	int tnum = 2;

	for(listIt broken = indiv1Tree->sprRang.begin();broken != indiv1Tree->sprRang.end();broken++){
		//try a reattachment point for the dummy root taxon
		indiv1Tree->SPRMutate(indiv1Tree->dataPart->NTax(), &(*broken), 0.01, 0);
		//optimize the result
		bestIndiv = 1;
		indiv[bestIndiv].CalcFitness(0);
		bestFitness = indiv[bestIndiv].Fitness();

		//this will make the various trees have different names when they are appended to the treelog
		//in BetterFinalOpt
		gen = tnum;
		//AppendTreeToTreeLog(0, 1);
		
		BetterFinalOptimization();

		outman.UserMessage("%d\tnode\t%d\tlnL\t%f", tnum, (*broken).nodeNum, indiv1Tree->lnL);

		//output the sitelikes
		indiv1Tree->ofprefix = conf->ofprefix;
		indiv1Tree->sitelikeLevel = -1;
		indiv1Tree->Score();

		//add the total score and the tree
		ordered.open(oname.c_str(), ios::app);
		ordered.precision(10);
		ordered << tnum << "\t" << -indiv1Tree->lnL << "\t";
		indiv[1].treeStruct->root->MakeNewick(treeString, false, true, false);
		ordered << treeString << "\n";
		ordered.close();

		//store the indiv and write the tree to file
		repResult = new Individual(&indiv[1]);
		storedTrees.push_back(repResult);

		indiv[1].CopySecByRearrangingNodesOfFirst(indiv1Tree, &indiv[0], true);
		tnum++;
		}
	bool coll = conf->collapseBranches;
	conf->collapseBranches = false;
	EvaluateStoredTrees(true);
	if(coll)
		outman.UserMessage("\nNOTE: collapsebranches setting ignored when writing and comparing optimized trees...");
	outman.UserMessage("\nWriting optimized trees and models to %s.all.tre", besttreefile.c_str());
	WriteStoredTrees(besttreefile.c_str());
	FinalizeOutputStreams(0);
	FinalizeOutputStreams(1);
	FinalizeOutputStreams(2);
	}

void Population::VariableStartingTreeOptimization(bool reducing){
	currentSearchRep = 1;
	SeedPopulationWithStartingTree(currentSearchRep);
	InitializeOutputStreams();
	
	string filename = conf->ofprefix + ".var.log";
	ofstream out(filename.c_str());
	out.precision(10);

	filename = conf->ofprefix + ".randblens.tre";
	ofstream randTrees(filename.c_str());
	dataPart->BeginNexusTreesBlock(randTrees);

	filename = conf->ofprefix + ".optblens.tre";
	ofstream optTrees(filename.c_str());
	dataPart->BeginNexusTreesBlock(optTrees);

	typedef vector<double> doubvec;
	//this is a vector of vectors, with each entry in the higher level vector being a vector
	//with all of the final rep scores for a given precision 
	vector<doubvec> finalScores;

	typedef vector<int> intvec;
	//vector of vectors, number of passes per rep per prec
	vector<intvec> numPasses;

	vector<intvec> numDerivCalcs;

	//a triple vector, with the branch lengths for each branch, rep and prec
	typedef vector<doubvec> doubdoubvec;
	vector<doubdoubvec> allBlens;

	vector<double> prec;

	//get the precision values to use from the arbitrarystring entry in the config file
	stringstream s;
	s.str(conf->arbitraryString);
	string p;
	while(!s.eof()){
		s >> p;
		double x = atof(p.c_str());
		prec.push_back(x);
		}

	int numReps = conf->searchReps;
	double prec1;
	int numNodes = indiv[0].treeStruct->getNumNodesTotal();
	for(int rep = 0;rep < numReps;rep++){
		//for each rep, rerandomize the branch lengths and output the tree
		indiv[0].treeStruct->RandomizeBranchLengthsExponential(conf->gammaShapeBrlen);
		indiv[0].treeStruct->root->MakeNewick(treeString, false, true, false);
		randTrees << "tree r" << rep << " = [&U] " << treeString << ";\n";
		//store this randomization
		indiv[1].CopySecByRearrangingNodesOfFirst(indiv[1].treeStruct, &indiv[0], true);

		indiv[0].SetDirty();
		indiv[0].CalcFitness(0);
		double imp = 999.9;
		double prevScore = indiv[0].Fitness();
//		outman.UserMessage("%f\t%f\t", prec[p], indiv[0].Fitness());

		for(int precNum=0;precNum < prec.size() && (!reducing || (reducing && precNum < 1)) ;precNum++){
			prec1 = prec[precNum];

			int pass=0;
			prevScore = indiv[0].Fitness();
			outman.UserMessage("%f\t%d", indiv[0].Fitness(), pass);
			do{
				if(reducing) prec1 = prec[min(pass, (int)prec.size()-1)];
				indiv[0].treeStruct->OptimizeAllBranches(prec1);
				indiv[0].SetDirty();
				indiv[0].CalcFitness(0);
//				indiv[0].treeStruct->OptimizeTreeScale(prec1);

				imp = indiv[0].Fitness() - prevScore;
				prevScore = indiv[0].Fitness();
				outman.UserMessage("%f\t%f\t%d", indiv[0].Fitness(), prec1, optCalcs);
				pass++;
				}while(imp > prec1 || pass < prec.size());
			outman.UserMessage("%f\t%d\n", indiv[0].Fitness(), pass);
			if(rep == 0){
				doubvec scoreTemp;
				scoreTemp.push_back(indiv[0].Fitness());
				finalScores.push_back(scoreTemp);
				intvec passTemp;
				passTemp.push_back(pass);
				numPasses.push_back(passTemp);
				intvec calcsTemp;
				calcsTemp.push_back(optCalcs);
				numDerivCalcs.push_back(calcsTemp);
				doubvec tempBlens;
				doubdoubvec tempBlens2;
				for(int b=1;b<numNodes;b++) tempBlens.push_back(indiv[0].treeStruct->allNodes[b]->dlen);
				tempBlens2.push_back(tempBlens);
				allBlens.push_back(tempBlens2);
				}
			else{
				finalScores[precNum].push_back(indiv[0].Fitness());
				numPasses[precNum].push_back(pass);
				numDerivCalcs[precNum].push_back(optCalcs);
				doubvec tempBlens;
				for(int b=1;b<numNodes;b++) tempBlens.push_back(indiv[0].treeStruct->allNodes[b]->dlen);
				allBlens[precNum].push_back(tempBlens);
				}
			indiv[0].treeStruct->root->MakeNewick(treeString, false, true, false);
			optTrees << "tree p" << prec1 << ".r" << rep << " = [&U] " << treeString << ";\n";
			//restore the randomization
			indiv[0].CopySecByRearrangingNodesOfFirst(indiv[0].treeStruct, &indiv[1], true);
			indiv[0].SetDirty();
			indiv[0].CalcFitness(0);
	//		scoresThisPrec.push_back(indiv[0].Fitness());
//			passesThisPrec.push_back(pass);
//			derivCalcsThisPrec.push_back(optCalcs);
//			for(int b=1;b<numNodes;b++) blensThisRep.push_back(indiv[0].treeStruct->allNodes[b]->dlen);
//			blensThisPrec.push_back(blensThisRep);
//			blensThisRep.clear();
			optCalcs = 0;
			//indiv[0].CopySecByRearrangingNodesOfFirst(indiv[0].treeStruct, &tempIndiv, true);
			}
//		finalScores.push_back(scoresThisPrec);
//		numPasses.push_back(passesThisPrec);
//		numDerivCalcs.push_back(derivCalcsThisPrec);
//		allBlens.push_back(blensThisPrec);
//		AppendTreeToTreeLog(0, 0);
//		scoresThisPrec.clear();
//		passesThisPrec.clear();
//		derivCalcsThisPrec.clear();
//		blensThisPrec.clear();
//		out << prec[p] << "\t";
		}
//	out << "\n";
	for(int precNum = 0;precNum < finalScores.size();precNum++){
		for(int rep = 0;rep < finalScores[precNum].size();rep++){
//			for(vector<doubvec>::iterator it = finalScores.begin();it != scores.end();it++){
			
			out << prec[precNum] << "\t" << rep << "\t" << finalScores[precNum][rep] << "\t" << numPasses[precNum][rep] << "\t" << numDerivCalcs[precNum][rep] << endl;

/*			for(vector<doubvec>::iterator it = scores.begin();it != scores.end();it++){
				out << (*it)[rep] << "\t";
				}
			for(vector<intvec>::iterator it = passes.begin();it != passes.end();it++){
				out << (*it)[rep] << "\t";
				}
			out << "\n";
*/			}
		}


	ofstream blens;
	for(int precNum = 0;precNum < finalScores.size();precNum++){
		char filename[100];
		if(reducing)
			sprintf(filename, "blens.%s.final.log", conf->ofprefix.c_str());
		else
			sprintf(filename, "blens.%s.%f.log", conf->ofprefix.c_str(), prec[precNum]);
		blens.open(filename);
		blens << "branch#\tfullyOpt\treps...\n";
		//careful here - the number of nodes includes the root, which has no blen and wasn't put into the 
		//blen vector. So, the indexing is [actualNodeNum - 1]
		for(int bnum=0;bnum<numNodes - 1;bnum++){
			blens << bnum+1 << "\t";
			//toss in the blens for one of the reps for the final prec, which we assume will be fully optimal
			blens << allBlens[allBlens.size()-1][0][bnum] << "\t";
			for(int rep = 0;rep < finalScores[precNum].size();rep++){
				blens << allBlens[precNum][rep][bnum] << "\t";
				}
			blens << endl;
			}
		blens.close();
		}

	randTrees << "end;\n";
	optTrees << "end;\n";
	randTrees.close();
	optTrees.close();

	FinalizeOutputStreams(0);
	FinalizeOutputStreams(1);
	FinalizeOutputStreams(2);
	}

void Population::QuickSort( FLOAT_TYPE **scoreArray, int top, int bottom ){

	int i = top;
	int j = bottom;
	FLOAT_TYPE x = scoreArray[ (top + bottom) / 2 ][1];
	do {
		while( scoreArray[i][1] < x  &&  i < bottom ) i++ ;
		while( x < scoreArray[j][1]  &&  j > top  ) j-- ;

		if( i <= j ) {
			for( int k = 0; k < 2; k++ ) {
				FLOAT_TYPE y = scoreArray[i][k];
				scoreArray[i][k] = scoreArray[j][k];
				scoreArray[j][k] = y;
			}
			i++;
			if(j) j--;
		}

	} while( i <= j );

	if( top  <    j    ) QuickSort( scoreArray, top, j );
	if(  i   <  bottom ) QuickSort( scoreArray, i, bottom );
}

FLOAT_TYPE Population::CalcAverageFitness(){
	FLOAT_TYPE total = ZERO_POINT_ZERO;
	
	for(unsigned i = 0; i < total_size; i++ ){
		// evaluate fitness
		if(indiv[i].IsDirty()){
			indiv[i].CalcFitness(subtreeNode);
			}
		assert(indiv[i].Fitness() != 1);
	
		total += indiv[i].Fitness();
		cumfit[i][0] = (FLOAT_TYPE)i;
		cumfit[i][1] = indiv[i].Fitness();
		}

	FLOAT_TYPE avg = total / (FLOAT_TYPE)total_size;

	// Sort fitnesses from low to high (bad to good)
	QuickSort( cumfit, 0, total_size-1 );

	// keep track of which individual is most fit each generation we've stored the 
	//fitnesses as ln-likelihoods in cumfit, so cumfit[0] will be the _least_ fit individual
	int mostFit = total_size-1;
#ifndef NO_EVOLUTION
	bestAccurateIndiv=bestIndiv = (int)cumfit[mostFit][0];
#else
	bestAccurateIndiv=bestIndiv = 0;
#endif
	//if subtree mode is active, we also want to find the best accurate indiv
	if(rank==0){
		while(paraMan->subtreeModeActive==true && indiv[bestAccurateIndiv].accurateSubtrees==false){
			mostFit--;
			bestAccurateIndiv=(int)cumfit[mostFit][0];
			}
		assert(mostFit>=0);
		}
	
	// keep track of all-time best
	if( indiv[bestIndiv].Fitness() > prevBestFitness ){
		prevBestFitness = bestFitness;
		globalBest = bestFitness = indiv[bestIndiv].Fitness();
		}

	if(memLevel>0){
		//if we are at some level of memory restriction, mark the clas of the old best
		//for reclamation, and protect those of the new best
		SetNewBestIndiv(bestIndiv);
		}

	CalculateReproductionProbabilies(cumfit, conf->selectionIntensity, total_size);
	return avg;
	
/*	Here's Paul's original selection criterion, based solely on rank
	//	
	// relative fitnesses are assigned based solely on position
	// of individual in sorted array - we forget the likelihoods (or treelengths)
	// at this point.  This allows the likelihoods to be close together
	// and still get a healthy distribution of relative fitnesses so that
	// there is real differential reproduction

	FLOAT_TYPE n = (FLOAT_TYPE)total_size;
	FLOAT_TYPE nn = n * ( n + 1.0 );
	FLOAT_TYPE incr = 2.0 / nn;
	FLOAT_TYPE cum = incr;

	cumfit[0][1] = cum;
	for( i = 1; i < total_size; i++ ) {
		cum += incr;
		cumfit[i][1] = cumfit[i-1][1] + cum;
	}
*/
}

void Population::CalculateReproductionProbabilies(FLOAT_TYPE **scoreArray, FLOAT_TYPE selectionIntensity, int indivsInArray){
	//DJZ 2-28-06 Generalizing this so that it can be used in multiple places with different 
	//subsets of individuals and selection intensities.  The 2-d array passed in (indivsInArray x 2)
	//has the scores in the [x][1] slots, and the indiv numbers in the [x][0] slots, and should already
	//be sorted from low to high (bad to good). The reproduction probs will be placed in the [x][1] before returning.

	//Probability of reproduction based on more or less on AIC weights, although
	//the strength of selection can be varied by changing the selectionIntensity
	//A selectionIntensity of 0.5 makes this equivalent to AIC weights, while 
	//smaller number makes the selection less severe
	FLOAT_TYPE *deltaAIC=new FLOAT_TYPE[indivsInArray];
	FLOAT_TYPE tot=ZERO_POINT_ZERO;

	for(int i=0;i<indivsInArray-1;i++){
		deltaAIC[i]=scoreArray[indivsInArray-1][1] - scoreArray[i][1];
		deltaAIC[i]=exp(-selectionIntensity * deltaAIC[i]);
		tot+=deltaAIC[i];
		}

	if(indivsInArray == total_size){
		deltaAIC[indivsInArray-1]=conf->holdoverPenalty;
		}
	else deltaAIC[indivsInArray-1]=ZERO_POINT_ZERO;

	deltaAIC[indivsInArray-1]=exp(-selectionIntensity * deltaAIC[indivsInArray-1]);
	tot+=deltaAIC[indivsInArray-1];


	for(int i=0;i<indivsInArray;i++)
		deltaAIC[i] /= tot;
	
	FLOAT_TYPE cum=deltaAIC[0];
	scoreArray[0][1] = cum;
	for(int i = 1; i < indivsInArray; i++ ) {
		cum += deltaAIC[i];
		scoreArray[i][1] = cum;
		}
	delete []deltaAIC;
	assert(abs(scoreArray[indivsInArray-1][1] - ONE_POINT_ZERO) < 0.001);
}
/*
void Population::CreateGnuPlotFile()
{
	char tmpstr[51];

	ofstream gnuf( params->gnufname );
	assert( gnuf );

	// set labels
	gnuf << "set xlabel \"Generation\"" << endl;
	gnuf << "set ylabel \"Fitness\"" << endl;
	gnuf << "set title \"" << params->plottitle << "\"" << endl;

	// alternate title containing config settings
	gnuf << "#set title \"";
	gnuf << "N=" << conf->nindivs;
	gnuf << " h=" << params->holdover;
	gnuf << " b=" << params->meanBrlenMuts;
	gnuf << " s=" << params->gammaShapeBrlen;
	gnuf << " seed=" << rnd.init_seed();
	gnuf << '\"' << endl;

	// make alternate output to mif file available
	gnuf << "#set terminal mif" << endl;
	gnuf << "#set output \"" << params->ofprefix << ".mif\"" << endl;

	// place legend on graph
	strcpy( tmpstr, "set key 5000,-45000.0" );
	gnuf << tmpstr << endl;

	// finally, the plot command
	gnuf << "plot \"" << params->logfname;
	if( params->fatlog ) {
		gnuf << "\" using 1:3 title \"best\" with lines";
		gnuf << ", \"" << params->logfname;
		gnuf << "\" using 1:2 title \"average\" with lines" << endl;
   }
   else {
		gnuf << "\" using 1:2 with lines" << endl;
   }

	gnuf << "pause -1" << endl;

	gnuf.close();
}
*/
void Population::DetermineParentage(){
	//determine each individual's parentage
	unsigned parent;
	FLOAT_TYPE r;

	for(unsigned i = 0; i < conf->nindivs; i++ ){

#ifndef NO_EVOLUTION
		if( i < conf->holdover ){// copy best individual's genotype to next generation
			if(rank==0){
				if(paraMan->subtreeModeActive==true) parent=bestAccurateIndiv;
				else parent=bestIndiv;
				}
			else parent=bestIndiv;
			}
		else if(rank==0)
			if(i==1 && rank==0 && indiv[bestIndiv].accurateSubtrees==true && paraMan->ReadyForSubtreeRecom(gen)){
			//if subtree mode is on and we haven't tried a subtreeRecom in a while, set up an individual for that
				parent=bestIndiv;
				newindiv[i].mutation_type=Individual::subtreeRecom;
				}
			else {// find a parent
				r = rnd.uniform();
				for( parent = 0; parent < total_size; parent++ ){
					if( r < cumfit[parent][1] ) break;
					}
				parent = (int)cumfit[parent][0];

#ifdef INPUT_RECOMBINATION
			
			paraMan->maxRecomIndivs = 3;
			paraMan->nremotes = NUM_INPUT;

			if(rank==0 && paraMan->subtreeModeActive==false && i>= (conf->nindivs - paraMan->maxRecomIndivs)){
/*				int *mates=new int[paraMan->nremotes];
				for(int j=0;j<paraMan->nremotes;j++) mates[j]=conf->nindivs+j;
				ScrambleArray(paraMan->nremotes, mates);
*/
				int foo=2;
				FLOAT_TYPE **recomSelect=new FLOAT_TYPE *[paraMan->nremotes];
				for(int q=0;q<paraMan->nremotes;q++)
					recomSelect[q]=new FLOAT_TYPE[2];
					
				int potentialPartners=0;
				for(int r=0;r<paraMan->nremotes;r++){
					int ind=conf->nindivs+r;
					recomSelect[r][0]=(FLOAT_TYPE)(ind);
					if(ind==parent //don't recombine with your parent
						|| (indiv[parent].topo == indiv[ind].topo) //don't recombine with another of the same topo	
						|| (indiv[ind].willrecombine == true))//don't recombine with someone who is already doing so		
						recomSelect[r][1]=-1e100;	
					else{
						recomSelect[r][1]=indiv[ind].Fitness();
						potentialPartners++;
						}
					}
				if(potentialPartners > 0){
					QuickSort(recomSelect, 0, paraMan->nremotes-1);
					CalculateReproductionProbabilies(recomSelect, 0.001, paraMan->nremotes);			
									
					int mateIndex;
					int curMate;
					// find someone else to recombine with

					FLOAT_TYPE r=rnd.uniform();
					for( mateIndex=0;mateIndex < paraMan->nremotes;mateIndex++)
						if( r < recomSelect[mateIndex][1]) break;
					curMate=recomSelect[mateIndex][0];


					newindiv[i].recombinewith=curMate;
					indiv[curMate].willrecombine=true;
					//this will be a new topology, so mark it as topo -1.  This will be dealt with when we update the topolist
					newindiv[i].topo=-1;
					}
				for(int q=0;q<paraMan->nremotes;q++)
					delete []recomSelect[q];
				delete []recomSelect;
				}
#endif

#ifdef MPI_VERSION
//new bipart recom conditions, 9-25-05

//DJZ 2-28-06 making recombination partner weakly tied to fitness (selctionIntensity of 0.01) rather than random

			if(rank==0 && paraMan->subtreeModeActive==false && i>= (conf->nindivs - paraMan->maxRecomIndivs)){
/*				int *mates=new int[paraMan->nremotes];
				for(int j=0;j<paraMan->nremotes;j++) mates[j]=conf->nindivs+j;
				ScrambleArray(paraMan->nremotes, mates);
*/
				int foo=2;
				FLOAT_TYPE **recomSelect=new FLOAT_TYPE *[paraMan->nremotes];
				for(int q=0;q<paraMan->nremotes;q++)
					recomSelect[q]=new FLOAT_TYPE[2];
					
				int potentialPartners=0;
				for(int r=0;r<paraMan->nremotes;r++){
					int ind=conf->nindivs+r;
					recomSelect[r][0]=(FLOAT_TYPE)(ind);
					if(ind==parent //don't recombine with your parent
						|| (indiv[parent].topo == indiv[ind].topo) //don't recombine with another of the same topo	
						|| (indiv[ind].willrecombine == true))//don't recombine with someone who is already doing so		
						recomSelect[r][1]=-1e100;	
					else{
						recomSelect[r][1]=indiv[ind].Fitness();
						potentialPartners++;
						}
					}
				if(potentialPartners > 0){
					QuickSort(recomSelect, 0, paraMan->nremotes-1);
					CalculateReproductionProbabilies(recomSelect, 0.01, paraMan->nremotes);			
									
					int mateIndex;
				int curMate;
				// find someone else to recombine with

					FLOAT_TYPE r=rnd.uniform();
					for( mateIndex=0;mateIndex < paraMan->nremotes;mateIndex++)
						if( r < recomSelect[mateIndex][1]) break;
					curMate=recomSelect[mateIndex][0];


					newindiv[i].recombinewith=curMate;
					indiv[curMate].willrecombine=true;
					//this will be a new topology, so mark it as topo -1.  This will be dealt with when we update the topolist
					newindiv[i].topo=-1;
					}
				for(int q=0;q<paraMan->nremotes;q++)
					delete []recomSelect[q];
				delete []recomSelect;
				}
#endif
		}
#else //ifdef NO_EVOLUTION
		parent = 0;
#endif
		
		newindiv[i].parent=parent;
		if(newindiv[i].mutation_type==Individual::subtreeRecom) newindiv[i].topo=-1; //VERIFY
		else newindiv[i].topo=indiv[parent].topo;
		indiv[ parent ].willreproduce=true;
		}
	}

void Population::FindTreeStructsForNextGeneration(){
	//find treestructs for all of the newindivs, either by getting an unused one from the previous
	//generation or by getting one from the unusedTree stack
	for(unsigned i = 0; i < total_size; i++ ){
		//see if the parent indiv has already been used in the new generation, or if it will recombine
		if( i < conf->nindivs && (indiv[newindiv[i].parent].reproduced||indiv[newindiv[i].parent].willrecombine )){	      
			//use a tree from the unused Indiv stack.  If it is empty, create an extra indiv that will 
			//eventually make it's way /back to that stack.  At most we should only ever have nindiv 
			//trees in the unused stack
			Tree *destPtr;
			if(unusedTrees.empty()){//create a new tree
				Tree *ttree=new Tree();
				destPtr=ttree;
				}
			else{
				destPtr=*(unusedTrees.end()-1);
				unusedTrees.pop_back();
				}
			newindiv[i].CopySecByRearrangingNodesOfFirst(destPtr,&indiv[newindiv[i].parent]);
			}
		else{
			//if the tree will not be used in recombination and has not already been used
			newindiv[i].CopyByStealingTree(&indiv[newindiv[i].parent]);
			indiv[ newindiv[i].parent].reproduced=true;
			if(i>conf->nindivs) newindiv[i].mutation_type=indiv[i].mutation_type;
			}
		}
	}
	
void Population::PerformMutation(int indNum){
	Individual *ind=&newindiv[indNum];
	Individual *par=&indiv[newindiv[indNum].parent];

	//FLOAT_TYPE beforeScore;
	bool recomPerformed;

	switch(ind->mutation_type){
/*		case Individual::exNNI: //exNNI and exlimSPR trump all other mutation types
			beforeScore=par->Fitness();
			NNIoptimization(indNum, 1);
			if(beforeScore==ind->Fitness()){
				topologies[ind->topo]->exNNItried=true;
				}
			//ind->accurateSubtrees=false;
			break;
		
		case Individual::exlimSPR:
			assert(0);
			SPRoptimization(indNum);
			ind->accurateSubtrees=false;
			break;
*/		

		case Individual::subtreeRecom:
			//perform subtree recom, which melds together the different subtrees worked on by the
			//remote nodes
			recomPerformed=SubtreeRecombination(indNum);
			if(recomPerformed==false) ind->mutation_type=0;
//			ind->treeStruct->calcs=calcCount;
//			calcCount=0;
			ind->CalcFitness(0);
			break;
			
		default:
			if(ind->recombinewith>-1){// perform recombination
				Individual *recompar=&indiv[ind->recombinewith];
				//don't want to standardize biparts anymore
				ind->treeStruct->CalcBipartitions(false);
				recompar->treeStruct->CalcBipartitions(false);
				ind->CrossOverWith( *recompar, adap->branchOptPrecision);
				ind->accurateSubtrees=false;
//				ind->treeStruct->calcs=calcCount;
//				calcCount=0;
				}
			if(ind->recombinewith==-1){//all types of "normal" mutation that occur at the inidividual level
				if(rank==0){//if we are the master
				 	if(ind->accurateSubtrees==false || paraMan->subtreeModeActive==false){

			       		ind->Mutate(adap->branchOptPrecision, adap);

						if(output_tree){
							treeLog << "  tree gen" << gen <<  "." << indNum << "= [&U] [" << ind->Fitness() << "][ ";
							string modstr;
							ind->modPart.FillGarliFormattedModelStrings(modstr);
							ind->treeStruct->root->MakeNewick(treeString, false, true);
							treeLog << modstr.c_str() << "]" << treeString << ";" << endl;
							output_tree=false;
							}

			       		//reclaim clas if the created tree has essentially no chance of reproducing
			       		if(((ind->Fitness() - BestFitness()) < (-11.5/conf->selectionIntensity))){
			       			ind->treeStruct->ReclaimUniqueClas();
			       			}
						}
					else{
						assert(0);//7/21/06 subtree mode would need to be updated to work again

						//if subtree mode is on and we are the master, mutate one of the nodes
						//that isn't in a subtree, or alternatively pick a subtree and mutate it
/*					#ifndef MASTER_DOES_SUBTREE
						if(paraMan->fewNonSubtreeNodes != true)
							ind->NonSubtreeMutate(paraMan, adap->branchOptPrecision, adap);
						else 
							ind->SubtreeMutate(subtreeNode, adap->branchOptPrecision, subtreeMemberNodes, adap);
					#else
						ind->SubtreeMutate(subtreeNode, adap->branchOptPrecision, subtreeMemberNodes, adap);
					#endif					
*/						}
					}
				else{//if we are a remote node
				 	if(subtreeNode==0) ind->Mutate(adap->branchOptPrecision, adap);
					else{
						assert(0);
						//ind->SubtreeMutate(subtreeNode, adap->branchOptPrecision, subtreeMemberNodes, adap);
						}
					}
				}
			}
		
		//check the accuracy of the subtrees
		#ifndef NDEBUG
		if(rank==0 && ind->accurateSubtrees==true)
			paraMan->CheckSubtreeAccuracy(ind->treeStruct);
		#endif
	}

void Population::NextGeneration(){

	DetermineParentage();

	FindTreeStructsForNextGeneration();

	//return any treestructs from the indivs that won't be used in recombination
	//and weren't used to make the newindivs.  This is necessary to keep from having
	//too many CLAs in use at any one time 
	for(unsigned j=0;j<conf->nindivs;j++){
		if(indiv[j].reproduced==false && indiv[j].willrecombine==false){
			//this reclaims all indiv's treestructs who have no offspring and no recombination partner
			indiv[j].treeStruct->RemoveTreeFromAllClas();
			unusedTrees.push_back(indiv[j].treeStruct);
			indiv[j].treeStruct=NULL;
			}
		}

	//to simplify all of the scoring that will be coming up (without passing 
	//a bunch of crap), set the models of the trees to correspond to that of the individuals
	UpdateTreeModels();
	
	//this loop is only for mutation and recom, so start from holdover
	for(unsigned indnum = conf->holdover; indnum < conf->nindivs; indnum++ ){
		PerformMutation(indnum);
		}

	UpdateTreeModels();

	//the only trees that we need to return at this point are ones that
	//did not reproduce AND were used in recom.  Those that weren't used 
	//in recom were already reclaimed above, and the treestructs set to NULL
	for(unsigned j=0;j<conf->nindivs;j++){
		if(indiv[j].reproduced==false && indiv[j].treeStruct!=NULL){
			indiv[j].treeStruct->RemoveTreeFromAllClas();
			unusedTrees.push_back(indiv[j].treeStruct);
			}
		//reset all of the individuals
		indiv[j].ResetIndiv();
		}

	// swap newindiv and indiv
	for(unsigned i=0;i<conf->nindivs;i++)
		indiv[i].ResetIndiv();
	for(unsigned i=conf->nindivs;i<total_size;i++)
		indiv[i].willrecombine=false;

	Individual* tmp = newindiv;
	newindiv = indiv;
	indiv = tmp;

	CalcAverageFitness(); //score individuals that need it
		
	#ifdef DEBUG_SCORES
	if(rank==0)	OutputFilesForScoreDebugging();
	#endif
	
	}

void Population::OutputFate(){
	//output everything that happened to each indiv in this generation to file

	for(unsigned i=0;i<total_size;i++){
		fate << 	gen << "\t" << i << "\t" << indiv[i].parent << "\t";

#ifdef MPI_VERSION
		fate << indiv[i].recombinewith << "\t";	
#endif
		fate << indiv[i].Fitness() << "\t" << indiv[i].mutation_type << "\t" << indiv[i].mutated_brlen << "\t";
#ifdef MPI_VERSION
	    fate  << indiv[i].accurateSubtrees << "\t";
#endif
		
	    fate << stopwatch.SplitTime() << "\t" << adap->branchOptPrecision;

//some extra debugging info
/*		fate << "\t" << indiv[i].topo << "\t";	
	    fate << indiv[i].treeStruct->calcs << "\t";
	    indiv[i].treeStruct->calcs=0;
	    int c, tr, r;
	    indiv[i].treeStruct->CountNumReservedClas(c, tr, r);
	    fate << c << "\t" << tr << "\t" << r << "\t";
//	    
*/	    fate << "\n";
		}
//	fate << claMan->NumFreeClas() << "\n";
	if(gen%20 ==0) fate.flush();
	}

void Population::OutputFilesForScoreDebugging(Individual *ind /*=NULL*/, int num){
	//create three files, one with all of the trees in each gen in nexus format
	//one with a paup block specifiying the scoring of the trees, and one containing 
	//a list of the scores from GAML

if(rank > 0) return;	

//ofstream outf;
//ofstream paupf;

#ifdef NNI_SPECTRUM

	char fname1[30];
	char fname2[30];
	sprintf(fname1, "toscore%d.tre", gen);
	sprintf(fname2, "toscore%d.nex", gen);
	if(num==1){
		outf.open(fname1);
		paupf.open(fname2);
		}
	else{
		outf.open(fname1, ios::app);
		paupf.open(fname2, ios::app);
		}
#endif

    //not sure what the original intended order of eval was here, but it caused warnings
    //if(gen==1 && ind==NULL || num==1){
	if(gen==1 && (ind==NULL || num==1)){
		dataPart->BeginNexusTreesBlock(outf);
		
		paupf << "#nexus\n\n";
		paupf << "begin paup;\n";
		paupf << "set warnreset=no incr=auto;\n";
		paupf << "execute " << conf->ofprefix.c_str() << ".nex;\n";
#ifndef NNI_SPECTRUM
		paupf << "gett file=toscore.tre storebr;" << endl;
#else
		paupf << "gett file=" << outf << " storebr;" << endl;
#endif
		}
		
	if(ind==NULL){
		for(unsigned i=0;i<total_size;i++){
			outf << "  utree " << gen << i << "= ";
			indiv[i].treeStruct->root->MakeNewick(treeString, false, true, true);
			outf << treeString << ";\n";
			
			paupf << "lset userbr ";
			//DEBUG PARTITION
/*			if(modSpec.Nst()==2) paupf << "nst=2 trat=" << indiv[i].mod->Rates(0) << " base=(" << indiv[i].mod->StateFreq(0) << " " << indiv[i].mod->StateFreq(1) << " " << indiv[i].mod->StateFreq(2) << ");\n" << "lsc " << (gen-1)*conf->nindivs+i+1;

			else paupf << "nst=6 rmat=(" << indiv[i].mod->Rates(0) << " " << indiv[i].mod->Rates(1) << " " << indiv[i].mod->Rates(2) << " " << indiv[i].mod->Rates(3) << " " << indiv[i].mod->Rates(4) << ") " << " base=(" << indiv[i].mod->StateFreq(0) << " " << indiv[i].mod->StateFreq(1) << " " << indiv[i].mod->StateFreq(2) << ") ";
			
#ifdef FLEX_RATES			
			paupf << "[FLEX RATES] ";
#else
			if(indiv[i].mod->NRateCats()>1) paupf << "rates=gamma shape=" << indiv[i].mod->Alpha() << " ";
			paupf << "pinv=" << indiv[i].mod->PropInvar() << " "; 
#endif

*/
			if(modSpecSet.GetModSpec(0)->Nst()==2) paupf << "nst=2 trat=" << indiv[i].modPart.GetModel(0)->Rates(0) << " base=(" << indiv[i].modPart.GetModel(0)->StateFreq(0) << " " << indiv[i].modPart.GetModel(0)->StateFreq(1) << " " << indiv[i].modPart.GetModel(0)->StateFreq(2) << ");\n" << "lsc " << (gen-1)*conf->nindivs+i+1;

			else paupf << "nst=6 rmat=(" << indiv[i].modPart.GetModel(0)->Rates(0) << " " << indiv[i].modPart.GetModel(0)->Rates(1) << " " << indiv[i].modPart.GetModel(0)->Rates(2) << " " << indiv[i].modPart.GetModel(0)->Rates(3) << " " << indiv[i].modPart.GetModel(0)->Rates(4) << ") " << " base=(" << indiv[i].modPart.GetModel(0)->StateFreq(0) << " " << indiv[i].modPart.GetModel(0)->StateFreq(1) << " " << indiv[i].modPart.GetModel(0)->StateFreq(2) << ") ";
			
#ifdef FLEX_RATES			
			paupf << "[FLEX RATES] ";
#else
			if(indiv[i].modPart.GetModel(0)->NRateCats()>1) paupf << "rates=gamma shape=" << indiv[i].modPart.GetModel(0)->Alpha() << " ";
			paupf << "pinv=" << indiv[i].modPart.GetModel(0)->PropInvar() << " "; 
#endif

			if(gen==1 && i==0) paupf << ";\n" << "lsc " << (gen-1)*total_size+i+1 << "/scorefile=paupscores.txt replace;\n";
			else paupf << ";\n" << "lsc " << (gen-1)*total_size+i+1 << "/scorefile=paupscores.txt append;\n";
			}
		}
	else{
		outf << "  utree " << num << "= ";
		ind->treeStruct->root->MakeNewick(treeString, false, true);
		outf << treeString << ";\n";

		//DEBUG PARTITION
/*		paupf << "lset userbr ";
		if(modSpec->Nst()==2) paupf << "nst=2 trat=" << ind->mod->Rates(0) << " base=(" << ind->mod->StateFreq(0) << " " << ind->mod->StateFreq(1) << " " << ind->mod->StateFreq(2) << ");\nlsc ";
		
		else paupf << "nst=6 rmat=(" << ind->mod->Rates(0) << " " << ind->mod->Rates(1) << " " << ind->mod->Rates(2) << " " << ind->mod->Rates(3) << " " << ind->mod->Rates(4) << ") " << " base=(" << ind->mod->StateFreq(0) << " " << ind->mod->StateFreq(1) << " " << ind->mod->StateFreq(2) << ") ";

#ifdef FLEX_RATES			
			paupf << "[FLEX RATES] ";
#else	
		if(ind->mod->NRateCats()>1) paupf << "rates=gamma shape=" << ind->mod->Alpha() << " ";
		paupf << "pinv=" << ind->mod->PropInvar() << " "; 
#endif
*/
		paupf << "lset userbr ";
		if(modSpecSet.GetModSpec(0)->Nst()==2) paupf << "nst=2 trat=" << ind->modPart.GetModel(0)->Rates(0) << " base=(" << ind->modPart.GetModel(0)->StateFreq(0) << " " << ind->modPart.GetModel(0)->StateFreq(1) << " " << ind->modPart.GetModel(0)->StateFreq(2) << ");\nlsc ";
		
		else paupf << "nst=6 rmat=(" << ind->modPart.GetModel(0)->Rates(0) << " " << ind->modPart.GetModel(0)->Rates(1) << " " << ind->modPart.GetModel(0)->Rates(2) << " " << ind->modPart.GetModel(0)->Rates(3) << " " << ind->modPart.GetModel(0)->Rates(4) << ") " << " base=(" << ind->modPart.GetModel(0)->StateFreq(0) << " " << ind->modPart.GetModel(0)->StateFreq(1) << " " << ind->modPart.GetModel(0)->StateFreq(2) << ") ";

#ifdef FLEX_RATES			
			paupf << "[FLEX RATES] ";
#else	
		if(ind->modPart.GetModel(0)->NRateCats()>1) paupf << "rates=gamma shape=" << ind->modPart.GetModel(0)->Alpha() << " ";
		paupf << "pinv=" << ind->modPart.GetModel(0)->PropInvar() << " "; 
#endif

#ifndef NNI_SPECTRUM
		if(num==1) paupf << ";\n" << "lsc " << num << "/scorefile=paupscores.txt replace;\n";
		else paupf << ";\n" << "lsc " << num << "/scorefile=paupscores.txt append;\n";
#else
		if(num==1) paupf << ";\n" << "lsc " << num << "/scorefile=paupscores" << gen << ".txt replace;\n";
		else paupf << ";\n" << "lsc " << num << "/scorefile=paupscores"  << gen << ".txt append;\n";
#endif
		}
#ifdef NNI_SPECTRUM
outf.close();
paupf.close();		
#endif
	}	

//this assumes that the tree to be appended is a member of the population
//if indNum is -1, then the bestIndiv from the pop is used
void Population::AppendTreeToTreeLog(int mutType, int indNum /*=-1*/){
	
	if(treeLog.is_open() == false || conf->outputTreelog==false) 
		return;

	const Individual *ind;
	int i = (indNum >= 0 ? indNum : bestIndiv);

	ind=&indiv[i];

//	if(Tree::outgroup != NULL) 
//		OutgroupRoot(ind, i);
		
	int num = 0;
	Individual tempInd;
	const Individual *theInd;
	if(Tree::outgroup != NULL || conf->collapseBranches){
		tempInd.DuplicateIndivWithoutCLAs(ind);
		tempInd.treeStruct->noCalcs = true;
		if(Tree::outgroup != NULL)
			OutgroupRoot(&tempInd, -1);
/*		//Can't decide if these should be collapsed or not here.  Thinking no.
		if(conf->collapseBranches){
			tempInd.treeStruct->root->CollapseMinLengthBranches(num);
			outman.UserMessage("%d COLLAPSED", num);
			}
*/
		theInd = &tempInd;
		}
	else
		theInd = ind;

	if(finishedRep)
		treeLog << "  tree final= [&U] [" << theInd->Fitness() << "][ ";
	else 
		treeLog << "  tree gen" << gen <<  "= [&U] [" << theInd->Fitness() << "\tmut=" << mutType << "][ ";
	string modstr;
	ind->modPart.FillGarliFormattedModelStrings(modstr);
	theInd->treeStruct->root->MakeNewick(treeString, false, true);
	treeLog << modstr.c_str() << "]" << treeString << ";" << endl;
	}


void Population::FinishBootstrapRep(const Individual *ind, int rep){

	if(bootLog.is_open() == false) 
		return;

	int num = 0;
	Individual tempInd;
	const Individual *theInd;
	if(Tree::outgroup != NULL || conf->collapseBranches){
		tempInd.DuplicateIndivWithoutCLAs(ind);
		if(Tree::outgroup != NULL)
			OutgroupRoot(&tempInd, -1);
		if(conf->collapseBranches){
			tempInd.treeStruct->root->CollapseMinLengthBranches(num);
//			outman.UserMessage("%d COLLAPSED", num);
			}
		theInd = &tempInd;
		}
	else
		theInd = ind;

	bootLog << "  tree bootrep" << rep <<  "= [&U] [" << theInd->Fitness() << " ";
	
	string modstr;
	theInd->modPart.FillGarliFormattedModelStrings(modstr);

	theInd->treeStruct->root->MakeNewick(treeString, false, true);
	bootLog << modstr.c_str() << "] " << treeString << ";" << endl;

	if(conf->outputPhylipTree)
		WritePhylipTree(bootLogPhylip);
	}

bool Population::OutgroupRoot(Individual *ind, int indnum){
	//if indnum != -1 the individual is in the indiv array, and a few extra things need to be done 
	ind->treeStruct->root->CheckforPolytomies();

	ind->treeStruct->CalcBipartitions(true);
	Bipartition b = *(Tree::outgroup);
	b.Standardize();
	TreeNode *r = ind->treeStruct->ContainsBipartitionOrComplement(b);
	
	if(r == NULL){
		//this means that there isn't a bipartition separating the outgroup and ingroup
		//so outgroup rooting is not possible
		return false;
		}
	
	TreeNode *temp = r;
	while(temp->IsTerminal() == false) 
		temp=temp->left;
	if(Tree::outgroup->ContainsTaxon(temp->nodeNum) == false || r->IsTerminal()) 
		r = r->anc;
	if(r->IsNotRoot()){
//		outman.UserMessage("REROOTED");
		ind->treeStruct->RerootHere(r->nodeNum);
		if(indnum != -1){
			ind->SetDirty();
			ind->CalcFitness(0);
			}
		return true;
		}
	else return false;
	}

void Population::WriteTreeFile( const char* treefname, int indnum, bool collapse /*=false*/ ){
	assert( treefname );
	string filename = treefname;
	filename += ".tre";

	//output an individual from the storedTrees if an indnum is passed in
	//otherwise the best in the population
	const Individual *ind;
	if(indnum == -1){
		ind = &indiv[bestIndiv];
		}
	else{
		assert(indnum < storedTrees.size());
		ind = storedTrees[indnum];
		}

	int num = 0;
	Individual tempInd;
	const Individual *theInd;
	if(Tree::outgroup != NULL || (conf->collapseBranches && collapse)){
		tempInd.DuplicateIndivWithoutCLAs(ind);
		if(Tree::outgroup != NULL)
			OutgroupRoot(&tempInd, -1);
		if(conf->collapseBranches && collapse){
			tempInd.treeStruct->root->CollapseMinLengthBranches(num);
//			outman.UserMessage("%d COLLAPSED", num);
			}
		theInd = &tempInd;
		}
	else
		theInd = ind;

#ifdef INCLUDE_PERTURBATION
	if(allTimeBest != NULL){
		if(best->Fitness() < allTimeBest->Fitness() || pertMan->ratcheted==true) return;
		}
#endif

#ifdef BOINC
	char physical_name[100];
	boinc_resolve_filename(filename.c_str(), physical_name, sizeof(physical_name));
	MFILE outf;
	outf.open(physical_name, "w");
#else
	ofstream outf;
	outf.open( filename.c_str() );
	outf.precision(8);
#endif
	string trans;
	string str;
	dataPart->BeginNexusTreesBlock(trans);
	//data->BeginNexusTreesBlock(outf);
	char temp[101];	
	if(indnum == -1)
		sprintf(temp, "tree best = [&U][!GarliScore %f][!GarliModel ", theInd->Fitness());
	else
		sprintf(temp, "tree bestREP%d = [&U][!GarliScore %f][!GarliModel ", indnum+1, theInd->Fitness());
	str += temp;
	string modstr;
	ind->modPart.FillGarliFormattedModelStrings(modstr);
	str += modstr;
	str += "]";

#ifdef BOINC
	const char *s = trans.c_str();
	outf.write(s, sizeof(char), trans.length());
	s = str.c_str();
	outf.write(s, sizeof(char), str.length());
	theInd->treeStruct->root->MakeNewick(treeString, false, true);
	size_t len = strlen(treeString);
	outf.write(treeString, sizeof(char), len);
	str = ";\nend;\n";
	s = str.c_str();
	outf.write(s, sizeof(char), str.length());
#else
	outf << trans;
	outf << str;
	outf.setf( ios::floatfield, ios::fixed );
	outf.setf( ios::showpoint );
	theInd->treeStruct->root->MakeNewick(treeString, false, true);
	outf << treeString << ";\n";
	outf << "end;\n";
#endif	
	//add a paup block setting the model params
	str = "";
	if(modSpecSet.GetModSpec(0)->IsNucleotide()){
		if(ind->modPart.NumModels() == 1){
			ind->modPart.GetModel(0)->FillPaupBlockStringForModel(str, filename.c_str());
			}
		else{
			str += "[\n";
			//modstr was already filled above
			//ind->modPart.FillGarliFormattedModelStrings(modstr);
			str += modstr;
			str += "\n]\n";
			modstr.clear();
			}
		}
#ifdef BOINC
	s = str.c_str();
	outf.write(s, sizeof(char), str.length());
	if((userTermination || timeTermination) && (indnum == storedTrees.size() - 1)){
		//str = "[!****NOTE: GARLI Run was terminated before termination condition was reached!\nLikelihood scores, topologies and model estimates obtained may not be fully optimal!****\n]";
		str = "[";
		str += TerminationWarningMessage();
		str += "]\n";
		s = str.c_str();
		outf.write(s, sizeof(char), str.length());
		}
#else
	outf << str; 
	//if(indnum < 0 && (userTermination || timeTermination))
	if((userTermination || timeTermination) && (indnum == storedTrees.size() - 1))
		outf << "[" << TerminationWarningMessage().c_str() << "]" << endl;
#endif
		
	outf.close();
	
	if(conf->outputPhylipTree){//output a phylip formatted tree if desired
		char phyname[85];
		sprintf(phyname, "%s.phy", treefname);
		ofstream phytree(phyname);
		phytree.precision(8);
		WritePhylipTree(phytree);
		phytree.close();
		}
	}

void Population::WriteStoredTrees( const char* treefname ){
	assert( treefname );

	string name;
	name = treefname;
	name += ".all.tre";
	ofstream outf( name.c_str() );
	outf.precision(8);

	dataPart->BeginNexusTreesBlock(outf);

	ofstream phytree;
	if(conf->outputPhylipTree){
		char phyname[85];
		sprintf(phyname, "%s.all.phy", treefname);
		phytree.open(phyname);
		phytree.precision(8);
		}

	int bestRep = EvaluateStoredTrees(false);
	
	Individual tempInd;
	for(unsigned r=0;r<storedTrees.size();r++){
		const Individual *curInd;
		if(Tree::outgroup != NULL || conf->collapseBranches){
			tempInd.DuplicateIndivWithoutCLAs(storedTrees[r]);
			//BMERGE TODO - still need to figure out rerooting
#ifndef USE_BEAGLE
			if(Tree::outgroup != NULL)
				OutgroupRoot(&tempInd, -1);
#endif
			if(conf->collapseBranches){
				int num = 0;
				tempInd.treeStruct->root->CollapseMinLengthBranches(num);
//				outman.UserMessage("%d COLLAPSED", num);
				}
			curInd = &tempInd;
			}
		else
			curInd = storedTrees[r];
		if(r == bestRep) 
			outf << "tree rep" << r+1 << "BEST = [&U][!GarliScore " << curInd->Fitness() << "][!GarliModel ";
		else 
			outf << "tree rep" << r+1 << " = [&U][!GarliScore " << curInd->Fitness() << "][!GarliModel ";
		string mods;
		curInd->modPart.FillGarliFormattedModelStrings(mods);
		outf << mods;
		outf << "]";

		outf.setf( ios::floatfield, ios::fixed );
		outf.setf( ios::showpoint );
		curInd->treeStruct->root->MakeNewick(treeString, false, true);
		outf << treeString << ";\n";

		if(conf->outputPhylipTree){//output a phylip formatted tree if requested
			WritePhylipTree(phytree);
			}
		}

/*
	for(unsigned r=0;r<storedTrees.size();r++){
		if(r == bestRep) 
			outf << "tree rep" << r+1 << "BEST = [&U][!GarliScore " << storedTrees[r]->Fitness() << "][!GarliModel ";
		else 
			outf << "tree rep" << r+1 << " = [&U][!GarliScore " << storedTrees[r]->Fitness() << "][!GarliModel ";
		storedTrees[r]->mod->OutputGarliFormattedModel(outf);
		outf << "]";

		outf.setf( ios::floatfield, ios::fixed );
		outf.setf( ios::showpoint );
		if(Tree::outgroup != NULL) 
			OutgroupRoot(storedTrees[r], -1);
		storedTrees[r]->treeStruct->root->MakeNewick(treeString, false, true);
		outf << treeString << ";\n";

		if(conf->outputPhylipTree){//output a phylip formatted tree if requested
			WritePhylipTree(phytree);
			}
		}
*/
	outf << "end;\n";

//	if(modSpecSet.GetModSpec(0)->IsNucleotide()){
		//add a paup block setting the model params
		//PARTITION
		//storedTrees[bestRep]->mod->OutputPaupBlockForModel(outf, name.c_str());
		if(storedTrees[bestRep]->modPart.NumModels() == 1 && storedTrees[bestRep]->modPart.GetModel(0)->IsNucleotide()){
			storedTrees[bestRep]->modPart.GetModel(0)->OutputPaupBlockForModel(outf, name.c_str());
			outf << "[!****NOTE: The model parameters loaded are the final model estimates****\n****from GARLI for the best scoring search replicate (#" << bestRep + 1 << ").****\n****The best model parameters for other trees may vary.****]" << endl;
			}
		else{
			for(int m = 0;m < storedTrees[bestRep]->modPart.NumModels();m++){
				//DEBUG
				if(storedTrees[bestRep]->modPart.GetModel(m)->IsNucleotide()){
					char mStr[20];
					sprintf(mStr, "[M%d\n", m + 1);
					outf << mStr;
					storedTrees[bestRep]->modPart.GetModel(m)->OutputPaupBlockForModel(outf, name.c_str());
					outf << "\n]\n";
					}
				}
			}
//		}
	if(userTermination || timeTermination){
		outf << "[" << TerminationWarningMessage().c_str() << "]";
		}
	outf.close();
	if(conf->outputPhylipTree) phytree.close();
	}

//CAREFUL HERE!  This function assumes the the treestring was just
//filled with MakeNewick, making a tree with taxon NUMBERS in the specification.
//This function then just reads that treestring and translates to taxon NAMES
//on the fly and outputs everything to the string passed in, which needs to 
//be already open
void Population::WritePhylipTree(ofstream &phytree){
	char *loc=treeString;
	NxsString temp;
	while(*loc){
		if(*loc == ':'){
			temp += *loc++;
			while(*loc != ',' && *loc != ')')
				temp += *loc++;
			phytree << temp.c_str();
			temp="";				
			}
		if(isdigit(*loc) == false) phytree << *loc++;
		else{
			while(isdigit(*loc))
				temp += *loc++;
			//The stored taxon names will have been gotten with GetEscaped, and thus might
			//have quotes around them if they have Nexus punctuation.  The quotes probably
			//shouldn't appear in the phylip output.  However, if the names have three single
			//quotes this corresponds to a single literal quote, in which case it will be output
			NxsString pname = dataPart->TaxonLabel(atoi(temp.c_str())-1);
			if(pname[0] == '\'' && pname[pname.size()-1] == '\''){
				pname.erase(pname.end()-1);
				pname.erase(pname.begin());
				}
			if(pname[0] == '\'' && pname[1] == '\'' ){
				pname.erase(pname.end()-1);
				pname.erase(pname.begin());
				}
			phytree << pname.c_str();
			//phytree << data->TaxonLabel(atoi(temp.c_str())-1);
			temp="";
			}
		}
	phytree << ";" << endl;
	}


char * Population::MakeNewick(int i, bool internalNodes)
{
	indiv[i].treeStruct->root->MakeNewick(treeString, internalNodes, true);
	assert(!treeString[stringSize-1]);
	return treeString;
}

//DZ 7-7 This function will get rid of multiple references to the same treeStruct
//from different individuals.  This keeps FLOAT_TYPE deletion from occuring in the destructor.
//Not the most elegant, but it works.
void Population::EliminateDuplicateTreeReferences(){

	bool dupe;
	vector<Tree *> tstructs;
	
	//go through the indiv array
	for(unsigned i=0;i<conf->nindivs;i++){
		//check if we have already encountered this treeStruct
		dupe=false;
		for(vector<Tree *>::iterator tit=tstructs.begin();tit!=tstructs.end();tit++){
			if(indiv[i].treeStruct==(*tit)){
				dupe=true;
				indiv[i].treeStruct=NULL;
				break;
				}
			}
		if(dupe==false){
			tstructs.push_back(indiv[i].treeStruct);
			}
		}
	
	//go through the newindiv array
	for(unsigned i=0;i<conf->nindivs;i++){
		//check if we have already encountered this treeStruct
		dupe=false;
		for(vector<Tree *>::iterator tit=tstructs.begin();tit!=tstructs.end();tit++){
			if(newindiv[i].treeStruct==(*tit)){
				dupe=true;
				newindiv[i].treeStruct=NULL;
				break;
				}
			}
		if(dupe==false){
			tstructs.push_back(newindiv[i].treeStruct);
			}
		}
		
	//go through the unusedTree vector
	for(vector<Tree*>::iterator vit=unusedTrees.begin();vit!=unusedTrees.end();vit++){
		dupe=false;
		for(vector<Tree *>::iterator tit=tstructs.begin();tit!=tstructs.end();tit++){
			if((*vit)==(*tit)){
				dupe=true;
				unusedTrees.erase(vit);
				vit--;
				break;
				}
			}
		}

	}

void Population::CheckAllTrees(){//debugging function
	for(unsigned i=0;i<conf->nindivs;i++){
		//check that trees are properly formed
		indiv[i].treeStruct->root->CheckforLeftandRight();
		indiv[i].treeStruct->root->CheckforPolytomies();
		indiv[i].treeStruct->root->CheckTreeFormation();
		//check that no individuals point to the same treeStruct
		for(unsigned j=i+1;j<conf->nindivs;j++)
			assert(!(indiv[i].treeStruct==indiv[j].treeStruct));
			}
		}

void Population::CheckTreesVsClaManager(){
	//go through each node for each tree and make sure that the numbers in the assignedClaArray are correct
/*	int numCopies=claMan->NumCopies();
	int numNodes=claMan->NumNodes();
	int count;
	for(int n=0;n<numNodes;n++){
		for(int c=0;c<numCopies;c++){
			count=0;
			for(int i=0;i<total_size;i++){	
				if(indiv[i].treeStruct->allNodes[claMan->ReverseConvertNodeIndex(n)]->claIndex==c) count++;
				}
			claMan->CheckAssignedNumber(count, n, c);
			}
		}
*/	}

/*		
int Population::SwapIndividuals(int n, const char* tree_strings_in, FLOAT_TYPE* kappa_probs_in, char** tree_strings_out_, FLOAT_TYPE** kappa_probs_out_)	{
	char*& tree_strings_out = *tree_strings_out_;
	FLOAT_TYPE*& kappa_probs_out = *kappa_probs_out_;
	
	int* indivs_to_send;
	GetNRandomIndivIndices(&indivs_to_send, n);
	GetSpecifiedTreeStrings(&tree_strings_out, n, indivs_to_send);
	GetSpecifiedKappas(&kappa_probs_out, n, indivs_to_send);

	
	
	// determine what to replace out (don't send out our best indiv!)
	int* indivs_to_replace = new int[n];
	for (int i = 0; i < n; ++i)	{
		if (indivs_to_send[i] == (int)cumfit[current_size-1][0])	{
			indivs_to_replace[i] = ((rand() % (current_size-1))+1 + indivs_to_send[i]) % current_size;
			assert(indivs_to_replace[i] != (int)cumfit[current_size-1][0]);
		}
		else
			indivs_to_replace[i] = indivs_to_send[i];
	}
	
	EliminateDuplicateTreeReferences();
	int x;
	const char *p = tree_strings_in;
	for (int i = 0; i < n; ++i)	{
		x = indivs_to_replace[i];
		assert(x != -1);  // make sure we're not at the end of the array
		assert(x != (int)cumfit[current_size-1][0]);  // make sure we're not replacing best
		// put in the new tree
		indiv[x].treeStruct->RemoveTreeFromAllClas();
		delete indiv[x].treeStruct;
		indiv[x].treeStruct = new Tree(p, params->data, sharedcl);
		indiv[x].treeStruct->AssignCLAsFromMaster();
		// put in the kappa prob
		indiv[x].kappa = kappa_probs_in[i];
		// set some other stuff
		indiv[x].SetDirty();
		indiv[x].parent = -1;
		p += strlen(p) + 1;
	}

	delete [] indivs_to_replace;
	delete [] indivs_to_send;
	return 0;
}
*/

/* This is all old parallel stuff not currently being used
int Population::ReplaceSpecifiedIndividuals(int count, int* which_array, const char* tree_strings, FLOAT_TYPE* model_string)	{
	//assert(count < CountTreeStrings(tree_strings)); // sanity check
	int which;
	for (int i = 0; i < count; ++i)	{
		which = which_array[i];
		Individual *ind=&indiv[which];
		ind->treeStruct->RemoveTreeFromAllClas();
		topologies[ind->topo]->RemoveInd(which);
		ind->topo=-1;
		ind->mutation_type=-1;
		
		delete ind->treeStruct;
		ind->treeStruct = new Tree(tree_strings, true);
		ind->treeStruct->AssignCLAsFromMaster();
		ind->mod->SetModel(model_string);
		ind->treeStruct->modPart=&ind->modPart;

		ind->SetDirty();
		tree_strings += strlen(tree_strings)+1;
		ind->treeStruct->modPart=&ind->modPart;
		}
	CompactTopologiesList();
	UpdateTopologyList(indiv);
	return 0;
}

int Population::GetNRandomIndivIndices(int** indiv_list, int n)	{
	int* ar = new int[total_size];
	for (unsigned i = 0; i < total_size; ++i)
		ar[i] = i;
	ScrambleArray<int>(total_size, ar);
	*indiv_list = new int[n];
	for (int i = 0; i < n; ++i)
		(*indiv_list)[i] = ar[i];
	delete [] ar;
	return 0;
}

int Population::GetNBestIndivIndices(int** indiv_list, int n)	{
	*indiv_list = new int[n];
	for (int i = 0; i < n; ++i)
		(*indiv_list)[i] = (int)cumfit[total_size-i-1][0];
	return 0;
}

int Population::GetSpecifiedTreeStrings(char** tree_strings_, int n, int* indiv_list)	{
	char*& tree_strings = *tree_strings_;
	int buf_size = 0;
	for (int i = 0; i < n; ++i)	// calc the buff size
		buf_size += (int) strlen(MakeNewick(indiv_list[i], true)) + 1;
	char* p = tree_strings = new char[buf_size+1];
	for (int i = 0; i < n; ++i)
		p += strlen(strcpy(p, MakeNewick(indiv_list[i], true))) + 1;
	*p = 0;
	assert(p-tree_strings == buf_size);
	return 0;
}

int Population::GetSpecifiedModels(FLOAT_TYPE** model_string, int n, int* indiv_list){
	FLOAT_TYPE *&model = *model_string;
	int string_size=0;
	//first calculate the appropriate size of the string and allocate it
	int nrates=modSpec->Nst()-1;
	string_size+=n*nrates;
	string_size+=n*4;//the pi's
	if(indiv[indiv_list[0]].mod->NRateCats()>1) string_size+=1*n;
#ifdef FLEX_RATES
	assert(0);
#else
	if(indiv[indiv_list[0]].mod->PropInvar()!=ZERO_POINT_ZERO) string_size+=1*n;
#endif
	model=new FLOAT_TYPE[string_size];
	
	int slot=0;
	for (int i = 0; i < n; ++i){
		//get the rates
		for(int r=0;r<nrates;r++)
			model[slot++] = indiv[indiv_list[i]].mod->Rates(r);	
		
		//get the pi's
		for(int b=0;b<4;b++)
			model[slot++] = indiv[indiv_list[i]].mod->StateFreq(b);
		
#ifdef FLEX_RATES
	assert(0);
#else
		//get alpha if we are using rate het
		if(indiv[indiv_list[0]].mod->NRateCats()>1)
			model[slot++] = indiv[indiv_list[i]].mod->Alpha();
		
		//get pinv if we are using invariant sites
		if(indiv[indiv_list[0]].mod->PropInvar()!=ZERO_POINT_ZERO)
			model[slot++] = indiv[indiv_list[i]].mod->PropInvar();		
#endif
		}
	return slot;
	}
*/

void Population::OutputLog()	{
	//log << gen << "\t" << bestFitness << "\t" << stopwatch.SplitTime() << "\t" << adap->branchOptPrecision << endl;
	if(!finishedRep) {
		log << gen << "\t" << BestFitness() << "\t" << stopwatch.SplitTime() << "\t" << adap->branchOptPrecision;
		if(conf->reportRunProgress)
			log << "\t" << 0.01 * (int) ceil(rep_fraction_done * 100) << "\t" << 0.01 * (int) ceil(tot_fraction_done * 100);

		log << endl;
#ifdef MAC_FRONTEND
		NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
		NSDictionary *progressDict = [NSDictionary dictionaryWithObjectsAndKeys:[NSNumber numberWithInt:gen], @"generation", [NSNumber numberWithDouble:BestFitness()], @"likelihood", [NSNumber numberWithInt:stopwatch.SplitTime()], @"time", [NSNumber numberWithDouble:adap->branchOptPrecision], @"precision", [NSNumber numberWithInt:lastTopoImprove], @"lastImprovement", nil];
		[[MFEInterfaceClient sharedClient] reportProgress:progressDict];
		[pool release];
#endif		
	}
	else{
		CalcAverageFitness();
		log << "Final\t" << BestFitness() << "\t" << stopwatch.SplitTime() << "\t" << adap->branchOptPrecision << endl;
		}
	}
/*
int Population::ReplicateSpecifiedIndividuals(int count, int* which, const char* tree_string, FLOAT_TYPE *model_string){
	assert(count > 0 && count <= (int)total_size);
	for (int i = 0; i < count; ++i)	{
		indiv[which[i]].treeStruct->RemoveTreeFromAllClas();
		delete indiv[which[i]].treeStruct;
		indiv[which[i]].treeStruct = new Tree(tree_string, true);
		indiv[which[i]].treeStruct->AssignCLAsFromMaster();
		indiv[which[i]].mod->SetModel(model_string);
		indiv[which[i]].treeStruct->modPart=&indiv[which[i]].modPart;
		indiv[which[i]].SetDirty();
		indiv[which[i]].treeStruct->modPart=&indiv[which[i]].modPart;
		}
	return 0;
}
*/
void Population::UpdateTreeModels(){
	for(unsigned ind=0;ind<total_size;ind++){
		newindiv[ind].treeStruct->modPart=&newindiv[ind].modPart;
//		indiv[ind].treeStruct->mod=indiv[ind].mod;
		}
	}

FLOAT_TYPE Population::IndivFitness(int i) {
	return indiv[i].Fitness();
	}

void Population::OutputModelAddresses(){
	ofstream mods("modeldeb.log", ios::app);
	
	for(unsigned i=0;i<total_size;i++){
		mods << "indiv " << i << "\t" << &indiv[i].modPart << "\t" << indiv[i].treeStruct->modPart << "\n";
		mods << "newindiv " << i << "\t" << &newindiv[i].modPart << "\t" << newindiv[i].treeStruct->modPart << "\n";
		}
	mods << endl;
	}

bool Population::NNIoptimization(unsigned indivIndex, int steps){
	Individual  currentBest;
	Individual  tempIndiv1, tempIndiv2, *best;
	int beginNode, endNode, optiNode;
	FLOAT_TYPE bestNNIFitness; 
	FLOAT_TYPE startingFitness;
	bool betterScore=false;
	
//	ofstream outf("nnidebug.tre");
//	ofstream scr("nniscores.tre");
	ofstream out;
	
	beginNode = newindiv[indivIndex].treeStruct->getNumTipsTotal() + 1;
	endNode = beginNode * 2 - 5;
	startingFitness = indiv[newindiv[indivIndex].parent].Fitness();
	bestNNIFitness = -FLT_MAX;
	
	steps = min(max(0,steps),newindiv[indivIndex].treeStruct->getNumTipsTotal()-3);
	indivIndex = min(max(0,(int)indivIndex),(int)conf->nindivs-1);

	//DJZ
	while(unusedTrees.size()<3){
		Tree *temp=new Tree();
		unusedTrees.push_back(temp);
		}
	
	tempIndiv1.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();
	tempIndiv2.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();
	currentBest.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();	
	//

	tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, &newindiv[indivIndex]);
	tempIndiv2.CopySecByRearrangingNodesOfFirst(tempIndiv2.treeStruct, &newindiv[indivIndex]);
	currentBest.CopySecByRearrangingNodesOfFirst(currentBest.treeStruct, &newindiv[indivIndex]);

	for(int i =0;i<steps;i++){
	  /*  debug information
	    cout <<"Begin exNNI mutation " <<subtreeM emberNodes.size()<<endl;
	    cout <<"Subtree node " << subtreeNode <<endl;
	    cout <<beginNode << "<->"  << endNode<<endl;
	    for(int i=0;i<subtreeM emberNodes.size() ;i++)cout << subtreeM emberNodes[i]<< "| " ;
	    cout <<endl;
	  */
	  assert(0);
/*	  
#ifdef MPI_VERSION
		for(int i=0;i<(subtreeMemberNodes.size()/2-1);i++)
		  {
		    optiNode = subtreeMemberNodes[i];
#else
		    for(optiNode=beginNode;optiNode<=endNode;optiNode++)
		      {
#endif
*/
		for(optiNode=beginNode;optiNode<=endNode;optiNode++){
			tempIndiv1.treeStruct->NNIMutate(optiNode,0, adap->branchOptPrecision, 0);
			tempIndiv2.treeStruct->NNIMutate(optiNode,1, adap->branchOptPrecision, 0);

			//newindiv[0].treeStruct->SetAllTempClasDirty();

			tempIndiv1.SetDirty();
			tempIndiv2.SetDirty();
			
			tempIndiv1.CalcFitness(0);
			tempIndiv2.CalcFitness(0);
			
			FLOAT_TYPE improvement = (FLOAT_TYPE)0.01;
			improvement = adap->recTopImproveSize;

			if(tempIndiv1.Fitness() > (bestNNIFitness) || tempIndiv2.Fitness() > (bestNNIFitness)){
				if(tempIndiv1.Fitness() > tempIndiv2.Fitness()) best=&tempIndiv1;
				else best=&tempIndiv2;
				
				bestNNIFitness = best->Fitness();
			
				currentBest.CopySecByRearrangingNodesOfFirst(currentBest.treeStruct, best, true);
				if(bestNNIFitness > startingFitness) betterScore=true;
				}
		
			//if the best tree we've found by NNI is better than what we started with, use it
			//for successive NNI attempts in this function
			if(bestNNIFitness > startingFitness){
			tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, &currentBest, true);
			tempIndiv2.CopySecByRearrangingNodesOfFirst(tempIndiv2.treeStruct, &currentBest, true);
				}
			else{//otherwise, revert to the starting tree
				tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, &newindiv[indivIndex], true);
				tempIndiv2.CopySecByRearrangingNodesOfFirst(tempIndiv2.treeStruct, &newindiv[indivIndex], true);								
				}
			} //end of loop through all possible NNIs 
		
	//copy the best tree that we found back into the population, whether or not it was better than what we
	//started with
	newindiv[indivIndex].CopySecByRearrangingNodesOfFirst(newindiv[indivIndex].treeStruct, &currentBest, true);
	} //end of loop through steps
	
	//Return the treestructs that we used temporarily back to the unused tree vector
	tempIndiv1.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(tempIndiv1.treeStruct);
	tempIndiv1.treeStruct=NULL;

	tempIndiv2.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(tempIndiv2.treeStruct);
	tempIndiv2.treeStruct=NULL;

	currentBest.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(currentBest.treeStruct);
	currentBest.treeStruct=NULL;
	
//	newindiv[indivIndex].treeStruct->SetAllTempClasDirty();
//	newindiv[indivIndex].mutation_type |= Individual::exNNI;

	return betterScore;
}

/* End of methods added by Yufeng Zhang*/

void Population::NNISpectrum(int sourceInd){
	Individual  tempIndiv1, tempIndiv2;
	int optiNode;
	FLOAT_TYPE previousFitness; 
	FLOAT_TYPE scorediff=ZERO_POINT_ZERO;
	//FLOAT_TYPE thresh=pertMan->nniAcceptThresh;

	int numNodes=indiv[sourceInd].treeStruct->getNumTipsTotal()-3;
	int *nodeArray=new int[numNodes];
	
	for(int i=0;i<numNodes;i++){
		//get all of the nodes, in order
		nodeArray[i]=numNodes+i+4;
		}

	previousFitness = indiv[sourceInd].Fitness();

	//DJZ
	while(unusedTrees.size()<3){
		Tree *temp=new Tree();
		unusedTrees.push_back(temp);
		}
	
	tempIndiv1.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();
	tempIndiv2.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();

	tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, &indiv[sourceInd]);
	tempIndiv2.CopySecByRearrangingNodesOfFirst(tempIndiv2.treeStruct, &indiv[sourceInd]);

	//make all the nodes dirty of all of the trees in the actual population, since they 
	//will be replaced by the perturbed individual and will take up valuable clas
	for(unsigned i=0;i<total_size;i++)
		indiv[i].treeStruct->MakeAllNodesDirty();

	tempGlobal=1;
	OutputFilesForScoreDebugging(&tempIndiv1, tempGlobal++);
	FLOAT_TYPE localprec;
	FLOAT_TYPE prec[7]={(FLOAT_TYPE).5, (FLOAT_TYPE).25, (FLOAT_TYPE).1, (FLOAT_TYPE).05, (FLOAT_TYPE).01, (FLOAT_TYPE).005, (FLOAT_TYPE).001};
	for(int q=0;q<7;q++){
		localprec=prec[q];

		char filename[50];
	//	sprintf(filename, "%d.%.4fscores.log", gen, localprec);
		ofstream temp(filename);
		temp.precision(12);
		temp << "start\t" << BestFitness() << "\n";


		for(int i=0;i<numNodes;i++){

			optiNode=nodeArray[i];

			tempIndiv1.treeStruct->NNIMutate(optiNode,0, localprec, 0);
			tempIndiv2.treeStruct->NNIMutate(optiNode,1, localprec, 0);

	//		tempIndiv1.SetDirty();
	//		tempIndiv2.SetDirty();
			
			tempIndiv1.SetFitness(tempIndiv1.treeStruct->lnL);
			tempIndiv2.SetFitness(tempIndiv2.treeStruct->lnL);


			temp << tempIndiv1.Fitness() << "\n" << tempIndiv2.Fitness() << "\n";
			if(q==0){
				OutputFilesForScoreDebugging(&tempIndiv1, tempGlobal++);
				OutputFilesForScoreDebugging(&tempIndiv2, tempGlobal++);
				}
			FLOAT_TYPE diff1=tempIndiv1.Fitness() - previousFitness;
			FLOAT_TYPE diff2=tempIndiv2.Fitness() - previousFitness;

			tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, &indiv[sourceInd], true);
			tempIndiv2.CopySecByRearrangingNodesOfFirst(tempIndiv2.treeStruct, &indiv[sourceInd], true);
			}		
		temp.close();
		}
	
	//Return the treestructs that we used temporarily back to the unused tree vector
	tempIndiv1.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(tempIndiv1.treeStruct);
	tempIndiv1.treeStruct=NULL;

	tempIndiv2.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(tempIndiv2.treeStruct);
	tempIndiv2.treeStruct=NULL;

	delete []nodeArray;
}

#ifdef INCLUDE_PERTURBATION
void Population::NNIPerturbation(int sourceInd, int indivIndex){
	Individual  currentBest;
	Individual  tempIndiv1, tempIndiv2, *best;
	int optiNode;
	FLOAT_TYPE previousFitness; 
//	bool betterScore=false;
	FLOAT_TYPE scorediff=ZERO_POINT_ZERO;
	FLOAT_TYPE thresh=pertMan->nniAcceptThresh;
	int nummoves=0;
	
	ofstream out;
	
//	int numNodes=indiv[indivIndex].treeStruct->getNumTipsTotal()-3;
//	int *nodeArray=new int[numNodes];
/*	
	for(int i=0;i<numNodes;i++){
		nodeArray[i]=indiv[indivIndex].treeStruct->GetRandomInternalNode();
		//get all of the nodes, in order
//		nodeArray[i]=numNodes+i+4;
		}
*/
//	ScrambleArray(numNodes, nodeArray);

	previousFitness = indiv[sourceInd].Fitness();

	//DJZ
	while(unusedTrees.size()<3){
		Tree *temp=new Tree();
		unusedTrees.push_back(temp);
		}
	
	tempIndiv1.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();
	tempIndiv2.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();
	currentBest.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();	
	//

	tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, &indiv[sourceInd]);
	tempIndiv2.CopySecByRearrangingNodesOfFirst(tempIndiv2.treeStruct, &indiv[sourceInd]);
	currentBest.CopySecByRearrangingNodesOfFirst(currentBest.treeStruct, &indiv[sourceInd]);

	int n=2;
	
	//make all the nodes dirty of all of the trees in the actual population, since they 
	//will be replaced by the perturbed individual and will take up valuable clas
	for(unsigned i=0;i<total_size;i++)
		indiv[i].treeStruct->MakeAllNodesDirty();

	char filename[50];
	if(rank < 10)
		sprintf(filename, "pertreport0%d.log", rank);
	else 
		sprintf(filename, "pertreport%d.log", rank);
	ofstream pert(filename, ios::app);
	pert.precision(10);
	pert << "gen\t" << gen << "\tstart\t" << BestFitness() << "\n";

	outman.UserMessage("Performing NNI Perturbation.  Starting score= %.4f", BestFitness());

	
/*	char filename[50];
	FLOAT_TYPE localprec=.5;
	sprintf(filename, "%d.%.4fscores.log", gen, localprec);
	ofstream temp(filename);
	temp.precision(12);
	temp << "start\t" << BestFitness() << "\n";
*/

//	for(int i=0;i<numNodes;i++){
	int attempts, accepts;
	for(accepts=0, attempts=0;(accepts<pertMan->nniTargetAccepts) && (attempts <= pertMan->nniMaxAttempts);){
		if(! (attempts++ % (pertMan->nniMaxAttempts/20))) outman.UserMessage(".");
		outman.flush();

		optiNode=indiv[indivIndex].treeStruct->GetRandomInternalNode();
//		optiNode=nodeArray[i];

//		tempIndiv1.treeStruct->NNIMutate(optiNode,0, localprec, 0);
//		tempIndiv2.treeStruct->NNIMutate(optiNode,1, localprec, 0);

		tempIndiv1.treeStruct->NNIMutate(optiNode,0, adap->branchOptPrecision, 0);
		tempIndiv2.treeStruct->NNIMutate(optiNode,1, adap->branchOptPrecision, 0);

//		tempIndiv1.SetDirty();
//		tempIndiv2.SetDirty();
		
		tempIndiv1.SetFitness(tempIndiv1.treeStruct->lnL);
		tempIndiv2.SetFitness(tempIndiv2.treeStruct->lnL);


//		temp << tempIndiv1.Fitness() << "\n" << tempIndiv2.Fitness() << "\n";

		FLOAT_TYPE diff1=tempIndiv1.Fitness() - previousFitness;
		FLOAT_TYPE diff2=tempIndiv2.Fitness() - previousFitness;

		//ignore NNI's that improve the fitness, because they are probably just undoing a previous NNI
//		if(((diff1 < ZERO_POINT_ZERO) && (diff1 + thresh > ZERO_POINT_ZERO)) || ((diff2 < ZERO_POINT_ZERO) && (diff2 + thresh > ZERO_POINT_ZERO))){
//			if((diff1 < ZERO_POINT_ZERO) && ((diff1 > diff2) || (diff2 >= ZERO_POINT_ZERO))) best=&tempIndiv1;
		if(diff1 < ZERO_POINT_ZERO || diff2 < ZERO_POINT_ZERO){
			if((diff1 < ZERO_POINT_ZERO) && ((diff1 > diff2) || (diff2 >= ZERO_POINT_ZERO))) best=&tempIndiv1;
			else best=&tempIndiv2;
			
			FLOAT_TYPE acceptanceProb=exp(-conf->selectionIntensity * (previousFitness - best->Fitness()));
			if(rnd.uniform() < acceptanceProb){
				FLOAT_TYPE thisdiff=best->Fitness() - previousFitness;
				assert(thisdiff < 0);
				scorediff += thisdiff;
				accepts++;
				previousFitness = best->Fitness();
				pert << accepts << "\t" << optiNode << "\t" << thisdiff << "\n";
			
				currentBest.CopySecByRearrangingNodesOfFirst(currentBest.treeStruct, best, true);
				}
			}

		tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, &currentBest, true);
		tempIndiv2.CopySecByRearrangingNodesOfFirst(tempIndiv2.treeStruct, &currentBest, true);
		}
		
	indiv[indivIndex].CopySecByRearrangingNodesOfFirst(indiv[indivIndex].treeStruct, &currentBest, true);
	
	//Return the treestructs that we used temporarily back to the unused tree vector
	tempIndiv1.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(tempIndiv1.treeStruct);
	tempIndiv1.treeStruct=NULL;

	tempIndiv2.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(tempIndiv2.treeStruct);
	tempIndiv2.treeStruct=NULL;

	currentBest.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(currentBest.treeStruct);
	currentBest.treeStruct=NULL;
	
  	UpdateTopologyList(indiv);
	indiv[indivIndex].SetDirty();
	indiv[indivIndex].CalcFitness(0);
	AssignNewTopology(indiv, indivIndex);
	UpdateTopologyList(indiv);
	
	SetNewBestIndiv(indivIndex);
	indiv[indivIndex].treeStruct->calcs=calcCount;
	calcCount=0;
	indiv[indivIndex].mutation_type=-1;
	pert << "end\t" << indiv[indivIndex].Fitness() << "\n";
	
	outman.UserMessage("Completed Perturbation.\n  %d NNI's accepted in %d attempts. Current score= %.4f", accepts, attempts, BestFitness());
	
//	delete []nodeArray;
}

void Population::TurnOffRatchet(){
	data->RestoreOriginalCounts();
	pertMan->ratcheted=false;
	
	claMan->MakeAllHoldersDirty();
	for(unsigned i=0;i<total_size;i++) indiv[i].SetDirty();
	CalcAverageFitness();
	bestFitness=BestFitness();
	pertMan->lastPertGeneration=gen;
	adap->reset=true;
	outman.UserMessage("Returning to normal character weighting...");
	char filename[50];
	if(rank < 10)
		sprintf(filename, "pertreport0%d.log", rank);
	else 
		sprintf(filename, "pertreport%d.log", rank);
	ofstream pert(filename, ios::app);
	pert << "Returning to normal character weighting..." << endl;
	pert.close();
	}

void Population::RestoreAllTimeBest(){
	UpdateTopologyList(indiv);
	topologies[indiv[0].topo]->RemoveInd(0);
	CompactTopologiesList();
	indiv[0].CopySecByRearrangingNodesOfFirst(indiv[0].treeStruct, allTimeBest, true);
	indiv[0].treeStruct->AssignCLAsFromMaster();
	AssignNewTopology(indiv, 0);
	indiv[0].CalcFitness(0);
	SetNewBestIndiv(0);
	FillPopWithClonesOfBest();
	CalcAverageFitness();
	}

void Population::RestoreBestForPert(){
	UpdateTopologyList(indiv);
	topologies[indiv[0].topo]->RemoveInd(0);
	CompactTopologiesList();
	indiv[0].CopySecByRearrangingNodesOfFirst(indiv[0].treeStruct, bestSinceRestart, true);
	indiv[0].treeStruct->AssignCLAsFromMaster();
	AssignNewTopology(indiv, 0);
	indiv[0].CalcFitness(0);
	SetNewBestIndiv(0);
	FillPopWithClonesOfBest();
	CalcAverageFitness();
	char filename[50];
	if(rank < 10)
		sprintf(filename, "pertreport0%d.log", rank);
	else 
		sprintf(filename, "pertreport%d.log", rank);
	ofstream pert(filename, ios::app);
	pert.precision(10);
	pert << "restoring best individual with score of " << indiv[0].Fitness() << "\n";
	pert.close();
	outman.UserMessage("restoring best individual with score of %.4f\n %d perturbation(s) performed without improvement.", bestSinceRestart->Fitness(), pertMan->numPertsNoImprove);
	}

void Population::StoreBestForPert(){
	if(BestFitness() > allTimeBest->Fitness()) StoreAllTimeBest();
	
	if(bestSinceRestart->treeStruct==NULL){
		if(unusedTrees.empty()){
			Tree *temp=new Tree();
			unusedTrees.push_back(temp);
			}					
		bestSinceRestart->treeStruct=*(unusedTrees.end()-1);
		unusedTrees.pop_back();
		}
	bestSinceRestart->CopySecByRearrangingNodesOfFirst(bestSinceRestart->treeStruct, &indiv[bestIndiv]);
	bestSinceRestart->topo=-1;
	//need to do this to be sure that the bestSinceRestart isn't tying up clas
	bestSinceRestart->treeStruct->RemoveTreeFromAllClas();
	
	char filename[50];
	if(rank < 10)
		sprintf(filename, "pertreport0%d.log", rank);
	else 
		sprintf(filename, "pertreport%d.log", rank);
	ofstream pert(filename, ios::app);
	pert.precision(10);
	pert << "storing best individual with score of " << bestSinceRestart->Fitness() << "\n";
	pert.close();
	outman.UserMessage("storing best individual with score of %.4f", bestSinceRestart->Fitness());
	}

void Population::StoreAllTimeBest(){
	WriteTreeFile(besttreefile);
	if(allTimeBest->treeStruct==NULL){
		if(unusedTrees.empty()){
			Tree *temp=new Tree();
			unusedTrees.push_back(temp);
			}					
		allTimeBest->treeStruct=*(unusedTrees.end()-1);
		unusedTrees.pop_back();
		}
	allTimeBest->CopySecByRearrangingNodesOfFirst(allTimeBest->treeStruct, &indiv[bestIndiv]);
	allTimeBest->topo=-1;
	//need to do this to be sure that the alltimebest isn't tying up clas
	allTimeBest->treeStruct->RemoveTreeFromAllClas();
	}
#endif

void Population::keepTrack(){

	if(((gen-1)%adap->intervalLength)==0){
		if(gen>1) adap->PrepareForNextInterval();
		}

	//remember that the indiv and newindiv arrays have already been switched, so the newindivs are the parents of the indivs

	//adap->reset will be true if we've Ratcheted, in which case the scores will be noncomparable
	//so reset these values
	if(adap->reset==true){
		adap->laststepscore=BestFitness();
		adap->lastgenscore=BestFitness();
		adap->reset=false;
		}
	
	if(gen==1)
		adap->lastgenscore = adap->laststepscore = newindiv[0].Fitness();
		
	adap->improvetotal[0] = BestFitness() - adap->laststepscore;

	for(unsigned i=0;i<conf->nindivs;i++){
//		FLOAT_TYPE scoreDif=indiv[i].Fitness() - newindiv[indiv[i].parent].Fitness();
		FLOAT_TYPE scoreDif=indiv[i].Fitness() - adap->lastgenscore;
		int typ = indiv[i].mutation_type;
		if(typ > 0){
			if(scoreDif>0){
	//			if(i==bestIndiv){
					//keep track of when the last significant beneficial topo mutation occured
					//this will be used for the stopping criterion, precision reduction and update reduction in the parallel version
#ifndef NO_EVOLUTION

#ifdef IGNORE_SMALL_TOPO_IMP
					if(typ&Individual::anyTopo){
						if(i == bestIndiv && scoreDif > significantTopoChange){
							indiv[0].treeStruct->attemptedSwaps.ClearAttemptedSwaps();
							indiv[0].treeStruct->CalcBipartitions(true);
							indiv[i].treeStruct->CalcBipartitions(true);
							if(indiv[0].treeStruct->IdenticalTopology(indiv[i].treeStruct->root)==false){
								lastTopoImprove=gen;
								if(i == bestIndiv){
									AppendTreeToTreeLog(indiv[bestIndiv].mutation_type);
									}
								}
							}
						else{//just ignore this small improvement.  Kill the individual's chance							
							//of reproducing
							FLOAT_TYPE scr=indiv[i].Fitness();
							indiv[i].SetFitness(-FLT_MAX);
							indiv[i].treeStruct->lnL = -FLT_MAX;
							CalcAverageFitness();
							bestIndiv=(int) cumfit[3][0];
							indiv[i].SetFitness(scr);
							}
						}

#else
					if(typ&Individual::anyTopo || adap->topoWeight==ZERO_POINT_ZERO){
						//clearing of the swaps records needs to be done for _any_ new best topo, not
						//just ones that are significantly better
						indiv[0].treeStruct->CalcBipartitions(true);
						indiv[i].treeStruct->CalcBipartitions(true);
						bool sameTopo = indiv[0].treeStruct->IdenticalTopologyAllowingRerooting(indiv[i].treeStruct->root);
						//if this is a new best and isn't the same topology, clear the swaplist
						if(i == bestIndiv && sameTopo == false)
							indiv[0].treeStruct->attemptedSwaps.ClearAttemptedSwaps();

						if(scoreDif > conf->significantTopoChange){
							//if this is a new best, it is a different topology and it is significantly better 
							//update the lastTopoImprove
							if(sameTopo == false){
								lastTopoImprove=gen;
								if(i == bestIndiv){
									AppendTreeToTreeLog(indiv[bestIndiv].mutation_type);
									}
								}
							else if(adap->topoWeight==ZERO_POINT_ZERO){
								if(i == bestIndiv){
									AppendTreeToTreeLog(indiv[bestIndiv].mutation_type);
									}
								}
							}
					}
#endif
#endif
					
					if(typ&(Individual::randNNI)){
						adap->randNNI[0] += scoreDif;
						}
//					if(typ&(Individual::exNNI)) 		adap->exNNI[0] += scoreDif;
					if(typ&(Individual::randSPR)) 	adap->randSPR[0] += scoreDif;
					if(typ&(Individual::limSPR)) 		adap->limSPR[0] += scoreDif;
//					if(typ&(Individual::exlimSPR)) 	 	adap->exlimSPR[0] += scoreDif;
	#ifdef GANESH
					if(typ&(Individual::randPECR)) 	 	adap->randPECR[0] += scoreDif;
	#endif
		//			if(typ&(Individual::taxonSwap)) 	adap->taxonSwap[0] += scoreDif; 
					if(typ == (Individual::brlen)) 		adap->onlyBrlen[0] += scoreDif;
					if(typ&(Individual::bipartRecom)) adap->bipartRecom[0] += scoreDif;
					if(typ&(Individual::randRecom)) 	adap->randRecom[0] += scoreDif;
					if(typ&(Individual::anyModel)) 	adap->anyModel[0] += scoreDif;
					
	#ifdef MPI_VERSION
					if(scoreDif > adap->branchOptPrecision){
						if(typ&(Individual::bipartRecom)) adap->bestFromRemote[0] += scoreDif;
						if(typ&(Individual::bipartRecom)) adap->bestFromRemoteNum[0] ++;
						if(typ&(Individual::subtreeRecom)) adap->bestFromRemote[0] += scoreDif;
						if(typ&(Individual::subtreeRecom)) adap->bestFromRemoteNum[0] ++;
						}
	#endif
	//				}
				}
			if(typ&(Individual::randNNI)) 	adap->randNNInum[0]++;
//			if(typ&(Individual::exNNI)) 	 	adap->exNNInum[0]++;
			if(typ&(Individual::randSPR)) 	adap->randSPRnum[0]++;
			if(typ&(Individual::limSPR))	 	adap->limSPRnum[0]++;
//			if(typ&(Individual::exlimSPR)) 	adap->exlimSPRnum[0]++;
	#ifdef GANESH
			if(typ&(Individual::randPECR)) 	adap->randPECRnum[0]++;
	#endif
		//	if(typ&(Individual::taxonSwap)) 	adap->taxonSwapnum[0]++;
			if(typ == (Individual::brlen)) 	adap->onlyBrlennum[0]++;
			if(typ&(Individual::bipartRecom)) adap->bipartRecomnum[0]++;
			if(typ&(Individual::randRecom)) 	adap->randRecomnum[0]++;
			if(typ&(Individual::anyModel))		adap->anyModelnum[0]++;
			}
		}

	adap->lastgenscore=BestFitness();

	//things to do on the final generation of an interval	
	if(gen%adap->intervalLength==0){
		//improveOverStoredIntervals is only used on generations that are multiples of intervalLength
		//so it won't contain the improvement in the latest interval until it's end
		adap->improveOverStoredIntervals=ZERO_POINT_ZERO;
		for(unsigned i=0;i<adap->intervalsToStore;i++)
			adap->improveOverStoredIntervals += adap->improvetotal[i];
		if(adap->improveOverStoredIntervals < ZERO_POINT_ZERO) adap->improveOverStoredIntervals = ZERO_POINT_ZERO;
		//update the mutation probailities
		if(gen>=(adap->intervalLength*adap->intervalsToStore*ZERO_POINT_FIVE)){
			adap->UpdateProbs();
			if(conf->outputMostlyUselessFiles) adap->OutputProbs(probLog, gen);
			}
		adap->laststepscore=BestFitness();
		}
	}

int ParallelManager::DetermineSubtrees(Tree *tr, ofstream &scr){
	//Determine what the best node we could choose to be the root would be in terms of 
	//the partitioning efficiency

	TreeNode *nd=tr->root;

	int bestRoot=0, orphans=ntax;
	bool done=false;
	FLOAT_TYPE bestScore=ZERO_POINT_ZERO;
	FLOAT_TYPE thisScore;

	ClearSubtrees();
	
	//trying new partitioning function
	int one=0, two=0, three=0;
	vector<Subtree *> sub1, sub2, sub3;
	if(nd->left->left) NewPartition(nd->left, one, sub1);
	if(nd->left->next->left) NewPartition(nd->left->next, two, sub2);
	if(nd->right->left) NewPartition(nd->right, three, sub3);
	for(vector<Subtree*>::iterator it = sub1.begin();it!=sub1.end();it++){
		bestScore += (*it)->score;
		(*it)->Log(scr);
		orphans -= (*it)->taxa;
		delete *it;
   		}
	for(vector<Subtree*>::iterator it = sub2.begin();it!=sub2.end();it++){
		bestScore += (*it)->score;
		(*it)->Log(scr);
		orphans -= (*it)->taxa;
		delete *it;
   		}
	for(vector<Subtree*>::iterator it = sub3.begin();it!=sub3.end();it++){
		bestScore += (*it)->score;
		(*it)->Log(scr);
		orphans -= (*it)->taxa;
		delete *it;
   		}
   	bestScore += pow((FLOAT_TYPE)(orphans), orphanFactor);

	scr << "root " << nd->nodeNum << " score " << bestScore << " orphans " << orphans << endl;
   		
   	sub1.clear();
   	sub2.clear();
	sub3.clear();
	
	nd=nd->left;

	while(!done){

		if(nd->CountTerminals(0) > minSubtreeSize){
			orphans=ntax;
			//ClearSubtrees();
			if(nd->left->left) NewPartition(nd->left, one, sub1);
			if(nd->right->left) NewPartition(nd->right, two, sub2);
			NewPartitionDown(nd->anc, nd, three, sub3);
			
			thisScore=ZERO_POINT_ZERO;
			for(vector<Subtree*>::iterator it = sub1.begin();it!=sub1.end();it++){
				thisScore += (*it)->score;
				(*it)->Log(scr);
				orphans -= (*it)->taxa;
				delete *it;
		   		}
			for(vector<Subtree*>::iterator it = sub2.begin();it!=sub2.end();it++){
				thisScore += (*it)->score;
				(*it)->Log(scr);
				orphans -= (*it)->taxa;
				delete *it;
		   		}
			for(vector<Subtree*>::iterator it = sub3.begin();it!=sub3.end();it++){
				thisScore += (*it)->score;
				(*it)->Log(scr);
				orphans -= (*it)->taxa;		
				delete *it;
		   		}
		   	thisScore+=pow((FLOAT_TYPE)(orphans), orphanFactor);
		   	scr << "root " << nd->nodeNum << " score " << thisScore << " orphans " << orphans << endl;
		   	sub1.clear();
		   	sub2.clear();
			sub3.clear();			
		
			if(thisScore < bestScore){
				bestScore=thisScore;
				bestRoot=nd->nodeNum;
				}
			}
		
		if(nd->left != NULL){
			nd=nd->left;
			}
		else if(nd->next != NULL){
			nd=nd->next;
			}
		else{
			while(nd->anc!=NULL){
				nd=nd->anc;
				if(nd->next != NULL){
					nd=nd->next;
					break;
					}
				}
			if(nd->anc==NULL){
				done=true;
				}
			}
		}
	return bestRoot;
/*	
	if(nd->left->CountBranches(0)>1)
		Partition(nd->left);
	if(nd->left->next->CountBranches(0)>1)
		Partition(nd->left->next);
	if(nd->right->CountBranches(0)>1)
		Partition(nd->right);
	FLOAT_TYPE bestScore=ScorePartitioning(nd->nodeNum, pscores);	
	FLOAT_TYPE thisScore;
	
	nd=nd->left;

	while(!done){
		if(nd->CountTerminals(0) > minSubtreeSize){
			ClearSubtrees();
			Partition(nd->left);		
			Partition(nd->right);
			PartitionDown(nd->anc, nd);
				
			thisScore=ScorePartitioning(nd->nodeNum, pscores);
			
			if(thisScore < bestScore){
				bestScore=thisScore;
				bestRoot=nd->nodeNum;
				}
			}
		
		if(nd->left != NULL){
			nd=nd->left;
			}
		else if(nd->next != NULL){
			nd=nd->next;
			}
		else{
			while(nd->anc!=NULL){
				nd=nd->anc;
				if(nd->next != NULL){
					nd=nd->next;
					break;
					}
				}
			if(nd->anc==NULL){
				done=true;
				}
			}
		}
	return bestRoot;
*/	}
	
	
void Population::StartSubtreeMode(){
	OutputFate();
	gen++;

	bool subtreesOK=false;
	int attempt=1;
	
	int origMinSubtreeSize=paraMan->minSubtreeSize;
	int origTargetSubtreeSize=paraMan->targetSubtreeSize;
	FLOAT_TYPE origOrphanFactor=paraMan->orphanFactor;

	ofstream pscores("partscores.log", ios::app);

	do{
		pscores << "gen " << gen << " attempt " <<  attempt << endl;
		int bestRoot=paraMan->DetermineSubtrees(indiv[bestIndiv].treeStruct, pscores);	
		//now we need to reroot to the best root found

		pscores << "best root=" << bestRoot << "\n";

		if(bestRoot!=0){
			indiv[bestIndiv].treeStruct->RerootHere(bestRoot);
			}

		//mark all of the trees inaccurate
		for(unsigned i=0;i<total_size;i++){
			indiv[i].accurateSubtrees=false;
			newindiv[i].accurateSubtrees=false;
			}

		indiv[bestIndiv].accurateSubtrees=true;
		FillPopWithClonesOfBest();

		TreeNode *nd=indiv[bestIndiv].treeStruct->root;

		int orphans=paraMan->ntax;
		bool done=false;
		FLOAT_TYPE bestScore=ZERO_POINT_ZERO;

		//trying new partitioning function
		int one=0, two=0, three=0;
		vector<Subtree *> sub1, sub2, sub3;
		if(nd->left->left) paraMan->NewPartition(nd->left, one, sub1);
		if(nd->left->next->left) paraMan->NewPartition(nd->left->next, two, sub2);
		if(nd->right->left) paraMan->NewPartition(nd->right, three, sub3);
		for(vector<Subtree*>::iterator it = sub1.begin();it!=sub1.end();it++){
			bestScore += (*it)->score;
			orphans -= (*it)->taxa;
			paraMan->subtrees.push_back(*it);
	   		}
		for(vector<Subtree*>::iterator it = sub2.begin();it!=sub2.end();it++){
			bestScore += (*it)->score;
			orphans -= (*it)->taxa;
			paraMan->subtrees.push_back(*it);
	   		}
		for(vector<Subtree*>::iterator it = sub3.begin();it!=sub3.end();it++){
			bestScore += (*it)->score;
			orphans -= (*it)->taxa;
			paraMan->subtrees.push_back(*it);
	   		}
	   		
	   	sub1.clear();
	   	sub2.clear();
		sub3.clear();
		
		if((int)paraMan->subtrees.size() > paraMan->nremotes){
			paraMan->targetSubtreeSize = (int) (paraMan->targetSubtreeSize * 1.05);
			attempt++;
			pscores << "too many subtrees, increasing target size to paraMan->targetSubtreeSize..." << endl;
			}
		else subtreesOK=true;
		}while(subtreesOK==false);
	
	//set these back to their original values
	paraMan->minSubtreeSize=origMinSubtreeSize;
	paraMan->targetSubtreeSize=origTargetSubtreeSize;
	paraMan->orphanFactor=origOrphanFactor;	
	
	paraMan->PrepareForSubtreeMode(&indiv[bestIndiv], gen);
	
	if(paraMan->fewNonSubtreeNodes==true) AssignSubtree(paraMan->ChooseSubtree(), bestIndiv);

	#ifdef MASTER_DOES_SUBTREE
	AssignSubtree(paraMan->subtrees[(int)(rnd.uniform()*paraMan->subtrees.size())]->nodeNum, bestIndiv);
	#endif
	}

void Population::StopSubtreeMode(){
	paraMan->subtreeModeActive=false;
	paraMan->subtreeDefGeneration = gen;
	for(unsigned i=0;i<total_size;i++){
		indiv[i].accurateSubtrees=false;
		newindiv[i].accurateSubtrees=false;
		}
		
	#ifdef MASTER_DOES_SUBTREE
	AssignSubtree(0, bestIndiv);
	#endif
	}


void ParallelManager::PrepareForSubtreeMode(Individual *ind, int gen){
	needUpdate=false;
	beforeFirstSubtree=false;
	
	for(int i=1;i<nremotes+1;i++){
		remoteSubtreeAssign[i]=0;
		localSubtreeAssign[i]=0;
		}
	subtreeDefNumber++;
	subtreeDefGeneration=lastFullRecom=gen;

	//put the nodes that aren't in any of the subtrees into a vector 
	//that is a member of subMan, so that the master can work on
	//them itself
	nonSubtreeNodesforNNI.clear();
	nonSubtreeNodesforSPR.clear();
	FindNonSubtreeNodes(ind->treeStruct->root);
	sort(nonSubtreeNodesforNNI.begin(), nonSubtreeNodesforNNI.end());
	sort(nonSubtreeNodesforSPR.begin(), nonSubtreeNodesforSPR.end());

	if(nonSubtreeNodesforSPR.size() < 20) fewNonSubtreeNodes=true;

	subtreeDefScore=ind->Fitness();
	}

//Partition with loose lower bound and upper bound

void ParallelManager::Partition(TreeNode *pointer){

  if(pointer->left==NULL) return;

  int largestOrphan=5;
//  int min=20;
//  int min = ntax < 100 ? 10 : 20;
  //int target=params->data->NTax()/subMan->nremotes;
  
	int n1 = pointer->left->CountTerminals(0);
	int n2 = pointer->right->CountTerminals(0);
	int n  = n1 + n2;
	if(n<minSubtreeSize) return;
	if( (n1>largestOrphan && n1<minSubtreeSize) || (n2>largestOrphan && n2<minSubtreeSize) || (n<targetSubtreeSize && n>=minSubtreeSize)){
 		Subtree *st = new Subtree(pointer->nodeNum, n, pointer->dlen, ZERO_POINT_ZERO);
		subtrees.push_back(st);
  		}
	else{
    	if(n1>=minSubtreeSize) Partition(pointer->left);
	    if(n2>=minSubtreeSize) Partition(pointer->right);
  		}
	}

void ParallelManager::NewPartition(TreeNode *pointer, int &orphans, vector<Subtree*> &subtreesAbove){
	vector<Subtree*> subtreesUpLeft, subtreesUpRight;
	FLOAT_TYPE scoreAbove=ZERO_POINT_ZERO;
	int orphansHere=0, orphansLeft=0, orphansRight=0;
	
	int n1 = pointer->left->CountTerminals(0);
	int n2 = pointer->right->CountTerminals(0);
	int n  = n1 + n2;
	if(n<minSubtreeSize) return;

   	if(n1>=minSubtreeSize){
		NewPartition(pointer->left, orphansLeft, subtreesUpLeft);
   		for(vector<Subtree*>::iterator it = subtreesUpLeft.begin();it!=subtreesUpLeft.end();it++){
   			subtreesAbove.push_back(*it);
   			scoreAbove += (*it)->score;
   			}
   		subtreesUpLeft.clear();
   		}
   	else orphansHere += n1;
    if(n2>=minSubtreeSize){
    	NewPartition(pointer->right, orphansRight, subtreesUpRight);
   		for(vector<Subtree*>::iterator it = subtreesUpRight.begin();it!=subtreesUpRight.end();it++){
   			subtreesAbove.push_back(*it);
   			scoreAbove += (*it)->score;
   			}   			   			
   		subtreesUpRight.clear();
   		}
   	else orphansHere += n2;
   	
   	orphans=orphansLeft+orphansRight+orphansHere;   	
	scoreAbove += pow((FLOAT_TYPE)orphans, orphanFactor);
	
	FLOAT_TYPE scoreHere = pow((FLOAT_TYPE)(targetSubtreeSize-n), 2);
	
	if(scoreAbove > scoreHere){
		for(vector<Subtree*>::iterator it = subtreesAbove.begin();it!=subtreesAbove.end();it++){
			Subtree *del=*it;
			delete del;
			}
		subtreesAbove.clear();
		Subtree *poo=new Subtree(pointer->nodeNum, n, pointer->dlen, pow((FLOAT_TYPE)(targetSubtreeSize-n), 2));
		subtreesAbove.push_back(poo);
		orphans=0;
		}
	}

void ParallelManager::NewPartitionDown(TreeNode *pointer, TreeNode *calledFrom, int &orphans, vector<Subtree*> &subtreesAbove){
	int n, n1, n2;
	TreeNode *sib, *anc;
	vector<Subtree*> subtreesUpLeft, subtreesUpRight;
	FLOAT_TYPE scoreAbove=ZERO_POINT_ZERO;
	int orphansHere=0, orphansLeft=0, orphansRight=0;

	if(pointer->nodeNum != 0){
	  if(pointer->left==calledFrom) sib=pointer->right;
	  else sib=pointer->left;
	  anc=pointer->anc;

	  n1 = sib->CountTerminals(0);
	  n2 = anc->CountTerminalsDown(0, pointer);
	  n  = n1 + n2;
	  if(n<minSubtreeSize) return;
	  
	   	if(n1>=minSubtreeSize){
			NewPartition(sib, orphansLeft, subtreesUpLeft);
	   		for(vector<Subtree*>::iterator it = subtreesUpLeft.begin();it!=subtreesUpLeft.end();it++){
	   			subtreesAbove.push_back(*it);
	   			scoreAbove += (*it)->score;
	   			}
	   		subtreesUpLeft.clear();
	   		}
	   	else orphansHere += n1;
	    if(n2>=minSubtreeSize){
	    	NewPartitionDown(anc, pointer, orphansRight, subtreesUpRight);
	   		for(vector<Subtree*>::iterator it = subtreesUpRight.begin();it!=subtreesUpRight.end();it++){
	   			subtreesAbove.push_back(*it);
	   			scoreAbove += (*it)->score;
	   			}   			   			
	   		subtreesUpRight.clear();
	   		}
	   	else orphansHere += n2;
	   	orphans=orphansLeft+orphansRight+orphansHere;   	
		scoreAbove += pow((FLOAT_TYPE)orphans, orphanFactor);
		
		}
	else{
		TreeNode *nd1, *nd2;
	
		if(pointer->left==calledFrom){
			nd1=pointer->left->next;
			nd2=pointer->right;
			}
		else if(pointer->left->next==calledFrom){
			nd1=pointer->left;
			nd2=pointer->right;
			}
		else {
			nd1=pointer->left;
			nd2=pointer->left->next;
			}
		
	  n1 = nd1->CountTerminals(0);
	  n2 = nd2->CountTerminals(0);
	  n  = n1 + n2;
	  if(n<minSubtreeSize) return;
	  	
	   	if(n1>=minSubtreeSize){
			NewPartition(nd1, orphansLeft, subtreesUpLeft);
	   		for(vector<Subtree*>::iterator it = subtreesUpLeft.begin();it!=subtreesUpLeft.end();it++){
	   			subtreesAbove.push_back(*it);
	   			scoreAbove += (*it)->score;
	   			}
	   		subtreesUpLeft.clear();
	   		}
	   	else orphansHere += n1;
	    if(n2>=minSubtreeSize){
	    	NewPartition(nd2, orphansRight, subtreesUpRight);
	   		for(vector<Subtree*>::iterator it = subtreesUpRight.begin();it!=subtreesUpRight.end();it++){
	   			subtreesAbove.push_back(*it);
	   			scoreAbove += (*it)->score;
	   			}   			   			
	   		subtreesUpRight.clear();
	   		}
	   	else orphansHere += n2;
	   	orphans=orphansLeft+orphansRight+orphansHere;   	
		scoreAbove += pow((FLOAT_TYPE)orphans, orphanFactor);
		}			

	FLOAT_TYPE scoreHere = (FLOAT_TYPE)(targetSubtreeSize-n)*(targetSubtreeSize-n);
	
	if(scoreAbove > scoreHere){
		for(vector<Subtree*>::iterator it = subtreesAbove.begin();it!=subtreesAbove.end();it++){
			delete *it;
			}
		subtreesAbove.clear();
		Subtree *poo=new Subtree(pointer->nodeNum, n, calledFrom->dlen, pow((FLOAT_TYPE)(targetSubtreeSize-n), 2));
		subtreesAbove.push_back(poo);
		orphans=0;
		}
	}

void ParallelManager::PartitionDown(TreeNode *pointer, TreeNode *calledFrom){

//  int min=20;
  int min = ntax < 100 ? 10 : 20;
  int largestOrphan=5;
//  int target=params->data->NTax()/subMan->nremotes;
#ifdef MPI_VERSION
	int target=2*ntax/nremotes;
#else
	int target=2*ntax/9;
#endif
  TreeNode *sib, *anc;
  
  if(pointer->nodeNum != 0){
	  
	  if(pointer->left==calledFrom) sib=pointer->right;
	  else sib=pointer->left;
	  anc=pointer->anc;
	    
	  int n1 = sib->CountTerminals(0);
	  int n2 = anc->CountTerminalsDown(0, pointer);
	  int n  = n1 + n2;
	  
	  if(n<min) return;
	  

	if( (n1>largestOrphan && n1<min) || (n2>largestOrphan && n2<min) || (n<target && n>=min)){
	 	Subtree *st = new Subtree(pointer->nodeNum, n, calledFrom->dlen, ZERO_POINT_ZERO);
 		subtrees.push_back(st);
		}
	  else{
	    if(n1>=min) Partition(sib);
	    if(n2>=min) PartitionDown(anc, pointer);
	  }
	}

	else{
		TreeNode *nd1, *nd2;
	
		if(pointer->left==calledFrom){
			nd1=pointer->left->next;
			nd2=pointer->right;
			}
		else if(pointer->left->next==calledFrom){
			nd1=pointer->left;
			nd2=pointer->right;
			}
		else {
			nd1=pointer->left;
			nd2=pointer->left->next;
			}
		
	  int n1 = nd1->CountTerminals(0);
	  int n2 = nd2->CountTerminals(0);
	  int n  = n1 + n2;
	  if(n<min) return;
	  
	if( (n1>largestOrphan && n1<min) || (n2>largestOrphan && n2<min) || (n<target && n>=min)){
	 	Subtree *st = new Subtree(pointer->nodeNum, n, calledFrom->dlen, ZERO_POINT_ZERO);
	 	subtrees.push_back(st);
	  }
	  else{
	    if(n1>=min) Partition(nd1);
	    if(n2>=min) Partition(nd2);
	  }	
	}
}

/* End of methods added by Yufeng Zhang*/

void Population::CheckSubtrees(){
	//this function will determine whether the subtree mode should be turned on or off
	//and whether the subtrees should be recalculated
	
	//if subtrees are currently active, see how many trees we have that have accurate subtrees
	//also include any remotes that we have assigned a subtree to, since the next time we 
	//communicate with them we will get a tree that has accurate subtrees		
	
	if(paraMan->subtreeModeActive==true){
		int count=0;
		for(unsigned i=0;i<total_size;i++){
			if(indiv[i].accurateSubtrees==true){
			    if(i<conf->nindivs) count++;
//				if(indiv[i].Fitness() > subMan->currentBest) subMan->currentBest=indiv[i].Fitness();
				}
			}
/*		for(int j=1;j<=paraMan->nremotes;j++)
			if(paraMan->remoteSubtreeAssign[j] > 0) count++;
		if(count==0){  //if nothing is accurate, we obviously need to recalc the subtrees
			paraMan->needUpdate=true;
			}
		if(count>=total_size) subMan->perturb=false;
*/	
		//other conditions for recalcing the subtrees can be put here.
/*		if((bestFitness - paraMan->subtreeDefScore) > paraMan->recalcThresh){
			paraMan->needUpdate=true;
			}
*/		}
	
	if(paraMan->subtreeModeActive==false){
		//determine some conditions for starting/restarting subtree mode here
		//this should probably depend at least partially on the master's score
		//having settled down and model mutations not helping much

		if((int)gen - paraMan->subtreeDefGeneration > paraMan->subtreeInterval && adap->improveOverStoredIntervals < paraMan->subtreeStartThresh){
			paraMan->subtreeModeActive=true;
			paraMan->needUpdate=true;
			}
		}
	
	if(paraMan->subtreeModeActive==true && paraMan->needUpdate==true){
		StartSubtreeMode();
		}

	if(paraMan->subtreeModeActive==true){
		//determine some conditions for stopping subtree mode here
		if((int)gen - paraMan->subtreeDefGeneration > paraMan->subtreeInterval){
			StopSubtreeMode();
			}
		}
	}

void Population::FillPopWithClonesOfBest(){
	Individual *best=&indiv[bestIndiv];
	best->treeStruct->modPart=&best->modPart;
	for(unsigned i=0;i<conf->nindivs;i++){
		if(&indiv[i]!=best){
			indiv[i].treeStruct->RemoveTreeFromAllClas();
			indiv[i].CopySecByRearrangingNodesOfFirst(indiv[i].treeStruct,best);
			indiv[i].mutation_type=-1;
			}
		indiv[i].treeStruct->modPart=&indiv[i].modPart;
		}
	CalcAverageFitness();
	}

void Population::AssignSubtree(int st, int indNum){
	subtreeNode=st;

	//we'll do all of this stuff if we are assigning a new subtree or if 
	//we are assigning 0 (turning off subtree mode)
	for(unsigned i=0;i<conf->nindivs;i++){
	    indiv[i].accurateSubtrees=false;
		newindiv[i].accurateSubtrees=false;
		indiv[i].treeStruct->UnprotectClas();
		}
	ResetMemLevel(dataPart->NTax()-2,claMan->NumClas());
	indiv[bestIndiv].treeStruct->ProtectClas();

	subtreeMemberNodes.clear();
	//add all of the nodenums in the subtree into the subtreeMemberNodes
	//note that the subtree node itself is not added
	if(subtreeNode!=0){
		if(rank==0) assert(indiv[indNum].accurateSubtrees==true);
		indiv[indNum].treeStruct->allNodes[subtreeNode]->left->AddNodesToList(subtreeMemberNodes);
	
		sort(subtreeMemberNodes.begin(),subtreeMemberNodes.end());
		reverse(subtreeMemberNodes.begin(),subtreeMemberNodes.end());
		for(unsigned i=0;i<conf->nindivs;i++){
		    indiv[i].accurateSubtrees=true;
		    newindiv[i].accurateSubtrees=true;
			}
			
		indiv[indNum].SetDirty();
		indiv[indNum].CalcFitness(subtreeNode);

		if(rank!=0){//if we are the master and are going to choose a subtree, don't do this
			indiv[indNum].treeStruct->SetupClasForSubtreeMode(subtreeNode);	

			int nodesNeedingClas=((int)subtreeMemberNodes.size())/2+2;//the nodes in the subtree, plus the subnode itself and it's anc
			ResetMemLevel(nodesNeedingClas,claMan->NumClas());
			}
		}
/*	else{
	//	subMan->active=false;
		for(int i=0;i<conf->nindivs;i++){
		    indiv[i].accurateSubtrees=false;
			newindiv[i].accurateSubtrees=false;
			indiv[i].treeStruct->UnprotectClas();
			}
		ResetMemLevel(params->data->NTax()-2,claMan->NumClas());
		indiv[bestIndiv].treeStruct->ProtectClas();
		}
*/	}

bool Population::SubtreeRecombination(int indivIndex){
	//try recombining the trees from the remotes working on different subtrees
	//in various combinations to give the optimal recombinant

	Individual  currentBest;
	Individual  tempIndiv1;
	bool betterScore=false;

	//DJZ
	while(unusedTrees.size()<2){
		Tree *temp=new Tree();
		unusedTrees.push_back(temp);
		}
	
	tempIndiv1.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();
	currentBest.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();	

	bool *recomEligable=new bool[total_size];

#undef FAKE_PARALLEL
int poo=1;
	
#ifndef FAKE_PARALLEL	
	int count=0;
	for(unsigned i=conf->nindivs;i<total_size;i++){
		if(indiv[i].accurateSubtrees==true && (paraMan->localSubtreeAssign[i-conf->nindivs+1] > 0)){
			paraMan->CheckSubtreeAccuracy(indiv[i].treeStruct);
	       	recomEligable[i]=true;
			if(i>=conf->nindivs)count++;
			}
		else recomEligable[i]=false;
		}
	if(count < (paraMan->nremotes)/2){
	    delete []recomEligable;
	    return false;
		}
#else
	int count=0;
	for(int i=0;i<total_size;i++){
		if(newindiv[i].accurateSubtrees==true && i!=indivIndex){
#ifndef NDEBUG
			paraMan->CheckSubtreeAccuracy(newindiv[i].treeStruct);
#endif
	       	recomEligable[i]=true;
			}
		else recomEligable[i]=false;
		}
#endif

	ofstream subrec("subrec.log", ios::app);
	subrec.precision(10);
	
	subrec << "Subdef " << paraMan->subtreeDefNumber <<  ", " << (int)paraMan->subtrees.size() << " subtrees, defined gen " << paraMan->subtreeDefGeneration << "\n";
	subrec << "Last full recom gen " << paraMan->lastFullRecom <<"\n";
	subrec << "nodenum\tsize\tpriority\tassigned\tbrlen\n";
	int totnode=0;
	for(vector<Subtree*>::iterator it=paraMan->subtrees.begin();it!=paraMan->subtrees.end();it++){
		subrec << (*it)->nodeNum << "\t" << (*it)->taxa << "\t" <<  (*it)->priority << "\t" << (*it)->numAssigned << "\t" << (*it)->blen << "\n";
		totnode+=(*it)->taxa;
		}
	subrec << totnode << " taxa are in a subtree\nremote\tlocalnode\tassignednode\tscore\n";
	for(int p=1;p<paraMan->nremotes+1;p++){
		subrec << p << "\t" << paraMan->localSubtreeAssign[p] << "\t" << paraMan->remoteSubtreeAssign[p] << "\t" << indiv[p+conf->nindivs-1].Fitness() << "\n";
		}

	subrec << "gen\t" << gen << "\nstart\t" << indiv[newindiv[indivIndex].parent].Fitness() << "\n";

	//this is important
	newindiv[indivIndex].CalcFitness(0);
	//we don't want to do this anymore
//	newindiv[indivIndex].treeStruct->ProtectClas();
	
	tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, &newindiv[indivIndex]);
	currentBest.CopySecByRearrangingNodesOfFirst(currentBest.treeStruct, &newindiv[indivIndex]);
	recomEligable[indivIndex]=false;
	
#ifndef FAKE_PARALLEL
	for(unsigned who=conf->nindivs;who<total_size;who++){
#else
	for(int who=0;who<total_size;who++){
#endif
		if(recomEligable[who]==true){

#ifndef FAKE_PARALLEL
			tempIndiv1.treeStruct->SubtreeBasedRecombination(indiv[who].treeStruct, paraMan->localSubtreeAssign[who - conf->nindivs + 1], false, adap->branchOptPrecision);
#else
			tempIndiv1.treeStruct->SubtreeBasedRecombination(newindiv[who].treeStruct, subtreeNode , tempIndiv1.mod->IsModelEqual(newindiv[who].mod), adap->branchOptPrecision);
#endif		

//			OutputFilesForScoreDebugging(&tempIndiv1, poo++);
		//	paupf.flush();
		//	outf.flush();

			tempIndiv1.SetDirty();
			tempIndiv1.CalcFitness(subtreeNode);
			
/*			ofstream poo("debug.log");
			poo.precision(10);
			tempIndiv1.treeStruct->OutputFirstClaAcrossTree(poo, tempIndiv1.treeStruct->root);
			poo.close();
*/			
			subrec << "with " << who << "\t(node " << paraMan->localSubtreeAssign[who - conf->nindivs + 1] << ")\t" << tempIndiv1.Fitness() << "\n";
			if(tempIndiv1.Fitness() > currentBest.Fitness()){
				//if the recombinant we create is better, make it the current best, mark it as 
				//ineligable so we don't try to add it again, and start back at the first eligable
				//recominant
				currentBest.CopySecByRearrangingNodesOfFirst(currentBest.treeStruct, &tempIndiv1, true);
				recomEligable[who]=false;
				//who=conf->nindivs-1;
				}
			else{
			    tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, &currentBest, true);
			    if(tempIndiv1.Fitness()==currentBest.Fitness()) recomEligable[who]=false;
				}
			}
		}
	newindiv[indivIndex].CopySecByRearrangingNodesOfFirst(newindiv[indivIndex].treeStruct, &currentBest, true);

	subrec << "end\t" << currentBest.Fitness() <<  endl;
	subrec.close();
	
	//Return the treestructs that we used temporarily back to the unused tree vector
	tempIndiv1.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(tempIndiv1.treeStruct);
	tempIndiv1.treeStruct=NULL;

	currentBest.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(currentBest.treeStruct);
	currentBest.treeStruct=NULL;
	
	delete []recomEligable;
	paraMan->lastFullRecom=gen;
	return true;
	}


FLOAT_TYPE ParallelManager::ScorePartitioning(int nodeNum, ofstream &pscores){
	
	int size=(int)subtrees.size();
	
	if(size<2 /*|| size>(nremotes-1)*/) return FLT_MAX;

	FLOAT_TYPE blenScore=ZERO_POINT_ZERO, subScore=ZERO_POINT_ZERO, fosterScore=ZERO_POINT_ZERO;
	int fosterTerms=ntax;
	
	int allots[1024];
	
#ifndef MPI_VERSION
nremotes=9;
#endif
	
	int a=0;
	for(vector<Subtree*>::iterator it=subtrees.begin();it!=subtrees.end();it++){
		blenScore -= log((FLOAT_TYPE)(*it)->blen);
		subScore += (*it)->taxa * (*it)->taxa;
		fosterTerms-=(*it)->taxa;
		allots[a++]=(*it)->taxa;
		(*it)->numAssigned=1;
		}
	
	int left=nremotes-size;

	while(left>0){
		int maxnum=0, max=0;
		for(int q=0;q<size;q++){
			if(allots[q]>maxnum){
				maxnum=allots[q];
				max=q;				
				}
			}
		subtrees[max]->numAssigned++;
		allots[max]=subtrees[max]->taxa/subtrees[max]->numAssigned;
		left--;
		}
	int maxallot=0;
	for(int q=0;q<size;q++){
		if(allots[q]>maxallot){
			maxallot=allots[q];
			}	
		}
	
	subScore = sqrt(subScore);
	if(fosterTerms> (ntax/20)) fosterScore = (FLOAT_TYPE)(fosterTerms-(ntax/20))*5;
	else fosterScore=0;
	blenScore*=2.0;
	
	FLOAT_TYPE tot= subScore + blenScore + fosterScore + maxallot;
	
	pscores << nodeNum << "\ts=" << size << "\tscr=" << tot << "\tfost=" << fosterTerms << "\tblenscr=" << blenScore << "\tsubscr=" << subScore << "\tallotscr" << maxallot << "\n";
	for(vector<Subtree*>::iterator it=subtrees.begin();it!=subtrees.end();it++){
		pscores << (*it)->nodeNum << "\t";
		pscores << (*it)->taxa << "\t";
		pscores << (*it)->blen << "\n";
		}
	pscores << "\n";
	return tot;
	}


int ParallelManager::ChooseSubtree(){
  /* subroutine to decide which sub tree to work on then.
     currently we select subtree randomly, and the chance to select a specific subtree is depends on the 
     size of the subtree and how many nodes it is assignged to
  */
  int totalsize = 0;
  int temp = 0;
  int size=(int)subtrees.size();
  if(size<=0) return (0);

  //DJZ see if any subtrees have not been assigned
  bool allassigned=false;
  int unassignedCount=0;

  for(int q=0;q<size;q++) subtrees[q]->numAssigned=0;

  for(int i=1;i<=nremotes;i++){
	if(remoteSubtreeAssign[i]>0){
		int j=0;
		while(subtrees[j]->nodeNum!=remoteSubtreeAssign[i]) j++;
		subtrees[j]->numAssigned++;
		}
	}

  for(int i=0;i<size;i++){
  	if(subtrees[i]->numAssigned==false) unassignedCount++;
	}
  if(unassignedCount==0) allassigned=true;
  
  int nd, max=0;
  for(int i = 0;i<size;i++){
	if(allassigned==true)
		subtrees[i]->priority = subtrees[i]->taxa / subtrees[i]->numAssigned;
	else if(subtrees[i]->numAssigned==0)
	subtrees[i]->priority = subtrees[i]->taxa;
	else subtrees[i]->priority = 0;
	totalsize += subtrees[i]->priority;

	if(subtrees[i]->priority>max){
		nd=i;
		max=subtrees[i]->priority;
		}
	}

  //assign the node with the highest priority overall
  //or the unassigned one with the highest priority
	subtrees[nd]->numAssigned++;
	return subtrees[nd]->nodeNum;

/*//assign a node randomly in proportion to it's size and assignedness
  for(int i=1;i<=size;i++)
    {
      if((randnum>=(temp*1.0/totalsize))&&(randnum<((temp+p[i])*1.0/totalsize))){
		assigned[i]+=1;
		return(node[i]);
		}
      else
	temp += p[i];
    }
*/ 
 //debug_mpi("problem in selectnode...");
  return (subtrees[0]->nodeNum);

}

void ParallelManager::FindNonSubtreeNodes(TreeNode *nd){
	bool subNode=false;
	for(int i=0;i<(int)subtrees.size();i++)
		if(nd->nodeNum==subtrees[i]->nodeNum) subNode=true;
	
	if(nd->nodeNum!=0){
		nonSubtreeNodesforSPR.push_back(nd->nodeNum);
		if(subNode==false && nd->nodeNum>ntax) nonSubtreeNodesforNNI.push_back(nd->nodeNum);
		}
	if(subNode==false){
		if(nd->left) FindNonSubtreeNodes(nd->left);
		if(nd->right) FindNonSubtreeNodes(nd->right);
		if(nd->anc==NULL) FindNonSubtreeNodes(nd->left->next);
		}
	}

void Population::GetRepNums(string &s){
	char buf[100];
	if(conf->bootstrapReps > 0){
		sprintf(buf, "Bootstrap rep %d (of %d) ", currentBootstrapRep, conf->bootstrapReps);
		s += buf;
		}
	if(conf->searchReps > 1){
		sprintf(buf, "Search rep %d (of %d)", currentSearchRep, conf->searchReps);
		s += buf;
		}
	}

void Population::OutputRepNums(ofstream &out){
	if(conf->bootstrapReps > 0 || conf->searchReps > 1){
		if(conf->bootstrapReps > 0) out << "Bootstrap rep " << currentBootstrapRep << " (of " << conf->bootstrapReps << ") ";
		if(conf->searchReps > 1) out << "Search rep " << currentSearchRep << " (of " << conf->searchReps <<  ")";
		out << "\n";
		}
	}

void Population::SetOutputDetails(){
	/*
		DONT_OUTPUT = 0,
		REPLACE = 1,
		APPEND = 2,
		NEWNAME = 4,

		WRITE_CONTINUOUS = 8,
		WRITE_REP_TERM = 16,
		WRITE_REPSET_TERM = 32,
		WRITE_PREMATURE = 64,
		
		FINALIZE_REP_TERM = 128,
		FINALIZE_REPSET_TERM = 256,
		FINALIZE_FULL_TERM = 512,
		FINALIZE_PREMATURE = 1024,

		WARN_PREMATURE = 2048,
		NEWNAME_PER_REP = 4096
	*/
	//not restarted from checkpoint
	//-replace all files that will be output
	if(conf->restart == false){
		screen_output = (output_details) (REPLACE | WRITE_CONTINUOUS | WARN_PREMATURE);
		log_output = (output_details) (REPLACE | WRITE_CONTINUOUS | WARN_PREMATURE);
		if(conf->outputMostlyUselessFiles)
			fate_output = problog_output = swaplog_output = (output_details) (REPLACE | WRITE_CONTINUOUS);
		else
			fate_output = problog_output = swaplog_output = (output_details) (DONT_OUTPUT);

		//not bootstrap
		if(conf->bootstrapReps == 0){
			bootlog_output = (output_details) (DONT_OUTPUT);
			best_output = (output_details) (REPLACE | WRITE_REPSET_TERM | WRITE_PREMATURE | WARN_PREMATURE);
			if(conf->outputCurrentBestTopology) 
				best_output = (output_details) (best_output | WRITE_CONTINUOUS);

			//normal 1 rep
			if(conf->searchReps == 1){
				all_best_output = (output_details) DONT_OUTPUT;
				treelog_output = (output_details) (conf->outputTreelog ? 
					(REPLACE | WRITE_CONTINUOUS | FINALIZE_REP_TERM | FINALIZE_PREMATURE | WARN_PREMATURE) 
					: DONT_OUTPUT);
				}
			//normal multirep run
			else if(conf->searchReps > 1){
				all_best_output = (output_details) (REPLACE | WRITE_REP_TERM | WRITE_PREMATURE | WARN_PREMATURE);
				//treelog_output = (output_details) (conf->outputTreelog ? (REPLACE | WRITE_CONTINUOUS | FINALIZE_REP_TERM | FINALIZE_PREMATURE | WARN_PREMATURE | NEWNAME_PER_REP) : DONT_OUTPUT);
				treelog_output = (output_details) (conf->outputTreelog ? 
					(REPLACE | WRITE_CONTINUOUS | FINALIZE_REP_TERM | FINALIZE_PREMATURE | WARN_PREMATURE | NEWNAME_PER_REP) 
					: DONT_OUTPUT);
				}
			if(conf->scoreOnly)
				log_output = treelog_output = fate_output = problog_output = swaplog_output = (output_details) DONT_OUTPUT;
			}
		//bootstrapping
		else {
			best_output = all_best_output = treelog_output = (output_details) DONT_OUTPUT;
			//bootstrap, 1 OR multiple search reps per bootstrap rep
			//bootlog_output = (output_details) (REPLACE | WRITE_REPSET_TERM | FINALIZE_FULL_TERM | FINALIZE_PREMATURE);
			//WORK - should this include WRITE)PREMATURE?
			bootlog_output = (output_details) (REPLACE | WRITE_REPSET_TERM | FINALIZE_FULL_TERM | FINALIZE_PREMATURE);
			}
		}
	else{//restarted 
		screen_output = (output_details) (APPEND | WRITE_CONTINUOUS | WARN_PREMATURE);
		log_output = (output_details) (APPEND | WRITE_CONTINUOUS | WARN_PREMATURE);
		if(conf->outputMostlyUselessFiles)
			fate_output = problog_output = swaplog_output = (output_details) (APPEND | WRITE_CONTINUOUS);
		else
			fate_output = problog_output = swaplog_output = (output_details) (DONT_OUTPUT);
		//restarted, not bootstrap
		if(conf->bootstrapReps == 0){
			bootlog_output = (output_details) (DONT_OUTPUT);
			best_output = (output_details) (REPLACE | WRITE_REPSET_TERM | WRITE_PREMATURE | WARN_PREMATURE);
			if(conf->outputCurrentBestTopology) 
				best_output = (output_details) (best_output | WRITE_CONTINUOUS);

			//restarted 1 rep search
			if(conf->searchReps == 1){
				all_best_output = (output_details) DONT_OUTPUT;
				treelog_output = (output_details) (conf->outputTreelog ? 
					(APPEND | WRITE_CONTINUOUS | FINALIZE_REP_TERM | FINALIZE_PREMATURE | WARN_PREMATURE) 
					: DONT_OUTPUT);
				}
			//restarted multirep run
			else if(conf->bootstrapReps == 0 && conf->searchReps > 1){
				all_best_output = (output_details) (REPLACE | WRITE_REP_TERM | WARN_PREMATURE);
				treelog_output = (output_details) (conf->outputTreelog ? 
					(APPEND | WRITE_CONTINUOUS | FINALIZE_REP_TERM | FINALIZE_PREMATURE | WARN_PREMATURE | NEWNAME_PER_REP) 
					: DONT_OUTPUT);
				}
			}
		else {//restarted bootstrap, 1 OR multiple search reps per bootstrap rep
			bootlog_output = (output_details) (APPEND | WRITE_REPSET_TERM | FINALIZE_FULL_TERM | FINALIZE_PREMATURE);
			best_output = all_best_output = treelog_output = (output_details) DONT_OUTPUT;
			}
		}
	}

void Population::DetermineFilename(output_details details, char *outname, string suffix){
	char restartString[20];
	char runString[20];

	if(conf->restart && (details & NEWNAME)) sprintf(restartString, ".restart%d", CheckRestartNumber(conf->ofprefix));
	else restartString[0]='\0';
	if(conf->searchReps > 1 && (details & NEWNAME_PER_REP)) sprintf(runString, ".rep%d", currentSearchRep);
	else runString[0]='\0';

	sprintf(outname, "%s%s%s.%s", conf->ofprefix.c_str(), runString, restartString, suffix.c_str());
	}

void Population::InitializeOutputStreams(){
	char temp_buf[100];

	char suffix[100];
	sprintf(suffix, "best");
	DetermineFilename(best_output, temp_buf, suffix);

	//sprintf(temp_buf, "%s%s.best", conf->ofprefix.c_str(), restartString);
	besttreefile = temp_buf;

	if(fate_output != DONT_OUTPUT){
		//initialize the fate file
		if(fate.is_open() == false){
			char suffix[100];
			sprintf(suffix, "fate0%d.log", rank);
			DetermineFilename(fate_output, temp_buf, suffix);
				
			if(fate_output & APPEND)
				fate.open(temp_buf, ios::app);
			else 
				fate.open(temp_buf);
			fate.precision(10);
			}
		#ifdef MPI_VERSION
		fate << "gen\tind\tparent\trecomWith\tscore\tMutType\t#brlen\taccurateSubtrees\tTime\tprecision\n";
		#else
		if(conf->restart) fate << "Restarting from checkpoint...\n";
		OutputRepNums(fate);
		fate << "gen\tind\tparent\tscore\tMutType\t#brlen\tTime\tprecision\n";
		#endif	
		}

	if(problog_output != DONT_OUTPUT){
		//initialize the problog
		if(probLog.is_open() == false){
			char suffix[100];
			sprintf(suffix, "problog0%d.log", rank);
			DetermineFilename(problog_output, temp_buf, suffix);

			if(problog_output & APPEND)
				probLog.open(temp_buf, ios::app);
			else 
				probLog.open(temp_buf);
			if(!probLog.good()) throw ErrorException("problem opening problog");
			}
		if(conf->restart) probLog << "Restarting from checkpoint...\n";
		OutputRepNums(probLog);
		adap->BeginProbLog(probLog, gen);
		}

		//initialize the swaplog
	if(swaplog_output != DONT_OUTPUT){
		if(conf->uniqueSwapBias != 1.0){
			if(swapLog.is_open() == false){
				char suffix[100];
				sprintf(suffix, "swaplog0%d.log", rank);
				DetermineFilename(swaplog_output, temp_buf, suffix);

				if(swaplog_output & APPEND)
					swapLog.open(temp_buf, ios::app);
				else 
					swapLog.open(temp_buf);				
				}
			if(conf->restart) swapLog << "Restarting from checkpoint...\n";
			OutputRepNums(swapLog);
			swapLog << "gen\tuniqueSwaps\ttotalSwaps\n";
			}
		}

	//initialize the log file
	if(log_output != DONT_OUTPUT){
		if(log.is_open() == false){
			char suffix[100];
			sprintf(suffix, "log0%d.log", rank);
			DetermineFilename(log_output, temp_buf, suffix);

			if(log_output & APPEND)
				log.open(temp_buf, ios::app);
			else 
				log.open(temp_buf);		
			log.precision(10);
			}
		OutputRepNums(log);
		if(conf->restart == false)
			log << "random seed = " << rnd.init_seed() << "\n";
		else{
			if(finishedRep == false) 
				if(finishedGenerations == true)
					log << "Restarting run before final optimization " << ", seed " << rnd.init_seed() << ", best lnL " << indiv[bestIndiv].Fitness() << endl;
				else
					log << "Restarting run at generation " << gen << ", seed " << rnd.init_seed() << ", best lnL " << indiv[bestIndiv].Fitness() << endl;
			else
				log << "Restarting from checkpoint...\n";
			}

		log << "gen\tbest_like\ttime\toptPrecision\n";
		}

	//initialize the treelog
	if(treelog_output != DONT_OUTPUT){
		if(treeLog.is_open() == false){
			char suffix[100];
			sprintf(suffix, "treelog0%d.tre", rank);
			DetermineFilename(treelog_output, temp_buf, suffix);

			if(treelog_output & APPEND)
				treeLog.open(temp_buf, ios::app);
			else 
				treeLog.open(temp_buf);		
			treeLog.precision(10);
			}
		treeLog.precision(10);
		if((conf->restart == false && conf->searchReps == 1) ||
			(conf->restart == false && (treelog_output & NEWNAME_PER_REP)) ||
			(conf->restart && (treelog_output & NEWNAME)))
			dataPart->BeginNexusTreesBlock(treeLog);
		AppendTreeToTreeLog(-1, -1);
		}
	
	//initialize the bootstrap tree file
	if(bootlog_output != DONT_OUTPUT){
		if(rank==0 && bootLog.is_open() == false){
			char suffix[100];
			sprintf(suffix, "boot.tre");
			DetermineFilename(bootlog_output, temp_buf, suffix);
			
			if(bootlog_output & APPEND){
				if(FileExists(temp_buf) && FileIsNexus(temp_buf)){
					//this will verify whether we previously started a trees block in the bootstrap file
					//if so, we shouldn't do so again.  If not, we need to start it now
					bootLog.open(temp_buf, ios::app);
					}
				else{
					bootLog.open(temp_buf, ios::app);
					dataPart->BeginNexusTreesBlock(bootLog);
					}
				}
			else{
				bootLog.open(temp_buf);
				dataPart->BeginNexusTreesBlock(bootLog);
				}
			bootLog.precision(10);
			}
		if(conf->outputPhylipTree){
			if(!(bootLogPhylip.is_open())){
				char suffix[100];
				sprintf(suffix, "boot.phy");
				DetermineFilename(bootlog_output, temp_buf, suffix);

				if(bootlog_output & APPEND)
					bootLogPhylip.open(temp_buf, ios::app);
				else 
					bootLogPhylip.open(temp_buf);		
				bootLogPhylip.precision(10);
				}	
			}
		}

	ClearDebugLogs();
	
	#ifdef DEBUG_SCORES
	outf.open("toscore.tre");
	paupf.open("toscore.nex");
	#endif

	#ifdef VARIABLE_OPTIMIZATION
	outf.open("toscore.tre");
	paupf.open("toscore.nex");
	#endif

	#ifdef PERIODIC_SCORE_DEBUG
	outf.open("toscore.tre");
	paupf.open("toscore.nex");
	#endif	
	}

/* OLD WAY
void Population::InitializeOutputStreams(){
	char temp_buf[100];
	char restart[12];

	if(conf->restart == true){
		//check if this run has been restarted before.  If so, increment the restart number
		sprintf(restart, ".restart%d", CheckRestartNumber(conf->ofprefix));
		}
	else restart[0]='\0';

	sprintf(temp_buf, "%s%s.best.tre", conf->ofprefix.c_str(), restart);
	besttreefile = temp_buf;

	if(conf->outputMostlyUselessFiles){
		//initialize the fate file
		if (rank > 9)
			sprintf(temp_buf, "%s%s.fate%d.log", conf->ofprefix.c_str(), restart, rank);
		else
			sprintf(temp_buf, "%s%s.fate0%d.log", conf->ofprefix.c_str(), restart, rank);
		
		fate.open(temp_buf);
		fate.precision(10);
		#ifdef MPI_VERSION
		fate << "gen\tind\tparent\trecomWith\tscore\tMutType\t#brlen\taccurateSubtrees\tTime\tprecision\n";
		#else
		fate << "gen\tind\tparent\tscore\tMutType\t#brlen\tTime\tprecision\n";
		#endif	

		//initialize the problog
		if (rank > 9)
			sprintf(temp_buf, "%s%s.problog%d.log", conf->ofprefix.c_str(), restart, rank);
		else
			sprintf(temp_buf, "%s%s.problog0%d.log", conf->ofprefix.c_str(), restart, rank);
		probLog.open(temp_buf);
		if(!probLog.good()) throw ErrorException("problem opening problog");
		adap->BeginProbLog(probLog, gen);

		//initialize the swaplog
		if(conf->uniqueSwapBias != ONE_POINT_ZERO){
			sprintf(temp_buf, "%s%s.swap.log", conf->ofprefix.c_str(), restart);
			swapLog.open(temp_buf);
			swapLog << "gen\tuniqueSwaps\ttotalSwaps\n";
			}
		}

	//initialize the log file
	if (rank > 9)
		sprintf(temp_buf, "%s%s.log%d.log", conf->ofprefix.c_str(), restart, rank);
	else
		sprintf(temp_buf, "%s%s.log0%d.log", conf->ofprefix.c_str(), restart, rank);
	log.open(temp_buf);
	log.precision(10);
	if(conf->restart == false)
		log << "random seed = " << rnd.init_seed() << "\n";
	else log << "Restarting run at generation " << gen << ", seed " << rnd.init_seed() << ", best lnL " << BestFitness() << endl;

	log << "gen\tbest_like\ttime\toptPrecision\n";

	//initialize the treelog
	if(conf->outputTreelog){
		if (rank > 9)
			sprintf(temp_buf, "%s%s.treelog%d.tre", conf->ofprefix.c_str(), restart, rank);
		else
			sprintf(temp_buf, "%s%s.treelog0%d.tre", conf->ofprefix.c_str(), restart, rank);

		treeLog.open(temp_buf);
		treeLog.precision(10);

		data->BeginNexusTreesBlock(treeLog);
		}
	
	//initialize the bootstrap tree file
	if(conf->bootstrapReps > 0 && rank==0){
		sprintf(temp_buf, "%s%s.boot.tre", conf->ofprefix.c_str(), restart);

		bootLog.open(temp_buf);
		bootLog.precision(10);

		data->BeginNexusTreesBlock(bootLog);

		if(conf->outputPhylipTree == true){
			sprintf(temp_buf, "%s%s.boot.phy", conf->ofprefix.c_str(), restart);

			bootLogPhylip.open(temp_buf);
			bootLogPhylip.precision(10);
			}	
		}	

	ClearDebugLogs();
	
	#ifdef DEBUG_SCORES
	outf.open("toscore.tre");
	paupf.open("toscore.nex");
	#endif

	#ifdef VARIABLE_OPTIMIZATION
	outf.open("toscore.tre");
	paupf.open("toscore.nex");
	#endif

	#ifdef PERIODIC_SCORE_DEBUG
	outf.open("toscore.tre");
	paupf.open("toscore.nex");
	#endif	
	
	}
*/

void Population::SetNewBestIndiv(int indivIndex){
	//this should be called when a new best individual is set outside of
	//CalcAverageFitness.  Particularly important for parallel.
	bestIndiv=indivIndex;
	globalBest=bestFitness=prevBestFitness=BestFitness();
	for(unsigned i=0;i<total_size;i++){
		if(i != bestIndiv){
			indiv[i].treeStruct->UnprotectClas();
			}
		}
	indiv[bestIndiv].treeStruct->ProtectClas();
	}

void Population::FinalizeOutputStreams(int type){
	/*the types are:
		repTerm = 0
		repsetTerm = 1
		fullTerm = 2
	*/

	//things should already have been finalized at the end of the previous execution
	if(restartedAfterTermination)
		return;

	bool prematureTermination = conf->checkpoint ? genTermination : (genTermination | timeTermination | userTermination);

	if(prematureTermination && type == 0){
		if(log_output & WARN_PREMATURE)
			log << TerminationWarningMessage().c_str() << endl;
		if(treelog_output & WARN_PREMATURE)
			if(treeLog.is_open()) 
				treeLog << TerminationWarningMessage().c_str() << endl;
		}

	//in these cases the termination is essentially a pause, so don't finalize or write anything
	if(workPhaseTermination || (conf->checkpoint && (timeTermination || userTermination)))
		return;

	bool repTerm, repsetTerm, fullTerm;
	if(prematureTermination){
		fullTerm = false;
		repsetTerm = false;
		repTerm = false;		
		}
	else if(type == 2){
		fullTerm = true;
		repsetTerm = false;
		repTerm = false;
		}
	else if(type == 1) {
		fullTerm = false;
		repsetTerm = true;
		repTerm = false;
		}
	else {
		fullTerm = false;
		repsetTerm = false;
		repTerm = true;
		}

	//if(((conf->bootstrapReps == 0 || currentBootstrapRep == conf->bootstrapReps) && (currentSearchRep == conf->searchReps)) || prematureTermination == true){
	if(log.is_open()){
		if((prematureTermination && (log_output & FINALIZE_PREMATURE)) ||
			(!prematureTermination &&
			   ((repTerm && (log_output & FINALIZE_REP_TERM)) || (repsetTerm && (log_output & FINALIZE_REPSET_TERM)) || (fullTerm && (log_output & FINALIZE_FULL_TERM))))
			) 
			log.close();
		}

	if(fate.is_open()){
		if((prematureTermination && (fate_output & FINALIZE_PREMATURE)) ||
			(!prematureTermination && 
			   ((repTerm && (fate_output & FINALIZE_REP_TERM)) || (repsetTerm && (fate_output & FINALIZE_REPSET_TERM)) || (fullTerm && (fate_output & FINALIZE_FULL_TERM))))
			) 
			fate.close();
		}

	if(probLog.is_open()){
		if((prematureTermination && (problog_output & FINALIZE_PREMATURE)) ||
			(!prematureTermination && 
				((repTerm && (problog_output & FINALIZE_REP_TERM)) || (repsetTerm && (problog_output & FINALIZE_REPSET_TERM)) || (fullTerm && (problog_output & FINALIZE_FULL_TERM))))
			)
			probLog.close();
		}

	if(swapLog.is_open()){
		if((prematureTermination && (swaplog_output & FINALIZE_PREMATURE)) || 
			(!prematureTermination &&
				((repTerm && (swaplog_output & FINALIZE_REP_TERM)) || (repsetTerm && (swaplog_output & FINALIZE_REPSET_TERM)) || (fullTerm && (swaplog_output & FINALIZE_FULL_TERM))))
			)
			swapLog.close();
		}

	if(treeLog.is_open()){
		if((prematureTermination && (treelog_output & FINALIZE_PREMATURE)) ||
			(!prematureTermination &&
				((repTerm && (treelog_output & FINALIZE_REP_TERM)) || (repsetTerm && (treelog_output & FINALIZE_REPSET_TERM)) || (fullTerm && (treelog_output & FINALIZE_FULL_TERM))))
			){
			treeLog << "end;\n";
			treeLog.close();
			}
		}

	if(bootLog.is_open()){
		if((prematureTermination && (bootlog_output & FINALIZE_PREMATURE)) ||
			(!prematureTermination && 
			   ((repTerm && (bootlog_output & FINALIZE_REP_TERM)) || (repsetTerm && (bootlog_output & FINALIZE_REPSET_TERM)) || (fullTerm && (bootlog_output & FINALIZE_FULL_TERM))))
			){
			bootLog << "end;\n";
			bootLog.close();
			}
		}

	if(bootLogPhylip.is_open()){
		if((prematureTermination && (bootlog_output & FINALIZE_PREMATURE)) ||
			(!prematureTermination && 
			   ((repTerm && (bootlog_output & FINALIZE_REP_TERM)) || (repsetTerm && (bootlog_output & FINALIZE_REPSET_TERM)) || (fullTerm && (bootlog_output & FINALIZE_FULL_TERM))))
			)
			bootLogPhylip.close();
		}

	#ifdef DEBUG_SCORES
	outf << "end;\n";
	outf.close();
	paupf << "end;\n";
	paupf.close();
	#endif
	}

/* OLD WAY
void Population::FinalizeOutputStreams(){
	if(prematureTermination == true){
		log << "***NOTE: Run was terminated before termination condition was reached!\nLikelihood scores, topologies and model estimates obtained may not be fully optimal!***" << endl;
		if(treeLog.is_open()) treeLog << "[! ***NOTE: GARLI run was terminated before termination condition was reached!\nLikelihood scores, topologies and model estimates obtained may not be fully optimal!***]" << endl;
		}
	fate.close();
	log.close();
	if(treeLog.is_open()){
		AppendTreeToTreeLog(-1);
		treeLog << "end;\n";
		treeLog.close();
		}
	if(bootLog.is_open()){
		bootLog << "end;\n";
		bootLog.close();
		}
	if(bootLogPhylip.is_open()) bootLogPhylip.close();
	probLog.close();
	if(swapLog.is_open()) swapLog.close();
	
	#ifdef DEBUG_SCORES
	outf << "end;\n";
	outf.close();
	paupf << "end;\n";
	paupf.close();
	#endif
	}
*/

/*BMERGE this needs to be updated
void Population::FindLostClas(){
	vector<CondLikeArraySet *> arr;
	
	for(unsigned i=0;i<total_size;i++){
		Tree *t=indiv[i].treeStruct;
		if(! (claMan->IsDirty(t->allNodes[0]->claIndexDown)))
			arr.push_back(claMan->GetClaSet(t->allNodes[0]->claIndexDown));
		if(! (claMan->IsDirty(t->allNodes[0]->claIndexUL)))
			arr.push_back(claMan->GetClaSet(t->allNodes[0]->claIndexUL));			
		if(! (claMan->IsDirty(t->allNodes[0]->claIndexUR)))
			arr.push_back(claMan->GetClaSet(t->allNodes[0]->claIndexUR));
		for(int n=t->getNumTipsTotal()+1;n<t->getNumNodesTotal();n++){
			if(! (claMan->IsDirty(t->allNodes[n]->claIndexDown)))
				arr.push_back(claMan->GetClaSet(t->allNodes[n]->claIndexDown));
			if(! (claMan->IsDirty(t->allNodes[n]->claIndexUL)))
				arr.push_back(claMan->GetClaSet(t->allNodes[n]->claIndexUL));			
			if(! (claMan->IsDirty(t->allNodes[n]->claIndexUR)))
				arr.push_back(claMan->GetClaSet(t->allNodes[n]->claIndexUR));
			}
		}
	sort(arr.begin(), arr.end());
	for(vector<CondLikeArraySet*>::iterator vit=arr.begin();vit!=arr.end();){
		if(*(vit)==*(vit+1)) arr.erase(vit+1);
		else vit++;
		if((vit+1)==arr.end()){
			break;
			}
		}
	assert(arr.size() + claMan->NumFreeClas() == claMan->NumClas());
	}
*/

void Population::LogNewBestFromRemote(FLOAT_TYPE scorediff, int ind){
#ifdef MPI_VERSION
        adap->bestFromRemoteNum[0]++;
        adap->bestFromRemote[0]+=scorediff;
        
        indiv[0].treeStruct->CalcBipartitions(true);
        indiv[ind].treeStruct->CalcBipartitions(true);
        
        //check if this is a new topology
        if(indiv[0].treeStruct->IdenticalTopology(indiv[ind].treeStruct->root) == false){
			debug_mpi("\ttopo different from prev best");
			AppendTreeToTreeLog(-1, ind);
			if(scorediff > significantTopoChange) lastTopoImprove=gen;
			}
		else{
			debug_mpi("\ttopo same as prev best");
//			AppendTreeToTreeLog(0, ind);
			}
        
#endif
        }

void Population::CheckRemoteReplaceThresh(){
#ifdef MPI_VERSION
	if(gen < adap->intervalLength * adap->intervalsToStore) return;
	FLOAT_TYPE totBestFromRemote=ZERO_POINT_ZERO;
	int totBestFromRemoteNum=0;
	for(int i=0;i<adap->intervalsToStore;i++){	
	         totBestFromRemoteNum += adap->bestFromRemoteNum[i];
        	 totBestFromRemote += adap->bestFromRemote[i];
		}
	if(totBestFromRemoteNum==0) paraMan->ReduceUpdateThresh();

#endif
}

/* 7/21/06 needs to be updated
void Population::SPRoptimization(int indivIndex){
//	for(int i=0;i<conf->nindivs;i++)
//		indiv[i].ResetIndiv();

//	bool topoChange=false;

	for(int reps=0;reps<1;reps++){
		int cutnum = newindiv[indivIndex].treeStruct->GetRandomNonRootNode();
		SPRoptimization(indivIndex, adap->limSPRrange, cutnum);
		}
	
/*	if(topoChange==true){
		if(topologies[indiv[bestIndiv].topo]->nInds>1){
			topologies[indiv[bestIndiv].topo]->RemoveInd(bestIndiv);
			indiv[bestIndiv].topo=ntopos++;
			topologies[indiv[bestIndiv].topo]->AddInd(bestIndiv);
			assert(topologies[indiv[bestIndiv].topo]->nInds==1);
			}
		topologies[indiv[bestIndiv].topo]->gensAlive=0;
		TopologyList::ntoposexamined++;
		UpdateTopologyList(indiv);
		}	
*/	
//	CalcAverageFitness();
//	OutputFilesForScoreDebugging();
       
//}

/* 7/21/06 needs to be updated
bool Population::SPRoptimization(int indivIndex, int range, int cutnum ){
	//DJZ 1/23/04 this is based on the NNI optimization function by Alan (ie, exhaustive nnis)
	//the main difference is that only one node will be used as the node to be cut off and
	//reattached, but then all reattachment points within a radius will be tried.
	//the marking of nodes as dirty is also necessarily different
	subset sprRange;
	
	Individual  currentBest;
	Individual  tempIndiv1;
	FLOAT_TYPE bestSPRFitness; 
	bool topoChange=false;
	
	ofstream outf("sprdebug.tre");
	ofstream scr("sprscores.tre");
	scr.precision(10);

	//DJZ
	while(unusedTrees.size()<2){
		Tree *temp=new Tree();
		unusedTrees.push_back(temp);
		}

	//bestSPRFitness = indiv[newindiv[indivIndex].parent].Fitness();
	bestSPRFitness = -1e100;


	tempIndiv1.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();
	currentBest.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();	
	//

	tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, &newindiv[indivIndex]);
	currentBest.CopySecByRearrangingNodesOfFirst(currentBest.treeStruct, &newindiv[indivIndex]);

	TreeNode **thenodes=tempIndiv1.treeStruct->allNodes;

	//choose the nodenum to be cut
//	int cutnum = params->rnd.random_int(tempIndiv1.treeStruct->numNodesTotal -1 ) +1;
	TreeNode *cutnode= thenodes[cutnum];
	
	//now determine the nodes that fall within the reattachment radius
	//this is Alan's code for putting together the subset, with a bit of my alteration
	//the subset will be centered on cutnode's anc, AKA connector
	sprRange.setseed(cutnode->anc->nodeNum);
	int connector=cutnode->anc->nodeNum;
	
    for(int i = 0;i<range;i++){
		int j = sprRange.total;
		for(int k=0; k < j; k++){
			if(sprRange.front[k]==i){
				TreeNode *cur=thenodes[sprRange.element[k]];
				if(cur->left!=NULL) 
						sprRange.addelement(cur->left->nodeNum, i+1, sprRange.pathlength[k]+cur->left->dlen);
				if(cur->right!=NULL)
						sprRange.addelement(cur->right->nodeNum, i+1, sprRange.pathlength[k]+cur->right->dlen);
				if(cur->anc!=NULL) 
						sprRange.addelement(cur->anc->nodeNum, i+1, sprRange.pathlength[k]+cur->dlen);
				}// end of loop through element of current subset
		    }// end of loop to findrange
		}

	if(cutnode->next != NULL) sprRange.elementremove(cutnode->next->nodeNum);
	if(cutnode->prev != NULL) sprRange.elementremove(cutnode->prev->nodeNum);
	sprRange.elementremove(connector); //connecting to the sib recreates the original tree
	//remove the nodes that are actually part of the subtree being cut, starting with connector
	thenodes[cutnum]->RemoveSubTreeFromSubset(sprRange, true);

	sprRange.compact();

	int broken=0;
	while(sprRange.element[broken]!=0){

		tempIndiv1.treeStruct->SPRMutate(cutnum, sprRange.element[broken++], adap->branchOptPrecision, 0, 0);
		if(sprRange.element[broken]==0 && sprRange.element[broken+1]!=0) broken++;

		//indiv[indivIndex].treeStruct->SetAllTempClasDirty();
		//Because a large section of the tree will be shared between the different attachment
		//points, we should see a decent savings by only making the temp clas dirty that we 
		//know might change, which should only be those that are considered as reattachments.
		//newindiv[indivIndex].treeStruct->SetSpecifiedTempClasDirty(sprRange.element);
				
		tempIndiv1.SetDirty();
		
		tempIndiv1.CalcFitness(0);
		
		//debug the scoring of the spr trees
/*		outf << "  utree " << gen << sprRange.element[broken] << "= ";
		tempIndiv1.treeStruct->root->MakeNewick(treeString);
		outf << treeString << ";" << endl;

		scr << tempIndiv1.Fitness() << endl;
		//
*/		
//		if(tempIndiv1.Fitness() > (bestSPRFitness + 0.01))
/*		if(tempIndiv1.Fitness() > bestSPRFitness)
			{
			bestSPRFitness = tempIndiv1.Fitness();
			currentBest.CopySecByRearrangingNodesOfFirst(currentBest.treeStruct, &tempIndiv1, true);
			topoChange=true;
			}
	
		//make the tempIndiv equal to the starting tree
		tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, &newindiv[indivIndex], true);
		} //end of loop through all possible NNIs 
		
	if(topoChange==true){
		newindiv[indivIndex].CopySecByRearrangingNodesOfFirst(newindiv[indivIndex].treeStruct, &currentBest, true);
		}
	
	//Return the treestructs that we used temporarily back to the unused tree vector
	tempIndiv1.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(tempIndiv1.treeStruct);
	tempIndiv1.treeStruct=NULL;
	currentBest.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(currentBest.treeStruct);
	currentBest.treeStruct=NULL;	
	
//	newindiv[indivIndex].treeStruct->SetAllTempClasDirty();
	newindiv[indivIndex].mutation_type |= Individual::exlimSPR;
	return topoChange;
}
*/
/* 7/21/06 needs to be update
void Population::SPRPerturbation(int sourceInd, int indivIndex){
	assert(0);
	//7/21/06 needs to be fixed to deal with changes made in
	//constraint implementation
	
	/*
	Individual  currentBest;
	Individual  tempIndiv1;
	Individual *source=&indiv[sourceInd];
	int range=pertMan->sprPertRange;
	FLOAT_TYPE thresh=10000.0;

	
//	ofstream outf("sprdebug.tre");
//	ofstream scr("sprscores.tre");
//	scr.precision(10);

	//DJZ
	while(unusedTrees.size()<2){
		Tree *temp=new Tree();
		unusedTrees.push_back(temp);
		}

	tempIndiv1.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();
	currentBest.treeStruct=*(unusedTrees.end()-1);
	unusedTrees.pop_back();	
	

	for(int cycle=0;cycle < pertMan->numSprCycles;cycle++){
		FLOAT_TYPE previousFitness=source->Fitness();
		FLOAT_TYPE bestDiff=-thresh;

		if(cycle==0){
			tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, source);
			currentBest.CopySecByRearrangingNodesOfFirst(currentBest.treeStruct, source);
			}
		else{
			tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, source, true);
			currentBest.CopySecByRearrangingNodesOfFirst(currentBest.treeStruct, source, true);			
			}

		TreeNode **thenodes=tempIndiv1.treeStruct->allNodes;

		//this is a little odd, but just call a normal SPRMutate with the proper range
		//so that the possible reattachment points within that range are gathered properly
		int cutnum=tempIndiv1.treeStruct->SPRMutate(range, adap->branchOptPrecision);

		tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, source, true);

		char filename[50];
		sprintf(filename, "pertreport%d.log", rank);
		ofstream pert(filename, ios::app);
		pert.precision(10);
	
		subset *sprRange=&(tempIndiv1.treeStruct->sprRange);
		
		pert.precision(10);
	//	pert2.precision(10);
		pert << "gen " << gen << " start " << source->Fitness() << "\t" << sprRange->total << " possible attachments\n";
		
		int bestDist=0;
		int broken=sprRange->total-1;
		while(broken>=0){
			if(! (broken % (int)ceil((FLOAT_TYPE)sprRange->total/5))) outman.UserMessageNoCR(".");
			outman.flush();

			tempIndiv1.treeStruct->SPRMutate(cutnum, sprRange->element[broken--], adap->branchOptPrecision, 0, 0);
			if(sprRange->element[broken]==0 && broken>0) broken--;

			tempIndiv1.SetFitness(tempIndiv1.treeStruct->lnL);;

			//divide the score difference by the square root of the node distance, to favor longer moves
			FLOAT_TYPE diff=(tempIndiv1.Fitness()-previousFitness)/sqrt(sprRange->front[broken]+1.0);
			if(diff>0) diff=(tempIndiv1.Fitness()-previousFitness)*(sprRange->front[broken]+1);
	//		pert2 << "node=\t" << sprRange->element[broken] << "\tdist=\t" << sprRange->front[broken] << "\tscore=\t" << tempIndiv1.Fitness() << "\t" << diff << "\t" << tempIndiv1.Fitness()-previousFitness << "\n";
			if(diff > bestDiff){
				bestDiff=diff;
				currentBest.CopySecByRearrangingNodesOfFirst(currentBest.treeStruct, &tempIndiv1, true);
				bestDist=sprRange->front[broken];
				pert << diff << "\t" << bestDist << "\n";
				}
			tempIndiv1.CopySecByRearrangingNodesOfFirst(tempIndiv1.treeStruct, source, true);
			} 
					
		if(bestDiff>-thresh){
			indiv[indivIndex].CopySecByRearrangingNodesOfFirst(indiv[indivIndex].treeStruct, &currentBest, true);
			}
		//set this tree up as the source for the next cycle
		source=&indiv[indivIndex];

		indiv[indivIndex].mutation_type |= Individual::exlimSPR;
		pert << "end score=" << currentBest.Fitness() << endl;
		
		outman.UserMessage("Accepted SPR with range of %d.  Current score= %.4f", bestDist, indiv[indivIndex].Fitness());
		}

	
	//Return the treestructs that we used temporarily back to the unused tree vector
	tempIndiv1.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(tempIndiv1.treeStruct);
	tempIndiv1.treeStruct=NULL;
	currentBest.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(currentBest.treeStruct);
	currentBest.treeStruct=NULL;
	}
*/

/*
void Population::CheckPerturbParallel(){
	if(paraMan->perturbModeActive==true){
		if(paraMan->allSent == false){
			for(int i=1;i<=paraMan->nremotes;i++){
				if(paraMan->needToSend[i]==true) break;
				if(i==paraMan->nremotes) paraMan->allSent=true;
				}		
			if(paraMan->allSent==true){
				//the pert generation is recorded as when we sent our last message
				pertMan->lastPertGeneration=gen;
				}
			}
		//keep the perturbModeActive flag true for a while so the master doesn't replace the remotes perturbed trees
		if(gen - pertMan->lastPertGeneration > pertMan->minPertInterval){
			paraMan->perturbModeActive=false;
//			pertMan->lastPertGeneration=gen;
			}
		}
	else if(pertMan->pertAbandoned==false && (gen - pertMan->lastPertGeneration) >= pertMan->minPertInterval
		&& (adap->improveOverStoredIntervals < pertMan->pertThresh)/* && (adap->branchOptPrecision == adap->minOptPrecision)*/ /*){
		for(int i=1;i<=paraMan->nremotes;i++){
			paraMan->needToSend[i]=true;
			}
		paraMan->perturbModeActive=true;
		paraMan->allSent=false;
		pertMan->lastPertGeneration=gen;
		}	
	}
*/

/* 7/21/06 needs to be updated
void Population::CheckPerturbSerial(){

	if(pertMan->pertType < 3 ){
	 	if(pertMan->pertType==1 && (gen - pertMan->lastPertGeneration) >= pertMan->minPertInterval/2 
	 		&& adap->randNNIweight != adap->origRandNNIweight){
			adap->randNNIweight=adap->origRandNNIweight;
//			pertMan->lastPertGeneration=gen;
			}


		if(pertMan->pertAbandoned==false && (gen - pertMan->lastPertGeneration) >= pertMan->minPertInterval 
			&& (adap->improveOverStoredIntervals < pertMan->pertThresh) /*&& (adap->branchOptPrecision == adap->minOptPrecision)*/ /*){
			if(pertMan->numPertsNoImprove <= pertMan->maxPertsNoImprove){
				if(BestFitness() > bestSinceRestart.Fitness()){
					StoreBestForPert();
					pertMan->numPertsNoImprove=0;
					}
				else{
					//if we haven't done better than the best we had before the previous perturbation, restore to 
					//that point and perturb again
					pertMan->numPertsNoImprove++;
					RestoreBestForPert();
					//abandoning perturbations
					if(pertMan->numPertsNoImprove > pertMan->maxPertsNoImprove){
						pertMan->pertAbandoned=true;
						pertMan->lastPertGeneration=gen;
						return;
						}
					}
												
				if(pertMan->pertType==1){
					int indToReplace = (bestIndiv==0 ? 1 : 0);
					NNIPerturbation(bestIndiv, indToReplace);
					SetNewBestIndiv(indToReplace);
					FillPopWithClonesOfBest();
					AppendTreeToTreeLog(-1, bestIndiv);
					//try disallowing NNIs immediately after the perturbation
					adap->randNNIweight=0.0;
					adap->randNNIprob=0.0;
					}
				else if(pertMan->pertType==0){
					//branch length perturbation
					int indToReplace = (bestIndiv==0 ? 1 : 0);
					indiv[indToReplace].treeStruct->PerturbAllBranches();
					indiv[indToReplace].SetDirty();
					indiv[indToReplace].CalcFitness(0);
					SetNewBestIndiv(indToReplace);
					FillPopWithClonesOfBest();
					}
				else{
					FLOAT_TYPE startscore=BestFitness();
					FLOAT_TYPE curscore;
//					do{
					int indToReplace = (bestIndiv==0 ? 1 : 0);
					int source=bestIndiv;
					outman.UserMessage("Performing SPR Perturbation.  Starting score=%.4f", BestFitness());
					SPRPerturbation(source, indToReplace);
					indiv[indToReplace].CalcFitness(0);
					SetNewBestIndiv(indToReplace);
					FillPopWithClonesOfBest();
					curscore=indiv[indToReplace].Fitness();
					AppendTreeToTreeLog(-1, bestIndiv);
					}
				pertMan->lastPertGeneration=gen;
				adap->reset=true;
				gen++;
				OutputFate();
				}
			}
		}

	else if(pertMan->pertType==3){
		if(pertMan->ratcheted==false){
			if(pertMan->pertAbandoned==false && (gen - pertMan->lastPertGeneration) >= pertMan->minPertInterval && adap->improveOverStoredIntervals < pertMan->pertThresh){
				if(BestFitness() > bestSinceRestart.Fitness()){
					StoreBestForPert();
					pertMan->numPertsNoImprove=0;
					}
				else{
					//if we haven't done better than the best we had before the previous perturbation, restore to 
					//that point and reweight again
					RestoreBestForPert();
					pertMan->numPertsNoImprove++;

					//abandoning perturbations
					if(pertMan->numPertsNoImprove > pertMan->maxPertsNoImprove){
						pertMan->pertAbandoned=true;
						return;
						}
					}
				pertMan->ratcheted=true;
				params->data->ReserveOriginalCounts();
				params->data->Reweight(pertMan->ratchetProportion);
			
				claMan->MakeAllHoldersDirty();
				for(int i=0;i<total_size;i++) indiv[i].SetDirty();
				CalcAverageFitness();
				bestFitness=BestFitness();
				pertMan->lastPertGeneration=gen;
				pertMan->scoreAfterRatchet=BestFitness();
				adap->reset=true;
				gen++;
				OutputFate();
				outman.UserMessage("Performing ratcheting: reweighting %.1f percent of characters.", pertMan->ratchetProportion*100);
				char filename[50];
				if(rank < 10)
					sprintf(filename, "pertreport0%d.log", rank);
				else 
					sprintf(filename, "pertreport%d.log", rank);
				ofstream pert(filename, ios::app);
				pert << "Performing ratcheting: reweighting " << pertMan->ratchetProportion*100 << " percent of characters." << endl; 
				pert.close();
				}
			}

		//turn ratchet off
		else{
			if((gen - pertMan->lastPertGeneration) >= pertMan->ratchetMaxGen || BestFitness() - pertMan->scoreAfterRatchet > pertMan->ratchetOffThresh){
				TurnOffRatchet();
				gen++;
				OutputFate();
				}
			}
		}
	}
*/

//SINGLE SITE FUNCTIONS
void Population::OptimizeSiteRates(){
	SeedPopulationWithStartingTree(1);

	const SequenceData *data = dataPart->GetSubset(0);

	//store a backup of the exisiting tree and blens
	Individual tempIndiv;
	tempIndiv.treeStruct=new Tree();
	tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, &indiv[0]);	

	char filename[100];
	sprintf(filename, "%s.siterates.log", conf->ofprefix.c_str());
	ofstream out(filename);

	Tree::min_brlen = 1.0e-20;
	Tree::max_brlen = 1.0e10;

	const int lastConst=data->LastConstant();

	typedef pair<FLOAT_TYPE, FLOAT_TYPE> float_pair;
	float_pair rateAndScore;
	vector<float_pair> allRates;

	for(int i=0;i<data->NChar();i++){
		if(i <= lastConst) 
			rateAndScore = pair<FLOAT_TYPE, FLOAT_TYPE>(ZERO_POINT_ZERO, indiv[0].modPart.GetModel(0)->StateFreq((data->GetConstStates())[i]));
		else{
			indiv[0].treeStruct->MakeAllNodesDirty();
			Tree::siteToScore = i;
			rateAndScore = indiv[0].treeStruct->OptimizeSingleSiteTreeScale(adap->branchOptPrecision);
			//restore the original blens
			indiv[0].CopySecByRearrangingNodesOfFirst(indiv[0].treeStruct, &tempIndiv, true);
			}
		allRates.push_back(rateAndScore);
		}
	out << "site#\tsiteRate\tsitelnL" << endl;
	for(int i=0;i<data->GapsIncludedNChar();i++){
		int packedColumn = data->Number(i);
		assert(packedColumn < (int)allRates.size());
		out << i+1 << "\t";
		if(packedColumn == -1) out << "NA\tNA" << endl;
		else out << allRates[packedColumn].first << "\t" << allRates[packedColumn].second << endl;
		}
	out.close();
	outman.UserMessage("Site-rate estimation complete.");
	}
