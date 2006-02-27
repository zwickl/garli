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

#include "memchk.h"
#ifdef MPI_VERSION
	#include <mpi.h>
#endif
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <cmath>
#include <iomanip>
#ifdef WINDOWS
#include <conio.h>
#endif

using namespace std;

#include <signal.h>
#include "population.h"
#include "individual.h"
#include "parameters.h"
#include "translatetable.h"
#include "mlhky.h"
#include "tree.h"
#include "topologylist.h"
#include "funcs.h"
#include "clamanager.h"
#include "stopwatch.h"
#include "bipartition.h"
#include "subset.h"
#include "adaptation.h"
#include "errorexception.h"

extern char programName[81];

int memLevel;
int calcCount=0;
int optCalcs;


ofstream outf, paupf;
int tempGlobal=1;

double globalBest;

#undef PERIODIC_SCORE_DEBUG

#undef DEBUG_SCORES

#undef NNI_SPECTRUM

#undef MASTER_DOES_SUBTREE

#undef VARIABLE_OPTIMIZATION

int debug_mpi(const char* fmt, ...);
int QuitNow();
void InterruptMessage( int );
void ClearDebugLogs();

#define EPSILON          1e-8

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

void InterruptMessage( int )
{
	askQuitNow = 1;
}

void Population::CatchInterrupt()
{
	if( signal( SIGINT, SIG_IGN ) != SIG_IGN ){
		signal( SIGINT, InterruptMessage );
		}
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
*/	#endif
	}

Population::~Population()
{
	EliminateDuplicateTreeReferences();  // TODO this be broken
	
	for (int i = 0; i < total_size; ++i)	{
		for (int j = 0; j < total_size; ++j)	{
			if (newindiv[i].treeStruct == indiv[j].treeStruct)	{
				newindiv[i].treeStruct = NULL;
				break;
			}
		}
	}
	
	if( indiv!=NULL )
		MEM_DELETE_ARRAY(indiv); // indiv has length params.nindivs

	if( newindiv!=NULL )
		MEM_DELETE_ARRAY(newindiv); // newindiv has length params.nindivs

	if( cumfit!=NULL ) {
		for( int i = 0; i < params->nindivs; i++ )
			MEM_DELETE_ARRAY(cumfit[i]); // cumfit[i] has length 2
		MEM_DELETE_ARRAY(cumfit); // cumfit has length params.nindivs
	}

	if( logfname!=NULL )
		MEM_DELETE_ARRAY(logfname); // logfname has length strlen(logfname)+1
		
	if( treeString!=NULL)
		MEM_DELETE_ARRAY(treeString);

	if(topologies!=NULL){
		for(int i=0;i<params->nindivs;i++){
			if(*(topologies+i)!=NULL){
				delete *(topologies+i);
				}
			}
		delete []topologies;
		}
		
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

void Population::Setup(const Parameters& params_, GeneralGamlConfig *conf, int nprocs, int r){
	rank=r;
	params = const_cast<Parameters*>(&params_); 

	stopwatch.Start();

	//this is where most of the allocation occurs
	int i, N;
	char nindivs_str[81];
	char ch;

	adap=new Adaptation(conf);
	refineStart=conf->refineStart;
	outputTreelog=conf->outputTreelog;
	outputMostlyUselessFiles = conf->outputMostlyUselessFiles;
	outputPhylipTree = conf->outputPhylipTree;
	
	enforceTermConditions = conf->enforceTermConditions;
	lastTopoImproveThresh = conf->lastTopoImproveThresh;
	improveOverStoredIntervalsThresh = conf->improveOverStoredIntervalsThresh;
	
	bootstrapReps = conf->bootstrapReps;
	inferInternalStateProbs = conf->inferInternalStateProbs;

	InitializeOutputStreams();
	
	subtreeNode=0;
	
	params->data->MakeAmbigStrings();
	//setup the bipartition statics
	Bipartition::blockBits=sizeof(int)*8;
	Bipartition::ntax=params->data->NTax();
	Bipartition::nBlocks=(int)ceil((double)params->data->NTax()/(double)Bipartition::blockBits);
	Bipartition::largestBlockDigit=1<<(Bipartition::blockBits-1);
	Bipartition::allBitsOn=(unsigned int)pow(2.0, Bipartition::blockBits)-1;
	Bipartition::str=new char[params->data->NTax()+1];
	Bipartition::str[params->data->NTax()] = '\0';

	//make sure that this is set before allocating any models
	Model::noPinvInModel = conf->dontInferProportionInvariant;

	//this is a really cheap hack
	Bipartition tmp;
	tmp.SetPartialBlockMask();

	if(rank == 0) 
		total_size = params->nindivs + nprocs-1;
	else
		total_size = params->nindivs;

	int max_indivs = total_size;

	//allocate the treeString
	//remember that we also encode internal node numbers sometimes
	double taxsize=log10((double) ((double)params->data->NTax())*params->data->NTax()*2);
	stringSize=(int)((params->data->NTax()*2)*(8+DEF_PRECISION)+taxsize);
	treeString=new char[stringSize];
	stringSize--;
	treeString[stringSize]='\0';

	//allocate the indiv array
	MEM_NEW_ARRAY(indiv,Individual,total_size);

	MEM_NEW_ARRAY(newindiv,Individual,total_size);

	for( i = 0; i < total_size; i++ ){
		newindiv[i].SetParams( params );
		newindiv[i].mod->SetParams( params );
		indiv[i].SetParams(params);
		indiv[i].mod->SetParams( params );
		}

	MEM_NEW_ARRAY(cumfit,double*,total_size);

	for( i = 0; !error && i < total_size; i++ ) {
		MEM_NEW_ARRAY(cumfit[i],double,2);
		}

	int max_topos = max_indivs;
	topologies=new TopologyList*[max_topos];
	TopologyList::SetIndL(indiv);
	ntopos=total_size;
	for(i=0;i<max_topos;i++)
		{topologies[i]=new TopologyList();
		topologies[i]->Allocate(100<(max_topos +1) ? 100 : (max_topos+1));
		if(i<total_size){
			topologies[i]->AddInd(i);
			TopologyList::ntoposexamined++;
			indiv[i].topo=i;
			}
		}

	//instantiate the ParallelManager
	if(rank==0){
		MasterGamlConfig *mastConf = (MasterGamlConfig*) (conf);
		paraMan = new ParallelManager(params->data->NTax(), nprocs, mastConf);
		}

#ifdef INCLUDE_PERTURBATION
	pertMan = new PerturbManager(conf);
#else
	pertMan = new PerturbManager();
#endif

	//instantiate the clamanager and figure out how much memory to snatch
	double memToUse;
	if(conf->availableMemory > 0){
		cout << "\nTotal system memory specified as " << conf->availableMemory << " megs" << endl;
		memToUse=0.8*conf->availableMemory;
		}
	else{
		cout << "\nMemory to be used for conditional likelihood arrays specified as " << conf->megsClaMemory << " megs" << endl;
		memToUse=conf->megsClaMemory;
		}
		
	const int KB = 1024;
	const int MB = KB*KB;
	int claSizePerNode = (4 * indiv[0].mod->NRateCats() * params->data->NChar() * sizeof(double)) + (params->data->NChar() * sizeof(int));
	int numNodesPerIndiv = params->data->NTax()-2;
	int sizeOfIndiv = claSizePerNode * numNodesPerIndiv;
	int idealClas =  3 * total_size * numNodesPerIndiv;
	int maxClas = (int)((memToUse*MB)/ claSizePerNode);
	int numClas;

	int L0=(int) (numNodesPerIndiv * total_size * 2);//a downward and one upward set for each tree
	int L1=(int) (numNodesPerIndiv * total_size + 2*total_size + numNodesPerIndiv); //at least a downward set and a full root set for every tree, plus one other set
	int L2=(int) (numNodesPerIndiv * 2.0 + 2*total_size);//a downward set for the best, one other full set and enough for each root direction
	int L3=(int) (numNodesPerIndiv * 1.5 - 2 + 2*total_size);//one full set, enough to reserve at least all of the full internals of the 
													 //best indiv and enough for each root
	if(maxClas >= L0){
		numClas=idealClas;
		memLevel = 0;		
		}
	else{
		numClas=maxClas;
	 	if(maxClas >= L1) memLevel = 1;
	 	else if(maxClas >= L2) memLevel = 2;
	 	else if(maxClas >= L3) memLevel = 3;
	 	else memLevel=-1;
		}

	cout.precision(4);
	cout << "allocating memory...\nusing " << (double)numClas*(double)claSizePerNode/(double)MB << " megs for conditional likelihood arrays.  Memlevel=" << memLevel << endl;
	cout << "For this dataset:\n";
	cout << "level 0: >= " << ceil(L0 * (claSizePerNode/(double)MB)) << " megs\n";
	cout << "level 1: " << ceil(L0 * ((double)claSizePerNode/MB))-1 << " to " << ceil(L1 * ((double)claSizePerNode/MB)) << " megs\n";
	cout << "level 2: " << ceil(L1 * ((double)claSizePerNode/MB))-1 << " to " << ceil(L2 * ((double)claSizePerNode/MB)) << " megs\n";
	cout << "level 3: " << ceil(L2 * ((double)claSizePerNode/MB))-1 << " to " << ceil(L3 * ((double)claSizePerNode/MB)) << " megs\n";
	cout << "not enough mem: <= " <<  ceil(L3 * ((double)claSizePerNode/MB))-1 << endl << endl;

	if(memLevel==-1){
		throw ErrorException("Not enough memory specified in config file (megsclamemory)!");
		}

	//increasing this more to allow for the possiblility of needing a set for all nodes for both the indiv and newindiv arrays
	//if we do tons of recombination 
	idealClas *= 2;
	//allocate an extra bunch of claHolders to allow for temporary trees that are necessary during perturbation
	claMan=new ClaManager(params->data->NTax()-2, numClas, idealClas, params->data->NChar(), indiv[0].mod->NRateCats());
	//set the trees static pointer to the clamanager
	Tree::claMan=claMan;
	Tree::data=params->data;
	Tree::rescaleEvery=16;
	Tree::meanBrlenMuts	= params->meanBrlenMuts;
	Tree::alpha		= params->gammaShapeBrlen;
	Tree::treeRejectionThreshold = params->treeRejectionThreshold;
	Model::mutationShape = params->gammaShapeModel;
	
	if( params->restart )
	{
		// restarting: read individuals from file having name params->statefname
		ifstream gmlf( params->statefname );

		// ignore first line this time; contains gen, time and seed and was already read by GetRestartParams in gamlmain.cpp
		
		gmlf.get( nindivs_str, 80, '\n' );
		gmlf.get(ch);
		assert( ch == '\n' );

		gmlf.get( nindivs_str, 80, '\n' );
		gmlf.get(ch);
		assert( ch == '\n' );

		N = atoi(nindivs_str);
		assert( N == params->nindivs );

		for( i = 0; i < params->nindivs; i++ ) {
			indiv[i].ReadTreeFromFile(gmlf);
			indiv[i].SetParams( params );
		}

		gmlf.close();
		ofstream tmpf( "restart_test.tre" );
		tmpf << "#nexus" << endl << endl;
		tmpf << "begin trees;" << endl;

		for( i = 0; i < params->nindivs; i++ ) {
			tmpf.setf( ios::floatfield, ios::fixed );
			tmpf.setf( ios::showpoint );
			tmpf << "tree gaml" << i << " = " << MakeNewick(i, false) << endl;
		}

		tmpf << "end;" << endl;
		tmpf.close();

	}
	else
	{
		if(bootstrapReps==0) SeedPopulationWithStartingTree();
	}

#ifndef UD_VERSION
	if(bootstrapReps==0) AppendTreeToTreeLog(-1, -1);
#endif

	if( !error ) is_setup = 1;

	current_size = original_size = total_size;
	
	for (i = params->nindivs; i < total_size; i++)	{
		indiv[i].reproduced = indiv[i].willreproduce = 1;
		newindiv[i].reproduced = newindiv[i].willreproduce = 1;
		indiv[i].parent=i;
		newindiv[i].parent=i;	
		}
}

void Population::ResetMemLevel(int numNodesPerIndiv, int numClas){
	const int KB = 1024;
	const int MB = KB*KB;
	
	int claSizePerNode = (4 * indiv[0].mod->NRateCats() * params->data->NChar() * sizeof(double)) + (params->data->NChar() * sizeof(int));
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
	}


void Population::SeedPopulationWithStartingTree(){
	
	if(strcmp(params->startfname, "random"))
		cout << "Obtaining starting conditions from file " << params->startfname << endl;
	
	//claMan->MakeAllHoldersDirty();
	for(int i=0;i<total_size;i++){
		if(indiv[i].treeStruct != NULL) indiv[i].treeStruct->RemoveTreeFromAllClas();
		if(newindiv[i].treeStruct != NULL) newindiv[i].treeStruct->RemoveTreeFromAllClas();
		}
	
	//if starting from a treefile, use the treestring
	//to create the first indiv, and then copy the tree and clas
	indiv[0].Randomize(params->startfname, rank);
	indiv[0].treeStruct->root->CheckforPolytomies();
	indiv[0].treeStruct->CheckBalance();
	indiv[0].treeStruct->mod=indiv[0].mod;

	cout.precision(10);
	cout << "Initial ln Likelihood: " << indiv[0].Fitness() << endl;
	
	indiv[0].treeStruct->CalcBipartitions();	
	
	if(refineStart==true){
		indiv[0].RefineStartingConditions(!(adap->modWeight==0.0), adap->branchOptPrecision);
		indiv[0].CalcFitness(0);
		cout << "lnL after optimization: " << indiv[0].Fitness() << endl;
		}	

	globalBest=bestFitness=prevBestFitness=indiv[0].Fitness();
	for(int i=1;i<total_size;i++){
		if(indiv[i].treeStruct==NULL) indiv[i].treeStruct=new Tree();
		indiv[i].CopySecByRearrangingNodesOfFirst(indiv[i].treeStruct, &indiv[0]);
		indiv[i].treeStruct->mod=indiv[i].mod;
		}
	UpdateTopologyList(indiv);
	CalcAverageFitness();
	}

int Population::Restart(int type, int rank, int nprocs, int count)	{
	//not working currently
	assert(0);
	if (type == 1)	{
		/* don't need this stuff
		EliminateDuplicateTreeReferences(); 			
		for (int i = 0; i < total_size; ++i)	{
			for (int j = 0; j < total_size; ++j)	{
				if (newindiv[i].treeStruct == indiv[j].treeStruct)	{
					newindiv[i].treeStruct = NULL;
					break;
				}
			}
		}
		*/
		bool readFromFile=false;
//		claMan->ClearAllClas();
		delete indiv[0].treeStruct;
		//fprintf(stderr, "entering Restart\n");
/*		if (!FileExists( params->startfname ))	{
			printf("starting tree file \"%s\" does not exist, using random start tree\n", params->startfname);
			indiv[0].Randomize("random", rank, sharedcl);
		}
		else{
			char temp[10000];
			ifstream stf( params->startfname, ios::in );
			if (!stf)	{
				printf("could not open file\"%s\", using random start tree\n", params->startfname);
				indiv[0].Randomize("random", rank, sharedcl);
			}
			else{
				for(int r=-1;r<nprocs*count+rank;r++)	{
					stf >> temp;//get the score
					stf >> temp;//get the tree or maybe kappa
					if(stf.eof())
						break;
					if(isdigit(temp[0])){
						assert(0);
						//indiv[0].kappa=atof(temp);
						stf >> temp;//get the tree
					}
				}
				if (stf.eof())	{
					printf("not enough trees in file, using random tree.\n");
					indiv[0].Randomize("random", rank, sharedcl); // carefull, this allocates a new treeStruct...
					params->starting_tree = "";
				}
				else{
					params->starting_tree=temp;
					readFromFile=true;
					}
			}
			stf.close();
		}
		for (int i = 0; i < total_size; ++i){
			if(i==0 && readFromFile){ // we don't want to be allocating a new treeStruct if we called Individual::Randomize()
				indiv[0].treeStruct=new Tree( params->starting_tree.c_str());
				indiv[0].treeStruct->AssignCLAsFromMaster();
				indiv[0].treeStruct->CheckBalance();		
				}		
			else indiv[i].CopySecByRearrangingNodesOfFirst(indiv[i].treeStruct, &indiv[0]);
			indiv[i].ResetIndiv();
			indiv[i].SetDirty();
			}

		allTimeBest=&indiv[0];
		allTimeBest->CalcFitness(0);
		bestFitness=prevBestFitness=allTimeBest->Fitness();
		
		//reset the topology list.  All individuals will be the same topology
		delete []topologies;
		int max_indivs=total_size;
		topologies=new TopologyList*[max_indivs];
		TopologyList::SetIndL(indiv);
		ntopos=total_size;
		for(int i=0;i<max_indivs;i++){
			topologies[i]=new TopologyList();
			topologies[i]->Allocate(100<(max_indivs +1) ? 100 : (max_indivs+1));
			topologies[0]->AddInd(i);
			indiv[i].topo=0;
			}
*//*		for(int i=0;i<max_indivs;i++){
			topologies[i]=new TopologyList();
			topologies[i]->Allocate(100<(max_indivs +1) ? 100 : (max_indivs+1));
			if(i<total_size){
				topologies[i]->AddInd(i);
				TopologyList::ntoposexamined++;
				indiv[i].topo=i;
				}
			}
*/		}	
	return 0;
}

void Population::Run(){

	calcCount=0;
	optCalcs=0;

	cout << "Running Genetic Algorithm with initial seed=" << rnd.init_seed() << endl;
	
	avgfit = CalcAverageFitness();

	#ifndef NO_OUTPUT
	cout.setf(ios::left);
	cout.precision(6);
	cout << "\t" << setw(10) <<  "gen" << setw(15) << "current lnL" << setw(10) << "precision" << setw(10) << "lastChange"<< endl;			
	#endif

	gen=0;
	OutputLog();

	CatchInterrupt();

	for (gen = 1; gen < params->stopgen+1; ++gen){
		
		NextGeneration();
		keepTrack();
		if(outputMostlyUselessFiles) OutputFate();		 
		if(!(gen % params->logEvery)) OutputLog();
		if(!(gen % params->saveEvery)){
			if(bootstrapReps==0) CreateTreeFile( params->treefname);
			#ifndef NO_OUTPUT
			cout.setf(ios::left);
			cout << "\t" << setw(10) << gen << setw(15) << setprecision(10) << indiv[bestIndiv].Fitness() << setprecision(4) << setw(10) << adap->branchOptPrecision << setw(10) << lastTopoImprove << endl;			
			#endif
			}
		if(askQuitNow == 1){
			cout << "Perform final branch-length optimization and terminate now? (y/n)" << endl;
			char c=getchar();
			cin.get();
			if(c=='y') break;
			else{
				askQuitNow = 0;
				CatchInterrupt();
				cout << "continuing ..." << endl;
				cin.get();
				}
			}

#ifdef PERIODIC_SCORE_DEBUG
		if(gen % 500 == 0 ||gen==1)
			OutputFilesForScoreDebugging(&indiv[bestIndiv], tempGlobal++);
#endif

#ifdef NNI_SPECTRUM
		if(gen % 1000 == 0 || gen==1)
			NNISpectrum(bestIndiv);
#endif
		if(!(gen%adap->intervalLength)){
			cout.precision(10);
			bool reduced=false;
			if(gen-lastTopoImprove > adap->intervalsToStore*adap->intervalLength){
				reduced=adap->ReducePrecision();
				}
			if(reduced){
				lastTopoImprove=gen;
				indiv[bestIndiv].treeStruct->OptimizeAllBranches(adap->branchOptPrecision);
				indiv[bestIndiv].SetDirty();
				CalcAverageFitness();
#ifndef UNIX
				cout << "optimization precision reduced, optimizing ...\t" << bestFitness << "->";
				cout << bestFitness << endl;
#endif
				}
			else if(!(gen%(2*adap->intervalLength*adap->intervalsToStore))){
#ifndef UNIX
//					cout << "optimizing ...\t" << bestFitness << "->";
#endif
					indiv[bestIndiv].treeStruct->OptimizeAllBranches(adap->branchOptPrecision);
					indiv[bestIndiv].SetDirty();
					CalcAverageFitness();
//					cout << bestFitness << endl;				
					}
			
			//termination conditions
			if(enforceTermConditions == true
				&& gen-lastTopoImprove > lastTopoImproveThresh 
				&& adap->improveOverStoredIntervals < improveOverStoredIntervalsThresh
				&& adap->branchOptPrecision == adap->minOptPrecision){
				cout << "Reached termination condition!\nlast topological improvement at gen " << lastTopoImprove << endl;
				cout << "Improvement over last " << adap->intervalsToStore*adap->intervalLength << " gen = " << adap->improveOverStoredIntervals << endl;
				break;
				}
#ifdef INCLUDE_PERTURBATION
			CheckPerturbSerial();
#endif
			}
		if(params->stoptime - stopwatch.SplitTime() < 120){
			cout << "time limit of " << params->stoptime << " seconds reached..." << endl;
			break;
			}
#ifdef INCLUDE_PERTURBATION
		if(pertMan->pertAbandoned==true && pertMan->restartAfterAbandon==true && (gen - pertMan->lastPertGeneration > pertMan->gensBeforeRestart)){
			params->starting_tree="";
			pertMan->lastPertGeneration=gen;
			pertMan->pertAbandoned=false;
			pertMan->numPertsNoImprove=0;
			bestSinceRestart.SetFitness(-1e100);
			if(indiv[bestIndiv].Fitness() > allTimeBest.Fitness()) StoreAllTimeBest();
			SeedPopulationWithStartingTree();
			cout << "restarting ...." << endl;
			}
#endif
		}
	FinalOptimization();
	gen=-1;
	OutputLog();
	if(bootstrapReps==0) FinalizeOutputStreams();
	
	if(inferInternalStateProbs == true){
		cout << "Inferring internal state probabilities...." << endl;
		indiv[bestIndiv].treeStruct->InferAllInternalStateProbs(params->ofprefix);
		}
		
	cout << "finished" << endl;
		
//	log << calcCount << " cla calcs, " << optCalcs << " branch like calls\n";
	}

void Population::FinalOptimization(){
	cout.setf(ios::fixed);
	cout.precision(5);
	cout << "Current score = " << indiv[bestIndiv].Fitness() << endl;
	
	if(pertMan->ratcheted) TurnOffRatchet();
	
	if(indiv[bestIndiv].Fitness() < allTimeBest.Fitness()){
		RestoreAllTimeBest();
		}
	
	for(int i=0;i<total_size;i++){
		if(i != bestIndiv) indiv[i].treeStruct->RemoveTreeFromAllClas();
		}
	
	cout << "Performing final branch optimization..." << endl;
	int pass=1;
	double incr;
	do{
		incr=indiv[bestIndiv].treeStruct->OptimizeAllBranches(adap->branchOptPrecision * pow(0.5, pass));
		indiv[bestIndiv].SetDirty();
		indiv[bestIndiv].CalcFitness(0);
		cout << "\tpass " << pass++  << " " << indiv[bestIndiv].Fitness() << endl;
		}while(incr > .00001 || pass < 10);
	cout << "Final score = " << indiv[bestIndiv].Fitness() << endl;
	log << "Score after final optimization: " << indiv[bestIndiv].Fitness() << endl;
	CreateTreeFile( params->treefname );
	cout.unsetf(ios::fixed);
	}

void Population::Bootstrap(){
	
	params->data->ReserveOriginalCounts();

	stopwatch.Start();
	CatchInterrupt();

	for(int rep=1;rep<=bootstrapReps;rep++){
		params->data->BootstrapReweight();
		
		cout << "bootstrap replicate " << rep << endl;
		SeedPopulationWithStartingTree();
		Run();

		doneWithRep = false;
		adap->branchOptPrecision = adap->startOptPrecision;
		AppendTreeToBootstrapLog(rep);
		cout << "finished with bootstrap rep " << rep << endl;
		}
	FinalizeOutputStreams();
	}


void Population::QuickSort( int top, int bottom )
{
	int i = top;
	int j = bottom;
	double x = cumfit[ (top + bottom) / 2 ][1];
	do {
		while( cumfit[i][1] < x  &&  i < bottom ) i++ ;
		while( x < cumfit[j][1]  &&  j > top  ) j-- ;

		if( i <= j ) {
			for( int k = 0; k < 2; k++ ) {
				double y = cumfit[i][k];
				cumfit[i][k] = cumfit[j][k];
				cumfit[j][k] = y;
			}
			i++;
			if(j) j--;
		}

	} while( i <= j );

	if( top  <    j    ) QuickSort( top, j );
	if(  i   <  bottom ) QuickSort( i, bottom );
}

double Population::CalcAverageFitness(){
	double total = 0.0;
	
	for(int i = 0; i < total_size; i++ ){
		// evaluate fitness
		if(indiv[i].IsDirty()){
			indiv[i].CalcFitness(subtreeNode);
			}
		assert(indiv[i].Fitness() != 1);
	
		total += indiv[i].Fitness();
		cumfit[i][0] = i;
		cumfit[i][1] = indiv[i].Fitness();
		}

	double avg = total / (double)total_size;

	//if we're in subtree mode, give the indivs without accurate subtrees a really
	//crappy lnL so they don't reproduce
/*	if(rank==0)
		if(paraMan->subtreeModeActive==true){
			for(int i=0;i<total_size;i++)
				if(indiv[i].accurateSubtrees==false) cumfit[i][1]=-1e100;	
			}
*/
	// Sort fitnesses from low to high (bad to good)
	QuickSort( 0, total_size-1 );

	// keep track of which individual is most fit each generation we've stored the 
	//fitnesses as ln-likelihoods in cumfit, so cumfit[0] will be the _least_ fit individual
	int mostFit = total_size-1;
	bestAccurateIndiv=bestIndiv = (int)cumfit[mostFit][0];
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
		globalBest=bestFitness = cumfit[mostFit][1];
		new_best_found = 1;
		numgensamebest=0;
		}
	else numgensamebest++;



//	allTimeBest = &indiv[bestIndiv];

	if(memLevel>0){
		//if we are at some level of memory restriction, mark the clas of the old best
		//for reclamation, and protect those of the new best
		SetNewBestIndiv(bestIndiv);
		}

	//probability of reproduction based on more or less on AIC weights, although
	//the strength of selection can be varied by changing the selectionIntensity
	//A selectionIntensity of 1.0 makes this equivalent to AIC weights, while 
	//smaller number make the selection less severe
	double *deltaAIC=new double[total_size];
	double tot=0.0;

	for(int i=0;i<total_size-1;i++){
		deltaAIC[i]=cumfit[total_size-1][1] - cumfit[i][1];
		deltaAIC[i]=exp(-params->selectionIntensity * deltaAIC[i]);
		tot+=deltaAIC[i];
		}

	deltaAIC[total_size-1]=params->holdoverPenalty;
	deltaAIC[total_size-1]=exp(-params->selectionIntensity * deltaAIC[total_size-1]);
	tot+=deltaAIC[total_size-1];

	for(int i=0;i<total_size;i++)
		deltaAIC[i] /= tot;
	
	double cum=deltaAIC[0];
	cumfit[0][1] = cum;
	for(int i = 1; i < total_size; i++ ) {
		cum += deltaAIC[i];
		cumfit[i][1] = cum;
		}
	delete []deltaAIC;

//only allow the best indiv to reproduce
/*	for(int i = 0; i < total_size; i++ ) {
		if(cumfit[i][0]==0) cumfit[i][1]=1.0;
		else cumfit[i][1]=0.0;
		}
	bestIndiv=0;		
*/

	//DEBUG
/*	if(gen > 1){
		bool same1, same2, same3;

		indiv[0].treeStruct->CalcBipartitions();
		indiv[1].treeStruct->CalcBipartitions();
		indiv[2].treeStruct->CalcBipartitions();
		indiv[3].treeStruct->CalcBipartitions();
		
		if(newindiv[indiv[1].parent].treeStruct){
			same1=newindiv[indiv[1].parent].treeStruct->IdenticalSubtreeTopology(indiv[1].treeStruct->root);
			if(same1 == false) assert(indiv[1].mutation_type & Individual::anyTopo);
			else assert(!(indiv[1].mutation_type & Individual::anyTopo));
			}

		if(newindiv[indiv[2].parent].treeStruct){
			same2=newindiv[indiv[2].parent].treeStruct->IdenticalSubtreeTopology(indiv[2].treeStruct->root);
			if(same2 == false) assert(indiv[2].mutation_type & Individual::anyTopo);
			else assert(!(indiv[2].mutation_type & Individual::anyTopo));
			}

		if(newindiv[indiv[3].parent].treeStruct){
			same3=newindiv[indiv[3].parent].treeStruct->IdenticalSubtreeTopology(indiv[3].treeStruct->root);
			if(same3 == false) assert(indiv[3].mutation_type & Individual::anyTopo);
			else assert(!(indiv[3].mutation_type & Individual::anyTopo));
			}
			
		ofstream deb("ident.log", ios::app);
		deb << gen << "\n" << same1 << "\n" << same2 << "\n" << same3 << "\n";
		deb.close();
		}
*/
	return avg;
	
/*	Here's Paul's original selection criterion, based solely on rank
	//	
	// relative fitnesses are assigned based solely on position
	// of individual in sorted array - we forget the likelihoods (or treelengths)
	// at this point.  This allows the likelihoods to be close together
	// and still get a healthy distribution of relative fitnesses so that
	// there is real differential reproduction

	double n = (double)total_size;
	double nn = n * ( n + 1.0 );
	double incr = 2.0 / nn;
	double cum = incr;

	cumfit[0][1] = cum;
	for( i = 1; i < total_size; i++ ) {
		cum += incr;
		cumfit[i][1] = cumfit[i-1][1] + cum;
	}
*/
}

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
	gnuf << "N=" << params->nindivs;
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

void Population::DetermineParentage(){
	//determine each individual's parentage
	int parent;
	double r;
	for(int i = 0; i < params->nindivs; i++ ){
		if( i < params->holdover ){// copy best individual's genotype to next generation
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


#ifdef MPI_VERSION
//new bipart recom conditions, 9-25-05

			if(rank==0 && paraMan->subtreeModeActive==false && i>= (params->nindivs-paraMan->maxRecomIndivs)){
				int *mates=new int[paraMan->nremotes];
				for(int j=0;j<paraMan->nremotes;j++) mates[j]=params->nindivs+j;
				ScrambleArray(paraMan->nremotes, mates);
				int mateIndex=0;
				int curMate;
				// find someone else to recombine with
				do{
					curMate=mates[mateIndex++];

					}while(mateIndex < paraMan->nremotes && 
						(curMate==parent //don't recombine with your parent
						|| (indiv[parent].topo == indiv[curMate].topo) //don't recombine with another of the same topo	
						|| (indiv[curMate].willrecombine == true)));//don't recombine with someone who is already doing so
			if(mateIndex < paraMan->nremotes){
					newindiv[i].recombinewith=curMate;
					indiv[curMate].willrecombine=true;
					//this will be a new topology, so mark it as topo -1.  This will be dealt with when we update the topolist
					newindiv[i].topo=-1;
					}
				delete []mates;
				}
#endif
			}
		
		newindiv[i].parent=parent;
		if(newindiv[i].mutation_type==Individual::subtreeRecom) newindiv[i].topo=-1; //VERIFY
		else newindiv[i].topo=indiv[parent].topo;
		indiv[ parent ].willreproduce=true;
		}
	}

void Population::FindTreeStructsForNextGeneration(){
	//find treestructs for all of the newindivs, either by getting an unused one from the previous
	//generation or by getting one from the unusedTree stack
	for(int i = 0; i < total_size; i++ ){
		//see if the parent indiv has already been used in the new generation, or if it will recombine
		if( i < params->nindivs && (indiv[newindiv[i].parent].reproduced||indiv[newindiv[i].parent].willrecombine )){	      
			//See if there is another ind with the same topology that also will not recombine
			int sot=-1;
			if(topologies[indiv[newindiv[i].parent].topo]->nInds>1)//if this isn't the only individual of this topo
				do{
					sot=topologies[indiv[newindiv[i].parent].topo]->GetNumberOfUnselectedInd();
					if( !(sot==newindiv[i].parent) && !(indiv[sot].reproduced) && !(indiv[sot].willrecombine)) break;
					}while(sot!=-1);
			if(sot>=0){
				newindiv[i].CopySecByStealingFirstTree(&indiv[sot],&indiv[newindiv[i].parent]);
				indiv[sot].reproduced=true;
				}
			else{
				//DZ 7-5 rewriting this.  If no unused tree with the same topology exists, use a tree from the 
				//unused Indiv stack.  If it is empty, create an extra indiv that will eventually make it's way 
				//back to that stack.  At most we should only ever have nindiv trees in the unused stack
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
			}
		else{
			//if the tree will not be used in recombination and has not already been used
			newindiv[i].CopyByStealingTree(&indiv[newindiv[i].parent]);
			indiv[ newindiv[i].parent].reproduced=true;
			if(i>params->nindivs) newindiv[i].mutation_type=indiv[i].mutation_type;
			}
		}
	}
	
void Population::PerformMutation(int indNum){
	Individual *ind=&newindiv[indNum];
	Individual *par=&indiv[newindiv[indNum].parent];

	double beforeScore;
	bool recomPerformed;

	switch(ind->mutation_type){
		case Individual::exNNI: //exNNI and exlimSPR trump all other mutation types
			beforeScore=par->Fitness();
			cout<< "\t EXNNI called at gen: "<< 
			  gen<<" : l="<<params->holdover<<endl;
			NNIoptimization(indNum, 1);
			if(beforeScore==ind->Fitness()){
				topologies[ind->topo]->exNNItried=true;
				}
			//ind->accurateSubtrees=false;
			break;
		
		case Individual::exlimSPR:
			SPRoptimization(indNum);
			ind->accurateSubtrees=false;
			break;
		
		case Individual::subtreeRecom:
			//perform subtree recom, which melds together the different subtrees worked on by the
			//remote nodes
			recomPerformed=SubtreeRecombination(indNum);
			if(recomPerformed==false) ind->mutation_type=0;
			ind->treeStruct->calcs=calcCount;
			calcCount=0;
			ind->CalcFitness(0);
			break;
			
		default:
			if(ind->recombinewith>-1){// perform recombination
				Individual *recompar=&indiv[ind->recombinewith];
				ind->treeStruct->CalcBipartitions();
				recompar->treeStruct->CalcBipartitions();
				ind->CrossOverWith( *recompar, adap->branchOptPrecision);
				ind->accurateSubtrees=false;
				ind->treeStruct->calcs=calcCount;
				calcCount=0;
				}
			if(ind->recombinewith==-1){//all types of "normal" mutation that occur at the inidividual level
				if(rank==0){//if we are the master
				 	if(ind->accurateSubtrees==false || paraMan->subtreeModeActive==false){
				 		//DEBUG
				 		//this is really cheesey, but ...
/*				 		int savedSeed=rnd.seed();
				 		double r=rnd.uniform();
				 		if(r < adap->topoMutateProb)
				 			OutputFilesForScoreDebugging(ind, tempGlobal++);
				 		rnd.set_seed(savedSeed);
*/				 	
			       		ind->Mutate(adap->branchOptPrecision, adap);
			       		//reclaim clas if the created tree has essentially no chance of reproducing
			       		if(((ind->Fitness() - indiv[bestIndiv].Fitness()) < (-11.5/params->selectionIntensity))){
			       			ind->treeStruct->ReclaimUniqueClas();
			       			}
			       		//DEBUG
			       		//if(ind->mutation_type & Individual::anyTopo) OutputFilesForScoreDebugging(ind, tempGlobal++);
						}
					else{
						//if subtree mode is on and we are the master, mutate one of the nodes
						//that isn't in a subtree, or alternatively pick a subtree and mutate it
					#ifndef MASTER_DOES_SUBTREE
						if(paraMan->fewNonSubtreeNodes != true)
							ind->NonSubtreeMutate(paraMan, adap->branchOptPrecision, adap);
						else 
							ind->SubtreeMutate(subtreeNode, adap->branchOptPrecision, subtreeMemberNodes, adap);
					#else
						ind->SubtreeMutate(subtreeNode, adap->branchOptPrecision, subtreeMemberNodes, adap);
					#endif					
						}
					}
				else{//if we are a remote node
				 	if(subtreeNode==0) ind->Mutate(adap->branchOptPrecision, adap);
					else{
						ind->SubtreeMutate(subtreeNode, adap->branchOptPrecision, subtreeMemberNodes, adap);
						}
					}
				}
			}
		
		//check the accuracy of the subtrees
		#ifndef NDEBUG
		if(rank==0 && ind->accurateSubtrees==true)
			paraMan->CheckSubtreeAccuracy(ind->treeStruct);
		#endif

		if((ind->mutation_type & Individual::anyTopo) || (ind->mutation_type & Individual::rerooted))
			AssignNewTopology(newindiv, indNum);		
		}

//note that we're passing the entire array of individuals here, not just a pointer to an individual
void Population::AssignNewTopology(Individual *indArray, int indNum){
	assert(indArray == topologies[0]->GetListOfInd());
	Individual *ind = &indArray[indNum];
	if(ind->topo == -1){
		ind->topo=ntopos++;
		topologies[ind->topo]->AddInd(indNum);
		}
	else if(topologies[ind->topo]->nInds>1){
		topologies[ind->topo]->RemoveInd(indNum);
		ind->topo=ntopos++;
		topologies[ind->topo]->AddInd(indNum);
		}
	topologies[ind->topo]->gensAlive=0;
	assert(topologies[ind->topo]->nInds==1);
	TopologyList::ntoposexamined++;
	}

void Population::NextGeneration(){

	new_best_found = 0;

	DetermineParentage();

	for(int i=0;i<ntopos;i++)
		topologies[i]->NewGeneration();

	FindTreeStructsForNextGeneration();

	UpdateTopologyList(newindiv);

	//return any treestructs from the indivs that won't be used in recombination
	//and weren't used to make the newindivs.  This is necessary to keep from having
	//too many CLAs in use at any one time 
	for(int j=0;j<params->nindivs;j++){
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
	for(int indnum = params->holdover; indnum < params->nindivs; indnum++ ){
		PerformMutation(indnum);
		}

	UpdateTopologyList(newindiv);
	UpdateTreeModels();

	//the only trees that we need to return at this point are ones that
	//did not reproduce AND were used in recom.  Those that weren't used 
	//in recom were already reclaimed above, and the treestructs set to NULL
	for(int j=0;j<params->nindivs;j++){
		if(indiv[j].reproduced==false && indiv[j].treeStruct!=NULL){
			indiv[j].treeStruct->RemoveTreeFromAllClas();
			unusedTrees.push_back(indiv[j].treeStruct);
			}
		//reset all of the individuals
		indiv[j].ResetIndiv();
		}

	// swap newindiv and indiv
	for(int i=0;i<params->nindivs;i++)
		indiv[i].ResetIndiv();
	for(int i=params->nindivs;i<total_size;i++)
		indiv[i].willrecombine=false;

	Individual* tmp = newindiv;
	newindiv = indiv;
	indiv = tmp;

	if( params->showProgress )
		cout << "Calculating average fitness" << endl;
	avgfit = CalcAverageFitness(); //score individuals that need it
		
	#ifdef DEBUG_SCORES
	if(rank==0)	OutputFilesForScoreDebugging();
	#endif
	
	}

void Population::OutputFate(){
	//output everything that happened to each indiv in this generation to file

	for(int i=0;i<total_size;i++){
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

void Population::CreateStateFile(){
assert(0);
	ofstream outf( params->statefname );
	outf<<setprecision(DEF_PRECISION);

   // save current generation (note: this is one less than the generation number
   // reported in the log file since this gen starts at 0 and the log file
   // generation count starts at 1).
   outf << gen << '\t';

   // save total elapsed time, which is sum of elapsed time from previous runs
   // and current elapsed time
   double seconds = stopwatch.SplitTime();

	outf << ( params->prev_time + seconds ) << '\t';

   // save current random number seed
   outf << rnd.seed() << '\t';
   outf << endl;

	// save current population data
//	outf << (*this) << endl;
	outf.close();
}

int Population::TimeToQuit()
{
	int quit = 0;

#if !defined( POWERMAC_VERSION )
	// try to open the file 'stop_now'
	// if this file does not exist (in the current directory)
	//	then GAML will continue running
	// if this file does exist, GAML will quit immediately
#if 0
	ifstream testf( "stop_now" );
	if( testf ) {
		quit = 1;
		cout << endl << programName << " is shutting down..." << endl;
	}
	else {
		cout << "  " << (gen+1) << " generations" << endl;
	}
	cout.flush();
	testf.close();
#else
	if( FileExists( "stop_now" ) ) {
		quit = 1;
		cout << endl << programName << " is shutting down..." << endl;
	}
	else {
		#ifndef UD_VERSION
		cout << "  " << (gen+1) << " generations" << endl;
		#endif
	}
#endif
#endif

	return quit;
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

	if(gen==1 && ind==NULL || num==1){
	
		outf << "#nexus" << endl << endl;
		outf << "begin trees;" << endl;
		TranslateTable tt( params->data );
		outf << tt << endl;
		
		paupf << "#nexus\n\n";
		paupf << "begin paup;\n";
		paupf << "set warnreset=no incr=auto;\n";
		paupf << "execute " << params->ofprefix << ".nex;\n";
#ifndef NNI_SPECTRUM
		paupf << "gett file=toscore.tre storebr;" << endl;
#else
		paupf << "gett file=" << outf << " storebr;" << endl;
#endif
		}
		
	if(ind==NULL){
		for(int i=0;i<total_size;i++){
			outf << "  utree " << gen << i << "= ";
			indiv[i].treeStruct->root->MakeNewick(treeString, false);
			outf << treeString << ";\n";
			
			paupf << "lset userbr ";
			if(indiv[i].mod->Nst()==2) paupf << "nst=2 trat=" << indiv[i].mod->Rates(0) << " base=(" << indiv[i].mod->Pi(0) << " " << indiv[i].mod->Pi(1) << " " << indiv[i].mod->Pi(2) << ");\n" << "lsc " << (gen-1)*params->nindivs+i+1;
			
			else paupf << "nst=6 rmat=(" << indiv[i].mod->Rates(0) << " " << indiv[i].mod->Rates(1) << " " << indiv[i].mod->Rates(2) << " " << indiv[i].mod->Rates(3) << " " << indiv[i].mod->Rates(4) << ") " << " base=(" << indiv[i].mod->Pi(0) << " " << indiv[i].mod->Pi(1) << " " << indiv[i].mod->Pi(2) << ") ";
			
			if(indiv[i].mod->NRateCats()>1) paupf << "rates=gamma shape=" << indiv[i].mod->Alpha() << " ";
			
			paupf << "pinv=" << indiv[i].mod->ProportionInvariant() << " "; 

			if(gen==1 && i==0) paupf << ";\n" << "lsc " << (gen-1)*total_size+i+1 << "/scorefile=paupscores.txt replace;\n";
			else paupf << ";\n" << "lsc " << (gen-1)*total_size+i+1 << "/scorefile=paupscores.txt append;\n";
			}
		}
	else{
		outf << "  utree " << num << "= ";
		ind->treeStruct->root->MakeNewick(treeString, false);
		outf << treeString << ";\n";
		
		paupf << "lset userbr ";
		if(ind->mod->Nst()==2) paupf << "nst=2 trat=" << ind->mod->Rates(0) << " base=(" << ind->mod->Pi(0) << " " << ind->mod->Pi(1) << " " << ind->mod->Pi(2) << ");\nlsc ";
		
		else paupf << "nst=6 rmat=(" << ind->mod->Rates(0) << " " << ind->mod->Rates(1) << " " << ind->mod->Rates(2) << " " << ind->mod->Rates(3) << " " << ind->mod->Rates(4) << ") " << " base=(" << ind->mod->Pi(0) << " " << ind->mod->Pi(1) << " " << ind->mod->Pi(2) << ") ";
		
		if(ind->mod->NRateCats()>1) paupf << "rates=gamma shape=" << ind->mod->Alpha() << " ";
		
		paupf << "pinv=" << ind->mod->ProportionInvariant() << " "; 
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

void Population::AppendTreeToTreeLog(int mutType, int indNum /*=-1*/){

	if(treeLog.is_open() == false || outputTreelog==false) return;

	Individual *ind;
	if(indNum==-1) ind=&indiv[bestIndiv];
	else ind=&indiv[indNum];

	treeLog << "  tree gen" << gen <<  "= [&U] [" << ind->Fitness() << "\tmut=" << mutType << "][ ";
	ind->mod->OutputGamlFormattedModel(treeLog);
	ind->treeStruct->root->MakeNewick(treeString, false);
	treeLog << "]" << treeString << ";" << endl;
	}


void Population::AppendTreeToBootstrapLog(int rep){

	if(bootLog.is_open() == false) return;

	Individual *ind = &indiv[bestIndiv];

	bootLog << "  tree bootrep" << rep <<  "= [&U] [" << ind->Fitness() << " ";
	
	ind->mod->OutputGamlFormattedModel(bootLog);
/*	
	const Model *m=ind->mod;
	if(m->Nst() == 2) bootLog << "nst=2 k=" << m->Rates(0) << "\t";
	else bootLog << "nst=6 rmat=(" << m->Rates(0) << " " << m->Rates(1) << " " << m->Rates(2) << " " << m->Rates(3) << "  " << m->Rates(4)  << ") ";
	bootLog << "base=(" << m->Pi(0) << " " << m->Pi(1) << " " << m->Pi(2) << ") ";

	if(m->NRateCats()>1) bootLog << "rates=g shape=" << m->Alpha() << " ";
	if(m->ProportionInvariant()!=0.0) bootLog << "pinv=" << m->ProportionInvariant();
*/	
	ind->treeStruct->root->MakeNewick(treeString, false);
	bootLog << "] " << treeString << ";" << endl;
	}


void Population::CreateTreeFile( const char* treefname, int fst /* = -1 */, int lst /* = -1 */ )
{
	int k;
	int first = ( fst < 0 ? 0 : fst );
	int last = ( lst < 0 ? params->nindivs : lst );
	assert( last <= params->nindivs );

	assert( treefname );
	ofstream outf;

	Individual *best=&indiv[bestIndiv];
	if(best->Fitness() < allTimeBest.Fitness() || pertMan->ratcheted==true) return;

	outf.open( treefname );
	outf.precision(8);
#ifndef UD_VERSION
	outf << "#nexus" << endl << endl;

	//rewritting this to output standard nexus tree files, not gamlviewer stuff
	int ntaxa = params->data->NTax();
	outf << "begin trees;\ntranslate\n";
	for(k=0;k<ntaxa;k++){
		outf << "  " << (k+1);
		NxsString tnstr = params->data->TaxonLabel(k);
		tnstr.blanks_to_underscores();
		outf << "  " << tnstr.c_str();
		if(k < ntaxa-1) 
			outf << ",\n";
		}		

	outf << ";\n";
	
	outf << "tree best = [&U][" << best->Fitness() << "][";
	best->mod->OutputGamlFormattedModel(outf);
	outf << "]";

	outf.setf( ios::floatfield, ios::fixed );
	outf.setf( ios::showpoint );
	best->treeStruct->root->MakeNewick(treeString, false);
	outf << treeString << ";\n";
	outf << "end;\n";
	
	//add a paup block setting the model params
	best->mod->OutputPaupBlockForModel(outf, treefname);
	outf.close();
	
	if(outputPhylipTree){//output a phylip formatted tree if desired
		char phyname[85];
		sprintf(phyname, "%s.phy", treefname);
		ofstream phytree(phyname);
		phytree.precision(8);
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
				phytree << params->data->TaxonLabel(atoi(temp.c_str())-1);
				temp="";
				}
			}
		phytree << ";";
		phytree.close();
		}
	
#else
//if using the UD serial version, just output the best tree in phylip format, with it's score before it
	best->treeStruct->root->MakeNewick(treeString, false);
	outf << best->treeStruct->lnL << "\t" << best->kappa << "\t" << treeString << ";";
	outf.close();
#endif
}

char * Population::MakeNewick(int i, bool internalNodes)
{
	indiv[i].treeStruct->root->MakeNewick(treeString, internalNodes);
	assert(!treeString[stringSize]);
	return treeString;
}
	
void Population::CompactTopologiesList(){
	for(int i=0;i<ntopos;i++)
		{if(topologies[i]->nInds==0)
				{topologies[i]->exNNItried=false;
				TopologyList *temp=topologies[i];
				for(int j=i;j<ntopos-1;j++)
					{topologies[j]=topologies[j+1];
					topologies[j]->DecrementTopoFieldOfInds();
					}
				topologies[ntopos-1]=temp;
				temp->gensAlive=0;
				ntopos--;
				i--;
				}
		}
}

//DZ 7-7 This function will get rid of multiple references to the same treeStruct
//from different individuals.  This keeps double deletion from occuring in the destructor.
//Not the most elegant, but it works.
void Population::EliminateDuplicateTreeReferences(){

	bool dupe;
	vector<Tree *> tstructs;
	
	//go through the indiv array
	for(int i=0;i<params->nindivs;i++){
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
	for(int i=0;i<params->nindivs;i++){
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

/*	this method extends the current population by n using the tree strings in tree_strs to build the
	new individuals.

	side affects:	members indiv and newindiv arrays are extended and new entries filled in.
					member topologies array is extended and new entries filled in.
					member nindiv is changed to reflect the new bigger indiv/newindiv arrays.
					member old_nindiv is changed to hold what nindiv previously held.
					member ntopos is updated to reflect the changes.

*/
/*
int Population::ExtendPopulation(int delta, char* tree_strs, double* kappa_probs)	{

	// sanity check
	assert(current_size == params->nindivs);

	EliminateDuplicateTreeReferences();	// cleans up ptr mess between indiv and newindiv
	
	prResizeIndividualArray(delta, tree_strs, kappa_probs);
	prResizeNewIndividualArray(delta);
	prResizeTopologyListArray(delta);
	prResizeCumFitArray(delta);

	current_size += delta;
	total_size = params->nindivs = current_size;

	return 0;
}
*/
/*
int Population::ShrinkPopulation(int delta, char** tree_strings_, double** kappa_probs_)	{
	char*& tree_strings = *tree_strings_;
	double*& kappa_probs = *kappa_probs_;

	EliminateDuplicateTreeReferences();	// // cleans up ptr mess between indiv and newindiv

	int* to_remove;
	GetNRandomIndivIndices(&to_remove, delta);
	GetSpecifiedTreeStrings(&tree_strings, delta, to_remove);
	GetSpecifiedKappas(&kappa_probs, delta, to_remove);
	delete [] to_remove;

	// TODO the follow funcs don't shrink the population by the indivs in to_remove[] (see above)
	prResizeIndividualArray(-delta);
	prResizeNewIndividualArray(-delta);
	prResizeTopologyListArray(-delta);
	prResizeCumFitArray(-delta);

	current_size -= delta;
	params->nindivs = current_size;
	 
	return 0;
}*/	
/*
int Population::prResizeIndividualArray(int n, char* tree_strs, double* kappa_probs)	{
	int new_size = current_size + n;
	assert(new_size >= 0);
	assert(current_size >= 0);
		
	if (new_size == current_size)	{
		return 0;
	}
	else if (new_size == 0)	{
		assert(indiv);
		for (int i = 0; i < current_size; ++i)
			indiv[i].treeStruct->RemoveTreeFromAllClas();
		delete [] indiv;
		indiv = NULL;
	}
	else if (new_size > current_size)	{
		assert(tree_strs != NULL);
		assert(kappa_probs != NULL);
		Individual* temp_indiv = new Individual[new_size];
		for (int i = 0; i < current_size; ++i)	{
			temp_indiv[i] = indiv[i];
			indiv[i].treeStruct = NULL; // so ~Individual won't destroy its tree
		}
		char* p = tree_strs;
		int j = 0;
		for (int i = current_size; i < new_size; ++i)	{
			temp_indiv[i].treeStruct = new Tree(p, params->data, sharedcl);
			temp_indiv[i].treeStruct->AssignCLAsFromMaster();
			temp_indiv[i].mod->CopyModel(mod);
			temp_indiv[i].params = params;
			p += strlen(p) + 1;
		}
		if (current_size > 0)	{
			assert(indiv);
			delete [] indiv;
		}
		indiv = temp_indiv;
	}
	else if (new_size < current_size)	{
		assert(tree_strs == NULL);
		Individual* temp_indiv = new Individual[new_size];
		for (int i = 0; i < new_size; ++i)	{
			temp_indiv[i] = indiv[i];
			indiv[i].treeStruct = NULL;
		}
		for (int i = new_size; i < current_size; ++i)
			indiv[i].treeStruct->RemoveTreeFromAllClas();
		delete [] indiv;
		indiv = temp_indiv;
	}
	else	{
		// shouldn't get here
		assert(false);
	}

	// sanity check:  make sure each indiv's params ptr is pointing to the right place.
	// sanity check:  make sure each indiv's treeStruct is not NULL.
	for (int i = 0; i < new_size; i++)	{
		assert(indiv[i].params == params);
		assert(indiv[i].treeStruct != NULL);
	}

	return 0;
}

int Population::prResizeNewIndividualArray(int n)	{
	int new_size = current_size + n;
	assert(new_size >= 0);
	assert(current_size >= 0);

	if (current_size == new_size)
		return 0;
	if (current_size > 0)	{
		assert(newindiv);
		delete [] newindiv;
	}
	else
		assert(newindiv == NULL);
	if (new_size > 0)
		newindiv = new Individual[new_size];
	else
		newindiv = NULL;

	return 0;
}
*/
/*	precondition:	indiv[] must be of size "current_size".
*/
/*
int Population::prResizeTopologyListArray(int n)	{
	int new_size = current_size + n;
	assert(new_size > 0);

	for (int i = 0; i < current_size; ++i)
		delete topologies[i];
	delete [] topologies;
	TopologyList::ntoposexamined = 0;

	topologies = new TopologyList*[new_size];
	TopologyList::SetIndL(indiv);
	for(int i = 0; i < new_size; i++)	{
		topologies[i] = new TopologyList();
		topologies[i]->Allocate( (100 < new_size+1) ? 100 : (new_size+1) );
		topologies[i]->AddInd(i);
		TopologyList::ntoposexamined++;
		indiv[i].topo=i;
	}
	
	ntopos = new_size;

	return 0;
}
*/

/*	this is a reincarnation of UpdateTopologyList().
	precondition:  indivs[] must be of size new_size
	question:  topology list grows only and never shrinks?  it seems like ntopos just gets bigger and bigger.
*/
int Population::prResizeTopologyListArray(int n)	{
	int new_size = current_size + n;
	int new_topos=0;
	assert(current_size >= 0);
	assert(new_size >= 0);
	for(int i=0; i < new_size; i++){
		if(indiv[i].topo < 0)	{
			indiv[i].topo=ntopos+new_topos++;
		}
	}
	ntopos+=new_topos;
	TopologyList::SetIndL(indiv);
	for(int i=0;i<ntopos;i++)
		topologies[i]->Clear();
	for(int i = 0; i < new_size; i++ )
		topologies[indiv[i].topo]->AddInd(i);
	CompactTopologiesList();

	return 0;
}

int Population::prResizeCumFitArray(int n)	{
	int new_size = current_size + n;
	assert(current_size >= 0);
	assert(new_size >= 0);

	if (current_size == new_size)
		return 0;
	if (current_size > 0)	{
		assert(cumfit);
		for (int i = 0; i < current_size; ++i)
			delete [] cumfit[i];
		delete [] cumfit;
	}
	else
		assert(cumfit == NULL);
	if (new_size > 0)	{
		cumfit = new double*[new_size];
		for (int i = 0; i < new_size; ++i)
			cumfit[i] = new double[2];
	}
	else
		cumfit = NULL;

	return 0;
}

void Population::CheckAllTrees(){//debugging function
	for(int i=0;i<indiv[0].params->nindivs;i++){
		//check that trees are properly formed
		indiv[i].treeStruct->root->CheckforLeftandRight();
		indiv[i].treeStruct->root->CheckforPolytomies();
		indiv[i].treeStruct->root->CheckTreeFormation();
		//check that no individuals point to the same treeStruct
		for(int j=i+1;j<indiv[0].params->nindivs;j++)
			assert(!(indiv[i].treeStruct==indiv[j].treeStruct));
			}
		}
	
void Population::UpdateTopologyList(Individual *inds){
	//bring topo list up to date
	//also checks if any individuals have not been assigned a topo (topo=-1), and does so.
	//BE SURE that params->nindivs has been updated to the new size before calling this.
	int new_topos=0;
	for(int i=0;i<total_size;i++){
		if(inds[i].topo<0){
			inds[i].topo=ntopos+new_topos++;
			}
		}
	ntopos+=new_topos;
	assert(ntopos<=total_size);
	TopologyList::SetIndL(inds);
	for(int i=0;i<ntopos;i++)
		topologies[i]->Clear();
	for(int i = 0; i <total_size; i++ )
		topologies[inds[i].topo]->AddInd(i);
	CompactTopologiesList();
	if(ntopos < total_size) assert(topologies[ntopos]->nInds == 0);
	}	

void Population::RemoveFromTopologyList(Individual *ind){
	topologies[ind->topo]->nInds--;
	ntopos--;
	ind->topo=-1;
	}			

void Population::CheckIndividuals(){
	for(int i=0;i<params->nindivs;i++){
		assert(!(indiv[i].topo>ntopos));
		}
	}
	
void Population::TopologyReport(){
	//this is for debugging purposes
	ofstream out("toporeport.log");
	if(!(out.good())){
		cout << "problem opening toporeport.log" << endl;
		exit(0);
		}
	for(int i=0;i<total_size;i++){
		out << "topo# " << i << "\t" << "ngen " << topologies[i]->gensAlive << "\t";
		out << "nInds " << topologies[i]->nInds << "\t" << "inds" << "\t";
		for(int j=0;j<topologies[i]->nInds;j++){
			out << topologies[i]->GetIndNums(j) << "\t";
			}
		out << endl;
		}
	out << "ind\ttopo" << endl;
	for(int i=0;i<params->nindivs;i++){
		out << i << "\t" << indiv[i].topo << endl;
		}
	out.close();
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
int Population::SwapIndividuals(int n, const char* tree_strings_in, double* kappa_probs_in, char** tree_strings_out_, double** kappa_probs_out_)	{
	char*& tree_strings_out = *tree_strings_out_;
	double*& kappa_probs_out = *kappa_probs_out_;
	
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
int Population::ReplaceSpecifiedIndividuals(int count, int* which_array, const char* tree_strings, double* model_string)	{
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
		ind->treeStruct = new Tree(tree_strings);
		ind->treeStruct->AssignCLAsFromMaster();
		ind->mod->SetModel(model_string);
		ind->treeStruct->mod=ind->mod;

		ind->SetDirty();
		tree_strings += strlen(tree_strings)+1;
		ind->treeStruct->mod=ind->mod;
		}
	CompactTopologiesList();
	UpdateTopologyList(indiv);
	return 0;
}

int Population::GetNRandomIndivIndices(int** indiv_list, int n)	{
	int* ar = new int[current_size];
	for (int i = 0; i < current_size; ++i)
		ar[i] = i;
	ScrambleArray<int>(current_size, ar);
	*indiv_list = new int[n];
	for (int i = 0; i < n; ++i)
		(*indiv_list)[i] = ar[i];
	delete [] ar;
	return 0;
}

int Population::GetNBestIndivIndices(int** indiv_list, int n)	{
	*indiv_list = new int[n];
	for (int i = 0; i < n; ++i)
		(*indiv_list)[i] = (int)cumfit[current_size-i-1][0];
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

int Population::GetSpecifiedModels(double** model_string, int n, int* indiv_list){
	double *&model = *model_string;
	int string_size=0;
	//first calculate the appropriate size of the string and allocate it
	int nrates=indiv[indiv_list[0]].mod->Nst()-1;
	string_size+=n*nrates;
	string_size+=n*4;//the pi's
	if(indiv[indiv_list[0]].mod->NRateCats()>1) string_size+=1*n;
	if(indiv[indiv_list[0]].mod->ProportionInvariant()!=0.0) string_size+=1*n;
	model=new double[string_size];
	
	int slot=0;
	for (int i = 0; i < n; ++i){
		//get the rates
		for(int r=0;r<nrates;r++)
			model[slot++] = indiv[indiv_list[i]].mod->Rates(r);	
		
		//get the pi's
		for(int b=0;b<4;b++)
			model[slot++] = indiv[indiv_list[i]].mod->Pi(b);
		
		//get alpha if we are using rate het
		if(indiv[indiv_list[0]].mod->NRateCats()>1)
			model[slot++] = indiv[indiv_list[i]].mod->Alpha();
			
		//get pinv if we are using invariant sites
		if(indiv[indiv_list[0]].mod->ProportionInvariant()!=0.0)
			model[slot++] = indiv[indiv_list[i]].mod->ProportionInvariant();		
		}
	return slot;
	}

void Population::OutputLog()	{
	//log << gen << "\t" << bestFitness << "\t" << stopwatch.SplitTime() << "\t" << adap->branchOptPrecision << endl;
	if(gen > -1)
		log << gen << "\t" << indiv[bestIndiv].Fitness() << "\t" << stopwatch.SplitTime() << "\t" << adap->branchOptPrecision << endl;
	else
		log << "Final\t" << indiv[bestIndiv].Fitness() << "\t" << stopwatch.SplitTime() << "\t" << adap->branchOptPrecision << endl;
	}

int Population::ReplicateSpecifiedIndividuals(int count, int* which, const char* tree_string, double *model_string){
	assert(count > 0 && count <= total_size);
	for (int i = 0; i < count; ++i)	{
		indiv[which[i]].treeStruct->RemoveTreeFromAllClas();
		delete indiv[which[i]].treeStruct;
		indiv[which[i]].treeStruct = new Tree(tree_string);
		indiv[which[i]].treeStruct->AssignCLAsFromMaster();
		indiv[which[i]].mod->SetModel(model_string);
		indiv[which[i]].treeStruct->mod=indiv[which[i]].mod;
		indiv[which[i]].SetDirty();
		indiv[which[i]].treeStruct->mod=indiv[which[i]].mod;
		}
	return 0;
}

void Population::UpdateTreeModels(){
	for(int ind=0;ind<current_size;ind++){
		newindiv[ind].treeStruct->mod=newindiv[ind].mod;
//		indiv[ind].treeStruct->mod=indiv[ind].mod;
		}
	}

double Population::IndivFitness(int i) {
	return indiv[i].Fitness();
	}

void Population::OutputModelAddresses(){
	ofstream mods("modeldeb.log", ios::app);
	
	for(int i=0;i<total_size;i++){
		mods << "indiv " << i << "\t" << indiv[i].mod << "\t" << indiv[i].treeStruct->mod << "\n";
		mods << "newindiv " << i << "\t" << newindiv[i].mod << "\t" << newindiv[i].treeStruct->mod << "\n";
		}
	mods << endl;
	}

/*  Methods added by Alan Zhang on Jan,09,2004     */

void Population::NNIoptimization(){
	for(int i=0;i<params->nindivs;i++)
		indiv[i].ResetIndiv();

	//	for(int i = 0;i<params->nindivs;i++){
		bool topoChange=NNIoptimization(bestIndiv, 1);
	//	}
	
	if(topoChange==true){
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
	
	CalcAverageFitness();
}

bool Population::NNIoptimization(int indivIndex, int steps){
	Individual  currentBest;
	Individual  tempIndiv1, tempIndiv2, *best;
	int beginNode, endNode, optiNode;
	double bestNNIFitness; 
	double startingFitness;
	bool betterScore=false;
	
//	ofstream outf("nnidebug.tre");
//	ofstream scr("nniscores.tre");
	ofstream out;
	
	beginNode = newindiv[indivIndex].treeStruct->getNumTipsTotal() + 1;
	endNode = beginNode * 2 - 5;
	startingFitness = indiv[newindiv[indivIndex].parent].Fitness();
	bestNNIFitness = -1e100;
	
	steps = min(max(0,steps),newindiv[indivIndex].treeStruct->getNumTipsTotal()-3);
	indivIndex = min(max(0,indivIndex),params->nindivs-1);

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
			
			double improvement = 0.01;
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
	newindiv[indivIndex].mutation_type |= Individual::exNNI;

	return betterScore;
}

/* End of methods added by Yufeng Zhang*/

void Population::NNIPerturbation(int sourceInd, int indivIndex){
	Individual  currentBest;
	Individual  tempIndiv1, tempIndiv2, *best;
	int optiNode;
	double previousFitness; 
//	bool betterScore=false;
	double scorediff=0.0;
	double thresh=pertMan->nniAcceptThresh;
	int nummoves=0;
	
	ofstream out;
	
//	int numNodes=indiv[indivIndex].treeStruct->getNumTipsTotal()-3;
//	int *nodeArray=new int[numNodes];
/*	
	for(int i=0;i<numNodes;i++){
		nodeArray[i]=indiv[indivIndex].treeStruct->GetRandomInternalNode();
		//DEBUG get all of the nodes, in order
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
	for(int i=0;i<total_size;i++)
		indiv[i].treeStruct->MakeAllNodesDirty();

	char filename[50];
	if(rank < 10)
		sprintf(filename, "pertreport0%d.log", rank);
	else 
		sprintf(filename, "pertreport%d.log", rank);
	ofstream pert(filename, ios::app);
	pert.precision(10);
	pert << "gen\t" << gen << "\tstart\t" << indiv[bestIndiv].Fitness() << "\n";

	#ifndef NO_OUTPUT
	cout << "Performing NNI Perturbation.  Starting score=" << indiv[bestIndiv].Fitness() << endl;
	#endif
	
	//DEBUG
/*	char filename[50];
	double localprec=.5;
	sprintf(filename, "%d.%.4fscores.log", gen, localprec);
	ofstream temp(filename);
	temp.precision(12);
	temp << "start\t" << indiv[bestIndiv].Fitness() << "\n";
*/

//	for(int i=0;i<numNodes;i++){
	int attempts, accepts;
	for(accepts=0, attempts=0;(accepts<pertMan->nniTargetAccepts) && (attempts <= pertMan->nniMaxAttempts);){
		#ifndef NO_OUTPUT
		if(! (attempts++ % (pertMan->nniMaxAttempts/20))) cout << ".";
		cout.flush();
		#endif

		optiNode=indiv[indivIndex].treeStruct->GetRandomInternalNode();
//		optiNode=nodeArray[i];
		//DEBUG
//		tempIndiv1.treeStruct->NNIMutate(optiNode,0, localprec, 0);
//		tempIndiv2.treeStruct->NNIMutate(optiNode,1, localprec, 0);

		tempIndiv1.treeStruct->NNIMutate(optiNode,0, adap->branchOptPrecision, 0);
		tempIndiv2.treeStruct->NNIMutate(optiNode,1, adap->branchOptPrecision, 0);

//		tempIndiv1.SetDirty();
//		tempIndiv2.SetDirty();
		
		tempIndiv1.SetFitness(tempIndiv1.treeStruct->lnL);
		tempIndiv2.SetFitness(tempIndiv2.treeStruct->lnL);

		//DEBUG
//		temp << tempIndiv1.Fitness() << "\n" << tempIndiv2.Fitness() << "\n";

		double diff1=tempIndiv1.Fitness() - previousFitness;
		double diff2=tempIndiv2.Fitness() - previousFitness;

		//ignore NNI's that improve the fitness, because they are probably just undoing a previous NNI
//		if(((diff1 < 0.0) && (diff1 + thresh > 0.0)) || ((diff2 < 0.0) && (diff2 + thresh > 0.0))){
//			if((diff1 < 0.0) && ((diff1 > diff2) || (diff2 >= 0.0))) best=&tempIndiv1;
		if(diff1 < 0.0 || diff2 < 0.0){
			if((diff1 < 0.0) && ((diff1 > diff2) || (diff2 >= 0.0))) best=&tempIndiv1;
			else best=&tempIndiv2;
			
			double acceptanceProb=exp(-params->selectionIntensity * (previousFitness - best->Fitness()));
			if(rnd.uniform() < acceptanceProb){
				double thisdiff=best->Fitness() - previousFitness;
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
	#ifndef NO_OUTPUT
	cout << "\nCompleted Perturbation.\n  " << accepts << " NNI's accepted in " << attempts << " attempts.  Current score=" << indiv[bestIndiv].Fitness() << endl;
	#endif
//	delete []nodeArray;
}

void Population::NNISpectrum(int sourceInd){
	Individual  tempIndiv1, tempIndiv2;
	int optiNode;
	double previousFitness; 
	double scorediff=0.0;
	double thresh=pertMan->nniAcceptThresh;

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
	for(int i=0;i<total_size;i++)
		indiv[i].treeStruct->MakeAllNodesDirty();

	tempGlobal=1;
	OutputFilesForScoreDebugging(&tempIndiv1, tempGlobal++);
	double localprec;
	double prec[7]={.5, .25, .1, .05, .01, .005, .001};
	for(int q=0;q<7;q++){
		localprec=prec[q];

		char filename[50];
		sprintf(filename, "%d.%.4fscores.log", gen, localprec);
		ofstream temp(filename);
		temp.precision(12);
		temp << "start\t" << indiv[bestIndiv].Fitness() << "\n";


		for(int i=0;i<numNodes;i++){

			optiNode=nodeArray[i];
			//DEBUG
			tempIndiv1.treeStruct->NNIMutate(optiNode,0, localprec, 0);
			tempIndiv2.treeStruct->NNIMutate(optiNode,1, localprec, 0);

	//		tempIndiv1.SetDirty();
	//		tempIndiv2.SetDirty();
			
			tempIndiv1.SetFitness(tempIndiv1.treeStruct->lnL);
			tempIndiv2.SetFitness(tempIndiv2.treeStruct->lnL);

			//DEBUG
			temp << tempIndiv1.Fitness() << "\n" << tempIndiv2.Fitness() << "\n";
			if(q==0){
				OutputFilesForScoreDebugging(&tempIndiv1, tempGlobal++);
				OutputFilesForScoreDebugging(&tempIndiv2, tempGlobal++);
				}
			double diff1=tempIndiv1.Fitness() - previousFitness;
			double diff2=tempIndiv2.Fitness() - previousFitness;

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

void Population::SPRoptimization(int indivIndex){
//	for(int i=0;i<params->nindivs;i++)
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
       
}

bool Population::SPRoptimization(int indivIndex, int range, int cutnum ){
	//DJZ 1/23/04 this is based on the NNI optimization function by Alan (ie, exhaustive nnis)
	//the main difference is that only one node will be used as the node to be cut off and
	//reattached, but then all reattachment points within a radius will be tried.
	//the marking of nodes as dirty is also necessarily different
	subset sprRange;
	
	Individual  currentBest;
	Individual  tempIndiv1;
	double bestSPRFitness; 
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
		if(tempIndiv1.Fitness() > bestSPRFitness)
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

void Population::SPRPerturbation(int sourceInd, int indivIndex){
	Individual  currentBest;
	Individual  tempIndiv1;
	Individual *source=&indiv[sourceInd];
	int range=pertMan->sprPertRange;
	double thresh=10000.0;

	
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
		double previousFitness=source->Fitness();
		double bestDiff=-thresh;

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
			#ifndef NO_OUTPUT
			if(! (broken % (int)ceil((double)sprRange->total/5))) cout << ".";
			cout.flush();
			#endif
			tempIndiv1.treeStruct->SPRMutate(cutnum, sprRange->element[broken--], adap->branchOptPrecision, 0, 0);
			if(sprRange->element[broken]==0 && broken>0) broken--;

			tempIndiv1.SetFitness(tempIndiv1.treeStruct->lnL);;

			//divide the score difference by the square root of the node distance, to favor longer moves
			double diff=(tempIndiv1.Fitness()-previousFitness)/sqrt(sprRange->front[broken]+1.0);
			if(diff>0) diff=(tempIndiv1.Fitness()-previousFitness)*(sprRange->front[broken]+1);
	//		pert2 << "node=\t" << sprRange->element[broken] << "\tdist=\t" << sprRange->front[broken] << "\tscore=\t" << tempIndiv1.Fitness() << "\t" << diff << "\t" << tempIndiv1.Fitness()-previousFitness << "\n";
			if(diff > bestDiff/* || diff > thresh*/){
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
		
		#ifndef NO_OUTPUT
		cout << "Accepted SPR with range of " << bestDist << ". Current score=" << indiv[indivIndex].Fitness() << endl;
		#endif
		}

	
	//Return the treestructs that we used temporarily back to the unused tree vector
	tempIndiv1.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(tempIndiv1.treeStruct);
	tempIndiv1.treeStruct=NULL;
	currentBest.treeStruct->RemoveTreeFromAllClas();
	unusedTrees.push_back(currentBest.treeStruct);
	currentBest.treeStruct=NULL;
	}

void Population::CheckPerturbSerial(){

	if(pertMan->pertType < 3 ){
	 	if(pertMan->pertType==1 && (gen - pertMan->lastPertGeneration) >= pertMan->minPertInterval/2 
	 		&& adap->randNNIweight != adap->origRandNNIweight){
			adap->randNNIweight=adap->origRandNNIweight;
//			pertMan->lastPertGeneration=gen;
			}


		if(pertMan->pertAbandoned==false && (gen - pertMan->lastPertGeneration) >= pertMan->minPertInterval 
			&& (adap->improveOverStoredIntervals < pertMan->pertThresh) /*&& (adap->branchOptPrecision == adap->minOptPrecision)*/){
			if(pertMan->numPertsNoImprove <= pertMan->maxPertsNoImprove){
				if(indiv[bestIndiv].Fitness() > bestSinceRestart.Fitness()){
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
					//DEBUG try disallowing NNIs immediately after the perturbation
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
					double startscore=indiv[bestIndiv].Fitness();
					double curscore;
//					do{
					int indToReplace = (bestIndiv==0 ? 1 : 0);
					int source=bestIndiv;
					#ifndef NO_OUTPUT
					cout << "Performing SPR Perturbation.  Starting score=" << indiv[bestIndiv].Fitness() << endl;
					#endif
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
				if(indiv[bestIndiv].Fitness() > bestSinceRestart.Fitness()){
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
				bestFitness=indiv[bestIndiv].Fitness();
				pertMan->lastPertGeneration=gen;
				pertMan->scoreAfterRatchet=indiv[bestIndiv].Fitness();
				adap->reset=true;
				gen++;
				OutputFate();
				cout << "Performing ratcheting: reweighting " << pertMan->ratchetProportion*100 << " percent of characters." << endl; 
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
			if((gen - pertMan->lastPertGeneration) >= pertMan->ratchetMaxGen || indiv[bestIndiv].Fitness() - pertMan->scoreAfterRatchet > pertMan->ratchetOffThresh){
				TurnOffRatchet();
				gen++;
				OutputFate();
				}
			}
		}
	}

void Population::TurnOffRatchet(){
	params->data->RestoreOriginalCounts();
	pertMan->ratcheted=false;
	
	claMan->MakeAllHoldersDirty();
	for(int i=0;i<total_size;i++) indiv[i].SetDirty();
	CalcAverageFitness();
	bestFitness=indiv[bestIndiv].Fitness();
	pertMan->lastPertGeneration=gen;
	adap->reset=true;
	cout << "Returning to normal character weighting..." << endl;
	char filename[50];
	if(rank < 10)
		sprintf(filename, "pertreport0%d.log", rank);
	else 
		sprintf(filename, "pertreport%d.log", rank);
	ofstream pert(filename, ios::app);
	pert << "Returning to normal character weighting..." << endl;
	pert.close();
	}


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
		&& (adap->improveOverStoredIntervals < pertMan->pertThresh)/* && (adap->branchOptPrecision == adap->minOptPrecision)*/){
		for(int i=1;i<=paraMan->nremotes;i++){
			paraMan->needToSend[i]=true;
			}
		paraMan->perturbModeActive=true;
		paraMan->allSent=false;
		pertMan->lastPertGeneration=gen;
		}	
	}

void Population::RestoreAllTimeBest(){
	UpdateTopologyList(indiv);
	topologies[indiv[0].topo]->RemoveInd(0);
	CompactTopologiesList();
	indiv[0].CopySecByRearrangingNodesOfFirst(indiv[0].treeStruct, &allTimeBest, true);
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
	indiv[0].CopySecByRearrangingNodesOfFirst(indiv[0].treeStruct, &bestSinceRestart, true);
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
	#ifndef UNIX
	cout << "restoring best individual with score of " << bestSinceRestart.Fitness() << ".\n" << pertMan->numPertsNoImprove << " perturbation(s) performed without improvement." << endl;
	#endif
	}

void Population::StoreBestForPert(){
	if(indiv[bestIndiv].Fitness() > allTimeBest.Fitness()) StoreAllTimeBest();
	
	if(bestSinceRestart.treeStruct==NULL){
		if(unusedTrees.empty()){
			Tree *temp=new Tree();
			unusedTrees.push_back(temp);
			}					
		bestSinceRestart.treeStruct=*(unusedTrees.end()-1);
		unusedTrees.pop_back();
		}
	bestSinceRestart.CopySecByRearrangingNodesOfFirst(bestSinceRestart.treeStruct, &indiv[bestIndiv]);
	bestSinceRestart.topo=-1;
	//need to do this to be sure that the bestSinceRestart isn't tying up clas
	bestSinceRestart.treeStruct->RemoveTreeFromAllClas();
	
	char filename[50];
	if(rank < 10)
		sprintf(filename, "pertreport0%d.log", rank);
	else 
		sprintf(filename, "pertreport%d.log", rank);
	ofstream pert(filename, ios::app);
	pert.precision(10);
	pert << "storing best individual with score of " << bestSinceRestart.Fitness() << "\n";
	pert.close();
	#ifndef UNIX
	cout << "storing best individual with score of " << bestSinceRestart.Fitness() << endl;
	#endif	
	}

void Population::StoreAllTimeBest(){
	CreateTreeFile(params->treefname);
	if(allTimeBest.treeStruct==NULL){
		if(unusedTrees.empty()){
			Tree *temp=new Tree();
			unusedTrees.push_back(temp);
			}					
		allTimeBest.treeStruct=*(unusedTrees.end()-1);
		unusedTrees.pop_back();
		}
	allTimeBest.CopySecByRearrangingNodesOfFirst(allTimeBest.treeStruct, &indiv[bestIndiv]);
	allTimeBest.topo=-1;
	//need to do this to be sure that the alltimebest isn't tying up clas
	allTimeBest.treeStruct->RemoveTreeFromAllClas();
	}

void Population::keepTrack(){

	if(((gen-1)%adap->intervalLength)==0){
		if(gen>1) adap->PrepareForNextInterval();
		}

	//remember that the indiv and newindiv arrays have already been switched, so the newindivs are the parents of the indivs

	//adap->reset will be true if we've Ratcheted, in which case the scores will be noncomparable
	//so reset these values
	if(adap->reset==true){
		adap->laststepscore=indiv[bestIndiv].Fitness();
		adap->lastgenscore=indiv[bestIndiv].Fitness();
		adap->reset=false;
		}
	
	if(gen==1)
		adap->lastgenscore = adap->laststepscore = newindiv[0].Fitness();
		
	adap->improvetotal[0] = indiv[bestIndiv].Fitness() - adap->laststepscore;

	for(int i=0;i<params->nindivs;i++){
//		double scoreDif=indiv[i].Fitness() - newindiv[indiv[i].parent].Fitness();
		double scoreDif=indiv[i].Fitness() - adap->lastgenscore;
		int typ = indiv[i].mutation_type;
		if(typ > 0){
			if(scoreDif>0){
	//			if(i==bestIndiv){
					//keep track of when the last significant beneficial topo mutation occured
					//this will be used for the stopping criterion, precision reduction and update reduction in the parallel version
					if(typ&Individual::anyTopo){
						if(scoreDif > significantTopoChange){
							indiv[0].treeStruct->CalcBipartitions();
							indiv[i].treeStruct->CalcBipartitions();
							if(indiv[0].treeStruct->IdenticalTopology(indiv[i].treeStruct->root)==false){
								lastTopoImprove=gen;
								if(i == bestIndiv) AppendTreeToTreeLog(indiv[bestIndiv].mutation_type);
								}
							}
						}
					
					if(typ&(Individual::randNNI)){
						adap->randNNI[0] += scoreDif;
						}
					if(typ&(Individual::exNNI)) 		adap->exNNI[0] += scoreDif;
					if(typ&(Individual::randSPR)) 	adap->randSPR[0] += scoreDif;
					if(typ&(Individual::limSPR)) 		adap->limSPR[0] += scoreDif;
					if(typ&(Individual::exlimSPR)) 	 	adap->exlimSPR[0] += scoreDif;
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
			if(typ&(Individual::exNNI)) 	 	adap->exNNInum[0]++;
			if(typ&(Individual::randSPR)) 	adap->randSPRnum[0]++;
			if(typ&(Individual::limSPR))	 	adap->limSPRnum[0]++;
			if(typ&(Individual::exlimSPR)) 	adap->exlimSPRnum[0]++;
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

	adap->lastgenscore=indiv[bestIndiv].Fitness();

	//things to do on the final generation of an interval	
	if(gen%adap->intervalLength==0){
		//improveOverStoredIntervals is only used on generations that are multiples of intervalLength
		//so it won't contain the improvement in the latest interval until it's end
		adap->improveOverStoredIntervals=0.0;
		for(int i=0;i<adap->intervalsToStore;i++)
			adap->improveOverStoredIntervals += adap->improvetotal[i];
		if(adap->improveOverStoredIntervals < 0.0) adap->improveOverStoredIntervals = 0.0;
		//update the mutation probailities
		if(gen>=(adap->intervalLength*adap->intervalsToStore*0.5)){
			adap->UpdateProbs();
			if(outputMostlyUselessFiles) adap->OutputProbs(probLog, gen);
			}
		adap->laststepscore=indiv[bestIndiv].Fitness();
		}
	}

int ParallelManager::DetermineSubtrees(Tree *tr, ofstream &scr){
	//Determine what the best node we could choose to be the root would be in terms of 
	//the partitioning efficiency

	TreeNode *nd=tr->root;

	int bestRoot=0, orphans=ntax;
	bool done=false;
	double bestScore=0.0;
	double thisScore;

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
   	bestScore += pow((double)(orphans), orphanFactor);

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
			
			thisScore=0.0;
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
		   	thisScore+=pow((double)(orphans), orphanFactor);
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
	double bestScore=ScorePartitioning(nd->nodeNum, pscores);	
	double thisScore;
	
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
	double origOrphanFactor=paraMan->orphanFactor;

	ofstream pscores("partscores.log", ios::app);

	do{
		pscores << "gen " << gen << " attempt " <<  attempt << endl;
		int bestRoot=paraMan->DetermineSubtrees(indiv[bestIndiv].treeStruct, pscores);	
		//now we need to reroot to the best root found

		pscores << "best root=" << bestRoot << "\n";

		if(bestRoot!=0){
			indiv[bestIndiv].treeStruct->RerootHere(bestRoot);
			AssignNewTopology(indiv, bestIndiv);
			}

		//mark all of the trees inaccurate
		for(int i=0;i<total_size;i++){
			indiv[i].accurateSubtrees=false;
			newindiv[i].accurateSubtrees=false;
			}

		indiv[bestIndiv].accurateSubtrees=true;
		FillPopWithClonesOfBest();

		TreeNode *nd=indiv[bestIndiv].treeStruct->root;

		int orphans=paraMan->ntax;
		bool done=false;
		double bestScore=0.0;

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
	for(int i=0;i<total_size;i++){
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
 		Subtree *st = new Subtree(pointer->nodeNum, n, pointer->dlen, 0.0);
		subtrees.push_back(st);
  		}
	else{
    	if(n1>=minSubtreeSize) Partition(pointer->left);
	    if(n2>=minSubtreeSize) Partition(pointer->right);
  		}
	}

void ParallelManager::NewPartition(TreeNode *pointer, int &orphans, vector<Subtree*> &subtreesAbove){
	vector<Subtree*> subtreesUpLeft, subtreesUpRight;
	double scoreAbove=0.0;
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
	scoreAbove += pow((double)orphans, orphanFactor);
	
	double scoreHere = pow((double)(targetSubtreeSize-n), 2);
	
	if(scoreAbove > scoreHere){
		for(vector<Subtree*>::iterator it = subtreesAbove.begin();it!=subtreesAbove.end();it++){
			Subtree *del=*it;
			delete del;
			}
		subtreesAbove.clear();
		Subtree *poo=new Subtree(pointer->nodeNum, n, pointer->dlen, pow((double)(targetSubtreeSize-n), 2));
		subtreesAbove.push_back(poo);
		orphans=0;
		}
	}

void ParallelManager::NewPartitionDown(TreeNode *pointer, TreeNode *calledFrom, int &orphans, vector<Subtree*> &subtreesAbove){
	int n, n1, n2;
	TreeNode *sib, *anc;
	vector<Subtree*> subtreesUpLeft, subtreesUpRight;
	double scoreAbove=0.0;
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
		scoreAbove += pow((double)orphans, orphanFactor);
		
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
		scoreAbove += pow((double)orphans, orphanFactor);
		}			

	double scoreHere = (targetSubtreeSize-n)*(targetSubtreeSize-n);
	
	if(scoreAbove > scoreHere){
		for(vector<Subtree*>::iterator it = subtreesAbove.begin();it!=subtreesAbove.end();it++){
			delete *it;
			}
		subtreesAbove.clear();
		Subtree *poo=new Subtree(pointer->nodeNum, n, calledFrom->dlen, pow((double)(targetSubtreeSize-n), 2));
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
	 	Subtree *st = new Subtree(pointer->nodeNum, n, calledFrom->dlen, 0.0);
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
	 	Subtree *st = new Subtree(pointer->nodeNum, n, calledFrom->dlen, 0.0);
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
		for(int i=0;i<total_size;i++){
			if(indiv[i].accurateSubtrees==true){
			    if(i<params->nindivs) count++;
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

		if(gen - paraMan->subtreeDefGeneration > paraMan->subtreeInterval && adap->improveOverStoredIntervals < paraMan->subtreeStartThresh){
			paraMan->subtreeModeActive=true;
			paraMan->needUpdate=true;
			}
		}
	
	if(paraMan->subtreeModeActive==true && paraMan->needUpdate==true){
		StartSubtreeMode();
		}

	if(paraMan->subtreeModeActive==true){
		//determine some conditions for stopping subtree mode here
		if(gen - paraMan->subtreeDefGeneration > paraMan->subtreeInterval){
			StopSubtreeMode();
			}
		}
	}

void Population::FillPopWithClonesOfBest(){
	UpdateTopologyList(indiv);
	Individual *best=&indiv[bestIndiv];
	best->treeStruct->mod=best->mod;
	for(int i=0;i<params->nindivs;i++){
		if(&indiv[i]!=best){
			indiv[i].treeStruct->RemoveTreeFromAllClas();
			topologies[indiv[i].topo]->RemoveInd(i);
			indiv[i].CopySecByRearrangingNodesOfFirst(indiv[i].treeStruct,best);
			topologies[indiv[bestIndiv].topo]->AddInd(i);
			indiv[i].mutation_type=-1;
			}
		indiv[i].treeStruct->mod=indiv[i].mod;
		}
	UpdateTopologyList(indiv);
	CalcAverageFitness();
	}

void Population::AssignSubtree(int st, int indNum){
	subtreeNode=st;

	//we'll do all of this stuff if we are assigning a new subtree or if 
	//we are assigning 0 (turning off subtree mode)
	for(int i=0;i<params->nindivs;i++){
	    indiv[i].accurateSubtrees=false;
		newindiv[i].accurateSubtrees=false;
		indiv[i].treeStruct->UnprotectClas();
		}
	ResetMemLevel(params->data->NTax()-2,claMan->NumClas());
	indiv[bestIndiv].treeStruct->ProtectClas();

	subtreeMemberNodes.clear();
	//add all of the nodenums in the subtree into the subtreeMemberNodes
	//note that the subtree node itself is not added
	if(subtreeNode!=0){
		if(rank==0) assert(indiv[indNum].accurateSubtrees==true);
		indiv[indNum].treeStruct->allNodes[subtreeNode]->left->AddNodesToList(subtreeMemberNodes);
	
		sort(subtreeMemberNodes.begin(),subtreeMemberNodes.end());
		reverse(subtreeMemberNodes.begin(),subtreeMemberNodes.end());
		for(int i=0;i<params->nindivs;i++){
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
		for(int i=0;i<params->nindivs;i++){
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
	for(int i=params->nindivs;i<total_size;i++){
		if(indiv[i].accurateSubtrees==true && (paraMan->localSubtreeAssign[i-params->nindivs+1] > 0)){
			paraMan->CheckSubtreeAccuracy(indiv[i].treeStruct);
	       	recomEligable[i]=true;
			if(i>=params->nindivs)count++;
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
		subrec << p << "\t" << paraMan->localSubtreeAssign[p] << "\t" << paraMan->remoteSubtreeAssign[p] << "\t" << indiv[p+params->nindivs-1].Fitness() << "\n";
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
	for(int who=params->nindivs;who<total_size;who++){
#else
	for(int who=0;who<total_size;who++){
#endif
		if(recomEligable[who]==true){

#ifndef FAKE_PARALLEL
			tempIndiv1.treeStruct->SubtreeBasedRecombination(indiv[who].treeStruct, paraMan->localSubtreeAssign[who - params->nindivs + 1], false, adap->branchOptPrecision);
#else
			tempIndiv1.treeStruct->SubtreeBasedRecombination(newindiv[who].treeStruct, subtreeNode , tempIndiv1.mod->IsModelEqual(newindiv[who].mod), adap->branchOptPrecision);
#endif		

			//DEBUG
//			OutputFilesForScoreDebugging(&tempIndiv1, poo++);
		//	paupf.flush();
		//	outf.flush();

			tempIndiv1.SetDirty();
			tempIndiv1.CalcFitness(subtreeNode);
			
			//DEBUG
/*			ofstream poo("debug.log");
			poo.precision(10);
			tempIndiv1.treeStruct->OutputFirstClaAcrossTree(poo, tempIndiv1.treeStruct->root);
			poo.close();
*/			
			subrec << "with " << who << "\t(node " << paraMan->localSubtreeAssign[who - params->nindivs + 1] << ")\t" << tempIndiv1.Fitness() << "\n";
			if(tempIndiv1.Fitness() > currentBest.Fitness()){
				//if the recombinant we create is better, make it the current best, mark it as 
				//ineligable so we don't try to add it again, and start back at the first eligable
				//recominant
				currentBest.CopySecByRearrangingNodesOfFirst(currentBest.treeStruct, &tempIndiv1, true);
				recomEligable[who]=false;
				//who=params->nindivs-1;
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


double ParallelManager::ScorePartitioning(int nodeNum, ofstream &pscores){
	
	int size=(int)subtrees.size();
	
	if(size<2 /*|| size>(nremotes-1)*/) return 9999999999.9;

	double blenScore=0.0, subScore=0.0, fosterScore=0.0;
	int fosterTerms=ntax;
	
	int allots[1024];
	
#ifndef MPI_VERSION
nremotes=9;
#endif
	
	int a=0;
	for(vector<Subtree*>::iterator it=subtrees.begin();it!=subtrees.end();it++){
		blenScore -= log((double)(*it)->blen);
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
	if(fosterTerms> (ntax/20)) fosterScore = (fosterTerms-(ntax/20))*5;
	else fosterScore=0;
	blenScore*=2.0;
	
	double tot= subScore + blenScore + fosterScore + maxallot;
	
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

void Population::InitializeOutputStreams(){
	char temp_buf[30];

	if(outputMostlyUselessFiles){
		//initialize the fate file
		if (rank > 9)
			sprintf(temp_buf, "%s.fate%d.log", params->ofprefix, rank);
		else
			sprintf(temp_buf, "%s.fate0%d.log", params->ofprefix, rank);
		
		fate.open(temp_buf);
		fate.precision(10);
		#ifdef MPI_VERSION
		fate << "gen\tind\tparent\trecomWith\tscore\tMutType\t#brlen\taccurateSubtrees\tTime\tprecision\n";
		#else
		fate << "gen\tind\tparent\tscore\tMutType\t#brlen\tTime\tprecision\n";
		#endif	

		//initialize the problog
		if (rank > 9)
			sprintf(temp_buf, "%s.problog%d.log", params->ofprefix, rank);
		else
			sprintf(temp_buf, "%s.problog0%d.log", params->ofprefix, rank);
		probLog.open(temp_buf);
		if(!probLog.good()) cout << "problem opening problog" << endl;
		adap->BeginProbLog(probLog);
		}

	//initialize the log file
	if (rank > 9)
		sprintf(temp_buf, "%s.log%d.log", params->ofprefix, rank);
	else
		sprintf(temp_buf, "%s.log0%d.log", params->ofprefix, rank);
	log.open(temp_buf);
	log.precision(10);
	log << "random seed = " << rnd.init_seed() << "\n";
	log << "gen\tbest_like\ttime\toptPrecision\n";

	//initialize the treelog
	if(outputTreelog){
		if (rank > 9)
			sprintf(temp_buf, "%s.treelog%d.log", params->ofprefix, rank);
		else
			sprintf(temp_buf, "%s.treelog0%d.log", params->ofprefix, rank);

		treeLog.open(temp_buf);
		treeLog.precision(10);

		params->data->BeginNexusTreesBlock(treeLog);
		}
	
	//initialize the bootstrap tree file
		if(bootstrapReps > 0 && rank==0){
			sprintf(temp_buf, "%s.boot.tre", params->ofprefix);
	
			bootLog.open(temp_buf);
			bootLog.precision(10);

			params->data->BeginNexusTreesBlock(bootLog);
			}	
/*	//initialize the modellog
	if (rank > 9)
		sprintf(temp_buf, "%s.modellog%d.log", params->ofprefix, rank);
	else
		sprintf(temp_buf, "%s.modellog0%d.log", params->ofprefix, rank);
	modelLog.open(temp_buf);
	
	modelLog << "gen\tr1\tr2\tr3\tr4\tr5\tpiA\tpiC\tpiG\tpiT\talpha\tpinv\n";
*/	

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

void Population::SetNewBestIndiv(int indivIndex){
	//this should be called when a new best individual is set outside of
	//CalcAverageFitness.  Particularly important for parallel.
	bestIndiv=indivIndex;
	globalBest=bestFitness=prevBestFitness=indiv[bestIndiv].Fitness();
	for(int i=0;i<total_size;i++){
		if(i != bestIndiv){
			indiv[i].treeStruct->UnprotectClas();
			}
		}
	indiv[bestIndiv].treeStruct->ProtectClas();
	}

void Population::FinalizeOutputStreams(){
	fate.close();
	log.close();
	treeLog << "end;\n";
	treeLog.close();
	if(bootLog.is_open()){
		bootLog << "end;\n";
		bootLog.close();
		}
	modelLog.close();
	probLog.close();
	
	#ifdef DEBUG_SCORES
	outf << "end;\n";
	outf.close();
	paupf << "end;\n";
	paupf.close();
	#endif
	}

void Population::FindLostClas(){
	vector<CondLikeArray *> arr;
	
	for(int i=0;i<total_size;i++){
		Tree *t=indiv[i].treeStruct;
		if(! (claMan->IsDirty(t->allNodes[0]->claIndexDown)))
			arr.push_back(claMan->GetCla(t->allNodes[0]->claIndexDown));
		if(! (claMan->IsDirty(t->allNodes[0]->claIndexUL)))
			arr.push_back(claMan->GetCla(t->allNodes[0]->claIndexUL));			
		if(! (claMan->IsDirty(t->allNodes[0]->claIndexUR)))
			arr.push_back(claMan->GetCla(t->allNodes[0]->claIndexUR));
		for(int n=t->getNumTipsTotal()+1;n<t->getNumNodesTotal();n++){
			if(! (claMan->IsDirty(t->allNodes[n]->claIndexDown)))
				arr.push_back(claMan->GetCla(t->allNodes[n]->claIndexDown));
			if(! (claMan->IsDirty(t->allNodes[n]->claIndexUL)))
				arr.push_back(claMan->GetCla(t->allNodes[n]->claIndexUL));			
			if(! (claMan->IsDirty(t->allNodes[n]->claIndexUR)))
				arr.push_back(claMan->GetCla(t->allNodes[n]->claIndexUR));
			}
		}
	sort(arr.begin(), arr.end());
	for(vector<CondLikeArray*>::iterator vit=arr.begin();vit!=arr.end();){
		if(*(vit)==*(vit+1)) arr.erase(vit+1);
		else vit++;
		if((vit+1)==arr.end()){
			break;
			}
		}
	assert(arr.size() + claMan->NumFreeClas() == claMan->NumClas());
	}

void Population::LogNewBestFromRemote(double scorediff, int ind){
#ifdef MPI_VERSION
        adap->bestFromRemoteNum[0]++;
        adap->bestFromRemote[0]+=scorediff;
        
        indiv[0].treeStruct->CalcBipartitions();
        indiv[ind].treeStruct->CalcBipartitions();
        
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
	double totBestFromRemote=0.0;
	int totBestFromRemoteNum=0;
	for(int i=0;i<adap->intervalsToStore;i++){	
	         totBestFromRemoteNum += adap->bestFromRemoteNum[i];
        	 totBestFromRemote += adap->bestFromRemote[i];
		}
	if(totBestFromRemoteNum==0) paraMan->ReduceUpdateThresh();

#endif
}
