// GARLI version 0.952b2 source code
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

#ifndef POPULATION_H
#define POPULATION_H

#include <iostream>
#include <vector>
#include <cfloat>
#include <cassert>

using namespace std;

#include "configoptions.h"
#include "individual.h"
#include "stopwatch.h"
#include "errorexception.h"

class CondLikeArray;
class TopologyList;
class Tree;
class ClaManager;
class Adaptation;

class Subtree{
	public:
	int nodeNum;
	int taxa;
	FLOAT_TYPE blen;
	int numAssigned;
	int priority;
	FLOAT_TYPE score;
	
	Subtree(int nn, int t, FLOAT_TYPE b, FLOAT_TYPE s){
		nodeNum=nn;
		taxa=t;
		priority=1;
		numAssigned=0;
		blen=b;
		score=s;
		}
	void Log(ofstream &out){
		out << nodeNum << "\t" << taxa << "\t" << score << "\n";
		}
	
	};

class ParallelManager{
	public:
	bool subtreeModeActive;
	bool perturbModeActive;

	int nremotes;
	int ntax;

	//variables related to subtree mode
	int subtreeDefNumber;
	int subtreeDefGeneration;
	FLOAT_TYPE subtreeDefScore;	//the score of the tree that the subtrees were initially calculated on
	int lastFullRecom; 		//the last time we tried a subtree based recombination
	vector<Subtree *> subtrees;
	bool fewNonSubtreeNodes;

	int *remoteSubtreeAssign;  // which node the remote was most recently assigned
	int *localSubtreeAssign;  // which node the remotes _were_ working on for the most recent tree received

	bool needUpdate;			//this means that the subtrees determined are no longer guaranteed to be valid for any trees in the population
	bool beforeFirstSubtree;

	FLOAT_TYPE updateThresh;
	FLOAT_TYPE startUpdateThresh;
	FLOAT_TYPE minUpdateThresh;
	FLOAT_TYPE updateReductionFactor;
	
	int subtreeInterval;
	FLOAT_TYPE subtreeStartThresh;
	int	minSubtreeSize;
	int targetSubtreeSize;
	FLOAT_TYPE orphanFactor;
	
	int maxRecomIndivs;

//	FLOAT_TYPE recalcThresh;		//the score threshold for forcing a recalc of the subtrees
//	FLOAT_TYPE subtreeThresh; 	//the score threshold for sending a remote a new tree if it 
								//is currently working on a subtree
//	FLOAT_TYPE nonSubtreeThresh;	//the score threshold for sending a remote a new tree if subtree
								//mode is off
 	

	bool perturb;
	bool *needToSend;
	bool allSent;

	public:
	vector<int> nonSubtreeNodesforNNI;
	vector<int> nonSubtreeNodesforSPR;

	void FindNonSubtreeNodes(TreeNode *nd);
  	
  	public:
	ParallelManager(int _ntax, int nproc, const MasterGamlConfig *mc){
		ntax=_ntax;
		needUpdate=true;
		subtreeModeActive=false;
		perturbModeActive=false;
		perturb=false;

		subtreeDefGeneration =-1;
		subtreeDefScore=-1;
		subtreeDefNumber=0;
		lastFullRecom=-1;

		beforeFirstSubtree=true;

		nremotes=nproc-1;
		remoteSubtreeAssign=new int[nproc];
		localSubtreeAssign=new int[nproc];
		needToSend=new bool[nproc];

		for(int i=0;i<nproc;i++){
			remoteSubtreeAssign[i]=0;
			localSubtreeAssign[i]=0;
			needToSend[i]=false;
			}
		allSent=true;
		
//		recalcThresh = mc->subtreeRecalcThresh;
		updateThresh = startUpdateThresh = mc->startUpdateThresh;
		minUpdateThresh = mc->minUpdateThresh;
		//DJZ 2/20/06
		//making the reduction factor depend on the min and max updateThreshes
		//and the number of reductions requested for the optprecision to 
		//go from its start to min
		//updateReductionFactor = mc->updateReductionFactor;
		
		updateReductionFactor=pow((FLOAT_TYPE) minUpdateThresh/startUpdateThresh, (FLOAT_TYPE) 1.0/ (mc->numPrecReductions));


//		subtreeThresh = mc->subtreeUpdateThresh;
//		nonSubtreeThresh = mc->nonsubtreeUpdateThresh;
		subtreeInterval = mc->subtreeInterval;
		subtreeStartThresh = mc->subtreeStartThresh;
		minSubtreeSize=mc->minSubtreeSize;
		targetSubtreeSize=mc->targetSubtreeSize;
		orphanFactor=mc->orphanFactor;
		
		maxRecomIndivs=mc->maxRecomIndivs;
		} 
	
	~ParallelManager(){
		if(remoteSubtreeAssign != NULL) delete []remoteSubtreeAssign;
		if(localSubtreeAssign != NULL) delete []localSubtreeAssign;
		if(needToSend != NULL) delete []needToSend;
		}
	
	bool ReadyForSubtreeRecom(int gen){
		return (subtreeModeActive==true && (gen - lastFullRecom >= subtreeInterval/2));
		}

	void ReduceUpdateThresh(){
		updateThresh *= updateReductionFactor;
		if(updateThresh < minUpdateThresh) updateThresh=minUpdateThresh;
		}

	void CheckSubtreeAccuracy(const Tree *tr){
		for(int i=0;i<(int)subtrees.size();i++){
			int countnum= tr->allNodes[subtrees[i]->nodeNum]->CountTerminals(0);
			assert(countnum== subtrees[i]->taxa);
			}	
		}
	void ClearSubtrees(){
		for(int i=0;i<(int)subtrees.size();i++){
			delete subtrees[i];
			}
		subtrees.clear();
		}
	int DetermineSubtrees(Tree *tr, ofstream &);
	void Partition(TreeNode *pointer);
	void NewPartition(TreeNode *pointer, int &orphans, vector<Subtree*> &subtreesAbove);
	void NewPartitionDown(TreeNode *pointer, TreeNode *calledFrom, int &orphans, vector<Subtree*> &subtreesAbove);
	void PartitionDown(TreeNode *pointer, TreeNode *calledFrom);
	FLOAT_TYPE ScorePartitioning(int nodeNum, ofstream &pscores);
	void PrepareForSubtreeMode(Individual *ind, int gen);
	int ChooseSubtree();
			
  };

#ifdef INCLUDE_PERTURBATION
class PerturbManager{
	public:
	int lastPertGeneration;
	bool pertAbandoned;
	int numPertsNoImprove;
	FLOAT_TYPE prevPertScore;
	FLOAT_TYPE scoreAfterRatchet;

	int pertType;
	FLOAT_TYPE pertThresh;
	int minPertInterval;
	int maxPertsNoImprove;
	bool restartAfterAbandon;
	int gensBeforeRestart;

	FLOAT_TYPE ratchetProportion;
	FLOAT_TYPE ratchetOffThresh;
	int ratchetMaxGen;
	bool ratcheted;

	FLOAT_TYPE nniAcceptThresh;
	int nniTargetAccepts;
	int nniMaxAttempts;
	
	int numSprCycles;
	int sprPertRange;
	
	public:
	PerturbManager(){
		pertAbandoned=true;
		ratcheted=false;
		pertThresh=0.0;
		}
	PerturbManager(const GeneralGamlConfig *conf){
		lastPertGeneration=-1;
		pertAbandoned=false;
		numPertsNoImprove=0;
		prevPertScore=-1;
		scoreAfterRatchet=-1;

		pertType = conf->pertType;
		if(pertType!=1 && pertType!=3){
			throw ErrorException("Sorry, only pertTypes 1 and 3 and currently supported!");
			}
		pertThresh = conf->pertThresh;
		minPertInterval = conf->minPertInterval;
		maxPertsNoImprove = conf->maxPertsNoImprove;
		restartAfterAbandon = conf->restartAfterAbandon;
		gensBeforeRestart = conf->gensBeforeRestart;
				
		nniTargetAccepts = conf->nniTargetAccepts;
		nniMaxAttempts = conf->nniMaxAttempts;

		ratchetProportion = conf->ratchetProportion;
		ratchetOffThresh = conf->ratchetOffThresh;
		ratchetMaxGen = conf->ratchetMaxGen;
		ratcheted=false;

		numSprCycles = conf->numSprCycles;
		sprPertRange = conf->sprPertRange;
		}
	};
#endif

class Population{

public:
	GeneralGamlConfig *conf;
	ClaManager *claMan;
	Adaptation *adap;
	unsigned gen;
	Individual* indiv;

private:
	unsigned total_size; //this will be equal to conf->nindiv, except in 
					//the case of the parallel master

	Individual* newindiv;
	FLOAT_TYPE bestFitness;
	int rank;//denotes which processor this is.  0 if serial
	int bestIndiv;
	int bestAccurateIndiv;

	int subtreeNode;
	int subtreeDefNumber;
	vector<int> subtreeMemberNodes;

#ifdef INCLUDE_PERTURBATION
	PerturbManager *pertMan;
#endif
	ParallelManager *paraMan;
	
	//termination related variables
	unsigned lastTopoImprove;
	unsigned lastPrecisionReduction;
	
	//DJZ adding these streams directly to the class so that they can be opened once and left open
	ofstream fate;
	ofstream log;
	ofstream treeLog;
	ofstream probLog;
	ofstream bootLog;
	ofstream bootLogPhylip;
	ofstream swapLog;

	string besttreefile;
	char *treeString;
	int stringSize;

	unsigned ntopos;
	bool prematureTermination;//if the user killed the run
		
	vector<Tree *> unusedTrees;
	FLOAT_TYPE prevBestFitness;
	FLOAT_TYPE avgfit;
	int new_best_found;
		
	int numgensamebest;

	public:
		enum { nomem=1, nofile, baddimen };
		int error;

	FLOAT_TYPE** cumfit;//allocated in setup, deleted in dest
		
	TopologyList **topologies;//allocated in Setup(), deleted in dest
			
	HKYData* data;
	Stopwatch stopwatch;
    FLOAT_TYPE starting_wtime;
    FLOAT_TYPE final_wtime;

#ifdef INCLUDE_PERTURBATION
	Individual *allTimeBest; //this is only used for perturbation or ratcheting
	Individual *bestSinceRestart;
#endif

	public:
		Population() : error(0), conf(NULL),
			bestFitness(-(FLT_MAX)), bestIndiv(0),
			prevBestFitness(-(FLT_MAX)),indiv(NULL), newindiv(NULL),
			cumfit(NULL), new_best_found(0), avgfit(0.0),
			gen(0), starting_wtime(0.0), final_wtime(0.0),
			paraMan(NULL), subtreeDefNumber(0), claMan(NULL), 
			treeString(NULL), adap(NULL),
			topologies(NULL), prematureTermination(false)
#ifdef INCLUDE_PERTURBATION			 
			pertMan(NULL), allTimeBest(NULL), bestSinceRestart(NULL),
#endif
			{
			lastTopoImprove = 0;
			lastPrecisionReduction = 0;
			}

		~Population();

		void QuickSort( FLOAT_TYPE **scoreArray, int top, int bottom );
		FLOAT_TYPE BestFitness() {
//			assert(bestFitness == indiv[bestIndiv].Fitness());
			return bestFitness; 
			}

#if !defined( PARALLEL_MPI_VERSION )
		void Run();
#endif

		int TimeToQuit();
		void CatchInterrupt();
		
		char *TreeStructToNewick(int i);
		char *MakeNewick(int, bool);
		void CreateGnuPlotFile();
		void WritePopulationCheckpoint(ofstream &out);
		void ReadPopulationCheckpoint();
		void WriteStateFiles();
		void ReadStateFiles();
		void GetConstraints();
		void WriteTreeFile( const char* treefname, int fst = -1, int lst = -1 );
		void WritePhylipTree(ofstream &phytree);

		void Setup(GeneralGamlConfig *conf, HKYData *, int nprocs = 1, int rank = 0);
		int Restart(int type, int rank, int nprocs, int restart_count);
		void SeedPopulationWithStartingTree();
		FLOAT_TYPE CalcAverageFitness();
		void CalculateReproductionProbabilies(FLOAT_TYPE **scoreArray, FLOAT_TYPE selectionIntensity, int indivsInArray);
		void NextGeneration();
		void DetermineParentage();
		void FindTreeStructsForNextGeneration();
		void PerformMutation(int indNum);

		int IsError() const { return error; }
		void ErrorMsg( char* msgstr, int len );
		void CompactTopologiesList();
		void EliminateDuplicateTreeReferences();
		friend istream& operator >>( istream& inf, Population& p );
		friend ostream& operator <<( ostream& outf, Population& p );

		int ExtendPopulation(int, char*, FLOAT_TYPE*);
		int ShrinkPopulation(int, char**, FLOAT_TYPE**);
		int SwapIndividuals(int, const char*, FLOAT_TYPE*, char**, FLOAT_TYPE**);
		int ReplaceSpecifiedIndividuals(int,  int*, const char*, FLOAT_TYPE*);
		int ReplicateSpecifiedIndividuals(int, int*, const char*, FLOAT_TYPE *);
		void FillPopWithClonesOfBest();
		
		int GetNRandomIndivIndices(int**, int);
		int GetNBestIndivIndices(int**, int);
		int GetSpecifiedTreeStrings(char**, int, int*);
		int GetSpecifiedRates(FLOAT_TYPE**, int, int*);
		int GetSpecifiedPis(FLOAT_TYPE**, int , int*);
		int GetSpecifiedModels(FLOAT_TYPE** model_string, int n, int* indiv_list);
		
		void UpdateTopologyList(Individual *inds);
		void CheckAllTrees();
		void CheckIndividuals();
		void TopologyReport();
		void RemoveFromTopologyList(Individual *ind);
		void SetupTopologyList(int maxNumTopos);
		void CheckTreesVsClaManager();
		FLOAT_TYPE IndivFitness(int i);
		
		void NNIPerturbation(int sourceInd, int indivIndex);
		void NNISpectrum(int sourceInd);

		void NNIoptimization();
//		void SPRoptimization(int indivIndex);
		bool NNIoptimization(unsigned IndivIndex, int steps);
//		bool SPRoptimization(int indivIndex, int range, int cutnum );
		void SPRPerturbation(int sourceInd, int indivIndex);
		void keepTrack();
		void DetermineSubsets(int);
		void Partition(TreeNode *pointer);
		void PartitionDown(TreeNode *pointer, TreeNode *calledFrom);

		void CheckPerturbSerial();
		void CheckPerturbParallel();
		void StoreBestForPert();
		void StoreAllTimeBest();
		void RestoreBestForPert();
		void RestoreAllTimeBest();
		void CheckSubtrees();
		void AssignSubtree(int st, int indNum);
		bool SubtreeRecombination(int);
		void StartSubtreeMode();
		void StopSubtreeMode();
		
	private:

		int prResizeIndividualArray(int, char* = NULL, FLOAT_TYPE* = NULL);
		int prResizeNewIndividualArray(int);
		int prResizeTopologyListArray(int);
		int prResizeCumFitArray(int);
				
	public:
		void InitializeOutputStreams();
		void FinalizeOutputStreams();

		void AppendTreeToTreeLog(int mutType, int indNum=-1);
		void FinishBootstrapRep(int rep);
		void UpdateTreeModels();
		
		void OutputFate();
		void OutputLog();
		void OutputModelReport();

		void OutputModelAddresses();
		void OutputClaReport(Individual *arr);
		void OutputFilesForScoreDebugging(Individual *ind=NULL, int num=0);

		void Bootstrap();
		void AssignNewTopology(Individual *ind, int indNum);
		void FindLostClas();
		void FinalOptimization();
		void ResetMemLevel(int numNodesPerIndiv, int numClas);
		void SetNewBestIndiv(int indivIndex);
		void LogNewBestFromRemote(FLOAT_TYPE, int);
		void CheckRemoteReplaceThresh();
		void TurnOffRatchet();
	};
#endif
