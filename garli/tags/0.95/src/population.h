// GARLI version 0.95b6 source code
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
	double blen;
	int numAssigned;
	int priority;
	double score;
	
	Subtree(int nn, int t, double b, double s){
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
	double subtreeDefScore;	//the score of the tree that the subtrees were initially calculated on
	int lastFullRecom; 		//the last time we tried a subtree based recombination
	vector<Subtree *> subtrees;
	bool fewNonSubtreeNodes;

	int *remoteSubtreeAssign;  // which node the remote was most recently assigned
	int *localSubtreeAssign;  // which node the remotes _were_ working on for the most recent tree received

	bool needUpdate;			//this means that the subtrees determined are no longer guaranteed to be valid for any trees in the population
	bool beforeFirstSubtree;

	double updateThresh;
	double startUpdateThresh;
	double minUpdateThresh;
	double updateReductionFactor;
	
	int subtreeInterval;
	double subtreeStartThresh;
	int	minSubtreeSize;
	int targetSubtreeSize;
	double orphanFactor;
	
	int maxRecomIndivs;

//	double recalcThresh;		//the score threshold for forcing a recalc of the subtrees
//	double subtreeThresh; 	//the score threshold for sending a remote a new tree if it 
								//is currently working on a subtree
//	double nonSubtreeThresh;	//the score threshold for sending a remote a new tree if subtree
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
		
		updateReductionFactor=pow(minUpdateThresh/startUpdateThresh, 1.0/ (mc->numPrecReductions));


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
	double ScorePartitioning(int nodeNum, ofstream &pscores);
	void PrepareForSubtreeMode(Individual *ind, int gen);
	int ChooseSubtree();
			
  };


class PerturbManager{
	public:
	int lastPertGeneration;
	bool pertAbandoned;
	int numPertsNoImprove;
	double prevPertScore;
	double scoreAfterRatchet;

	int pertType;
	double pertThresh;
	int minPertInterval;
	int maxPertsNoImprove;
	bool restartAfterAbandon;
	int gensBeforeRestart;

	double ratchetProportion;
	double ratchetOffThresh;
	int ratchetMaxGen;
	bool ratcheted;

	double nniAcceptThresh;
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

class Population
{
public:
	int total_size; //this will be equal to conf->nindiv, except in 
					//the case of the parallel master
					
	Individual* indiv;
	Individual* newindiv;
	double bestFitness;
	int rank;//denotes which processor this is.  0 if serial
	
	int subtreeNode;
	int subtreeDefNumber;
	vector<int> subtreeMemberNodes;
	
	PerturbManager *pertMan;
	ParallelManager *paraMan;
	
	ClaManager *claMan;
	Adaptation *adap;
	bool refineStart;
	bool outputTreelog;
	bool outputMostlyUselessFiles;
	bool outputPhylipTree;
	
	bool bootstrap;
	int bootstrapReps;
	bool inferInternalStateProbs;
	int bestIndiv;
	int bestAccurateIndiv;

	//termination related variables
	bool enforceTermConditions;
	int lastTopoImprove;
	int lastPrecisionReduction;
	int lastTopoImproveThresh;
	double improveOverStoredIntervalsThresh;
	double significantTopoChange;//the score difference from the current best required for 
								 //a new topology to really be considered "better"
	
private:
	//DJZ adding these streams directly to the class so that they can be opened once and left open
	ofstream fate;
	ofstream log;
	ofstream treeLog;
	ofstream probLog;
	ofstream bootLog;
	ofstream swapLog;
	char besttreefile[100];

	int ntopos;
	bool prematureTermination;//if the user killed the run

	char *treeString;
	long stringSize;
		
	vector<Tree *> unusedTrees;
	double prevBestFitness;
	double avgfit;
	int new_best_found;
		
	int numgensamebest;
	

	char* logfname;
	
	ofstream logf;

	private:
		void QuickSort( double **scoreArray, int top, int bottom );

	public:
		enum { nomem=1, nofile, baddimen };
		int error;

		double** cumfit;
		//allocated in setup, deleted in dest
		TopologyList **topologies;
			//allocated in Setup(), deleted in dest
		long gen;
		GeneralGamlConfig *conf;
		HKYData* data;
		Individual *allTimeBest; //this is only used for perturbation or ratcheting
		Individual *bestSinceRestart;
		Stopwatch stopwatch;
        	double starting_wtime;
            double final_wtime;

	public:
		Population() : error(0), conf(NULL),
			bestFitness(-(DBL_MAX)), bestIndiv(0),
			prevBestFitness(-(DBL_MAX)), logfname(0),
			indiv(NULL), newindiv(NULL),
			cumfit(NULL), new_best_found(0), avgfit(0.0),
			gen(0), starting_wtime(0.0), final_wtime(0.0),
#ifdef INCLUDE_PERTURBATION			 
			 pertMan(NULL),
#endif
			 paraMan(NULL), subtreeDefNumber(0), claMan(NULL), 
			 inferInternalStateProbs(0), bootstrapReps(0),
			 outputMostlyUselessFiles(0), outputPhylipTree(0),
			 significantTopoChange(0.01), allTimeBest(NULL),
			 bestSinceRestart(NULL), treeString(NULL), adap(NULL),
			 topologies(NULL), prematureTermination(false)
			{
			//allTimeBest.SetFitness(-1e100);
			//bestSinceRestart.SetFitness(-1e100);
			
			lastTopoImprove = 0;
			lastPrecisionReduction = 0;
			lastTopoImproveThresh = 1000;
			improveOverStoredIntervalsThresh = 0.1;
			}

		~Population();

		double BestFitness() { return bestFitness; }

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

		void Setup(GeneralGamlConfig *conf, HKYData *, int nprocs = 1, int rank = 0);
		int Restart(int type, int rank, int nprocs, int restart_count);
		void SeedPopulationWithStartingTree();
		double CalcAverageFitness();
		void CalculateReproductionProbabilies(double **scoreArray, double selectionIntensity, int indivsInArray);
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

		int ExtendPopulation(int, char*, double*);
		int ShrinkPopulation(int, char**, double**);
		int SwapIndividuals(int, const char*, double*, char**, double**);
		int ReplaceSpecifiedIndividuals(int,  int*, const char*, double*);
		int ReplicateSpecifiedIndividuals(int, int*, const char*, double *);
		void FillPopWithClonesOfBest();
		
		int GetNRandomIndivIndices(int**, int);
		int GetNBestIndivIndices(int**, int);
		int GetSpecifiedTreeStrings(char**, int, int*);
		int GetSpecifiedRates(double**, int, int*);
		int GetSpecifiedPis(double**, int , int*);
		int GetSpecifiedModels(double** model_string, int n, int* indiv_list);
		
		void UpdateTopologyList(Individual *inds);
		void CheckAllTrees();
		void CheckIndividuals();
		void TopologyReport();
		void RemoveFromTopologyList(Individual *ind);
		void SetupTopologyList(int maxNumTopos);
		void CheckTreesVsClaManager();
		double IndivFitness(int i);
		
		void NNIPerturbation(int sourceInd, int indivIndex);
		void NNISpectrum(int sourceInd);

		void NNIoptimization();
//		void SPRoptimization(int indivIndex);
		bool NNIoptimization(int IndivIndex, int steps);
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

		int prResizeIndividualArray(int, char* = NULL, double* = NULL);
		int prResizeNewIndividualArray(int);
		int prResizeTopologyListArray(int);
		int prResizeCumFitArray(int);
				
	public:
		void InitializeOutputStreams();
		void FinalizeOutputStreams();

		void AppendTreeToTreeLog(int mutType, int indNum=-1);
		void AppendTreeToBootstrapLog(int rep);
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
		void LogNewBestFromRemote(double, int);
		void CheckRemoteReplaceThresh();
		void TurnOffRatchet();
	};
#endif
