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


#ifndef __TREE_H
#define __TREE_H

#include <vector>
#include <list>

using namespace std;

#include "defs.h"
#include "rng.h"
#include "treenode.h"
#include "clamanager.h"
#include "hashdefines.h"
#include "model.h"
#include "sequencedata.h"
#include "reconnode.h"


#undef BRENT

class SequenceData;
class NucleotideData;
class ClaManager;
class GeneralGamlConfig;
class Model;
class Individual;
extern rng rnd;

#define RESCALE_ARRAY_LENGTH 60

class Tree{
	protected:
		int numTipsTotal;
		int numTipsAdded;
		int numNodesAdded;
		int numBranchesAdded;
		int numNodesTotal;
		int *taxtags;//int[ntax+1] used in tagging terminals in recombine function
			//allocated in SharedTreeConstruction, deleted in dest

	public:
		FLOAT_TYPE lnL;	// holds likelihood score
		Model *mod;
		TreeNode *root;
		TreeNode **allNodes;
		ReconList sprRang;

#ifdef EQUIV_CALCS
		bool dirtyEQ;
#endif
		//a bunch of statics
		static FLOAT_TYPE meanBrlenMuts;
		static FLOAT_TYPE alpha; //alpha shape of blen mutation, not gamma rate het
		static FLOAT_TYPE min_brlen;
		static FLOAT_TYPE max_brlen;
		static FLOAT_TYPE exp_starting_brlen;
		static ClaManager *claMan;
		static const SequenceData *data;
		static FLOAT_TYPE treeRejectionThreshold;
		static vector<Constraint> constraints;
		static AttemptedSwapList attemptedSwaps;
		static FLOAT_TYPE uniqueSwapBias;
		static FLOAT_TYPE distanceSwapBias;
		static unsigned rescaleEvery;
		static FLOAT_TYPE rescaleBelow;
		static list<TreeNode *> nodeOptVector;
		
		static FLOAT_TYPE uniqueSwapPrecalc[500];
		static FLOAT_TYPE distanceSwapPrecalc[1000];
//DEBUG		
//		static FLOAT_TYPE rescalePrecalcThresh[30];
//		static FLOAT_TYPE rescalePrecalcMult[30];
//		static int rescalePrecalcIncr[30];

		static FLOAT_TYPE rescalePrecalcThresh[RESCALE_ARRAY_LENGTH];
		static FLOAT_TYPE rescalePrecalcMult[RESCALE_ARRAY_LENGTH];
		static int rescalePrecalcIncr[RESCALE_ARRAY_LENGTH];

		static Bipartition *outgroup;

		static int siteToScore;

		int calcs;

	enum{//the directions for sweeping of CLAs
		DOWN = 1,
		UPLEFT = 2,
		UPRIGHT = 3,
		ROOT = 4
		};

	enum{
		DIRTY = 0,
		CLEAN_STANDARDIZED = 1,
		CLEAN_UNSTANDARDIZED = 2,
		TEMP_ADJUSTED = 3
		}bipartCond;
		
	public: 
		//construction and allocation functions
		Tree();
		Tree(NucleotideData*,CondLikeArray **sharedcl);
		Tree(const char*, bool numericalTaxa, bool allowPolytomies=false, bool allowMissingTaxa=false);
		void AllocateTree();
		void AssignCLAsFromMaster();
	
		//destructor
		~Tree();

		//functions for manipulating and making trees
		void AddRandomNode(int nodenum, int & );
		void AddRandomNodeWithConstraints(int nodenum, int &placeInAllNodes, Bipartition *mask);
		void MakeTrifurcatingRoot(bool reducenodes, bool clasAssigned);
		bool ArbitrarilyBifurcate();
		void SortAllNodesArray();
		void EliminateNode(int nn);
		int FindUnusedNode(int start);
		inline void SetBranchLength(TreeNode *nd, FLOAT_TYPE len);
		bool IdenticalSubtreeTopology(const TreeNode *other);
		bool IdenticalTopology(const TreeNode *other);
		bool IdenticalTopologyAllowingRerooting(const TreeNode *other);
		void RotateNodesAtRoot(TreeNode *newroot);
		void RerootHere(int newroot);
		void SwapNodeDataForReroot(TreeNode *nroot);
		void CheckBalance();
		void SwapAndFreeNodes(TreeNode *cop);
		void OutputBinaryFormattedTree(OUTPUT_CLASS &) const;
		void ReadBinaryFormattedTree(FILE *);

		//functions for copying trees
		void MimicTopologyButNotInternNodeNums(TreeNode *copySource,TreeNode *replicate,int &placeInAllNodes);
		void MimicTopo(TreeNode *, bool firstNode, bool sameModel);	
     	void MimicTopo(const Tree *source);
        void CopyBranchLens(const Tree *s);
		void CopyClaIndeces(const Tree *source, bool remove);

		// mutation functions
		int TopologyMutator(FLOAT_TYPE optPrecision, int range, int subtreeNode);
		void DeterministicSwapperByDist(Individual *source, double optPrecision, int range, bool furthestFirst);
		void DeterministicSwapperByCut(Individual *source, double optPrecision, int range, bool furthestFirst);
		void DeterministicSwapperRandom(Individual *source, double optPrecision, int range);
		void GatherValidReconnectionNodes(ReconList &list, int maxDist, TreeNode *cut, const TreeNode *subtreeNode, Bipartition *partialMask=NULL);
		void GatherValidReconnectionNodes(int maxRange, TreeNode *cut, const TreeNode *subtreeNode, Bipartition *partialMask=NULL);
		void FillAllSwapsList(ReconList *cuts, int reconLim);
		unsigned FillWeightsForAllSwaps(ReconList *cuts, double *);
		bool AssignWeightsToSwaps(TreeNode *cut);
		int SPRMutate(int cutnum, ReconNode *broke, FLOAT_TYPE optPrecision, int subtreeNode);
		int SPRMutateDummy(int cutnum, ReconNode *broke, FLOAT_TYPE optPrecision, int subtreeNode);
		void ReorientSubtreeSPRMutate(int oldRoot, ReconNode *newRoot, FLOAT_TYPE optPrecision);
		void ReorientSubtreeSPRMutateDummy(int oldRoot, ReconNode *newRoot, FLOAT_TYPE optPrecision);
		int BrlenMutate();
		int BrlenMutateSubset(const vector<int> &subtreeList);
		void ScaleWholeTree(FLOAT_TYPE factor=-1.0);
		//deprecated mutation functions
		int VariableSPRMutate(int range, FLOAT_TYPE optPrecision);
		void SPRMutate(int cutnum, int broknum, FLOAT_TYPE optPrecision, const vector<int> &nonSubNodes);
		void NNIMutate(int node, int branch, FLOAT_TYPE optPrecision, int subtreeNode);
		void VariableNNIMutate(int node, int branch, FLOAT_TYPE optPrecision, int subtreeNode);
		void LocalMove();

		//recombination
		int BipartitionBasedRecombination( Tree *t, bool sameModel, FLOAT_TYPE optPrecision);
		int SubtreeBasedRecombination( Tree *t, int recomNodeNum, bool sameModel, FLOAT_TYPE optPrecision);
		void RecombineWith( Tree *t, bool sameModel , FLOAT_TYPE optPrecision );

		//functions for dealing with constraints and bipartitions
		static void LoadConstraints(ifstream &con, int nTaxa);
		bool SwapAllowedByConstraint(const Constraint &constr, TreeNode *cut, ReconNode *broken, const Bipartition &proposed, const Bipartition *partialMask);
		
		//functions for determining if adding a particular taxon to a particular place in a growing tree is allowed by any constraints
		//this is used if the taxon to be added is NOT already in the tree (so not for testing the allowability of swaps)
		//these depend on the bipartitions across the tree NOT being standardized (not all having the taxon 1 bit on)
		bool TaxonAdditionAllowedByPositiveConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *toAdd, TreeNode *broken);
		bool TaxonAdditionAllowedByNegativeConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *toAdd, TreeNode *broken);
		bool TaxonAdditionAllowedByPositiveBackboneConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *toAdd, TreeNode *broken);
		bool TaxonAdditionAllowedByNegativeBackboneConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *toAdd, TreeNode *broken);
		
		bool RecursiveAllowedByConstraintWithMask(const Constraint &constr, const Bipartition *partialMask, const TreeNode *nd);
		//bool RecursiveAllowedByNegativeConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *nd);
		void CalcBipartitions(bool standardize);
		void OutputBipartitions();
		TreeNode *ContainsBipartition(const Bipartition &bip);
		TreeNode *ContainsBipartitionOrComplement(const Bipartition &bip);
		TreeNode *ContainsMaskedBipartitionOrComplement(const Bipartition &bip, const Bipartition &mask);
		void AdjustBipartsForSwap(int cut, int broken);

		// functions for computing likelihood
		bool ConditionalLikelihood(int direction, TreeNode* nd);	
		int ConditionalLikelihoodRateHet(int direction, TreeNode* nd, bool fillFinalCLA=false);
		FLOAT_TYPE GetScorePartialTerminalRateHet(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const char *Ldata);
		FLOAT_TYPE GetScorePartialTerminalNState(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const char *Ldata);
		FLOAT_TYPE GetScorePartialInternalRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat);
		FLOAT_TYPE GetScorePartialInternalNState(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat);
		int Score(int rootNodeNum =0);
	
		//functions to optimize blens and params
		pair<FLOAT_TYPE, FLOAT_TYPE> CalcDerivativesRateHet(TreeNode *nd1, TreeNode *nd2);
		FLOAT_TYPE NewtonRaphsonOptimizeBranchLength(FLOAT_TYPE precision1, TreeNode *nd, bool goodGuess);
#ifdef OPT_DEBUG
		FLOAT_TYPE NewtonRaphsonSpoof(FLOAT_TYPE precision1, TreeNode *nd, bool goodGuess);
#endif
		void GetDerivsPartialTerminal(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, const char *Ldata, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot, const unsigned *ambigMap =NULL);
		void GetDerivsPartialTerminalNState(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, const char *Ldata, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot);
		void GetDerivsPartialTerminalNStateRateHet(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, const char *Ldata, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot);
		void GetDerivsPartialInternal(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, FLOAT_TYPE &d1, FLOAT_TYPE &d2);
		void GetDerivsPartialInternalNState(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, FLOAT_TYPE &d1, FLOAT_TYPE &d2);
		void GetDerivsPartialInternalNStateRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot);
		void GetDerivsPartialInternalEQUIV(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, FLOAT_TYPE &d1, FLOAT_TYPE &d2, char *equiv);
		void CalcFullCLAInternalInternal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr);
		void CalcFullCLATerminalTerminal(CondLikeArray *destCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, const char *Ldata, const char *Rdata);
		void CalcFullCLAInternalTerminal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const FLOAT_TYPE *pr1, const FLOAT_TYPE *pr2, char *data2, const unsigned *ambigMap);

		void CalcFullCLAInternalInternalEQUIV(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr,const char *leftEQ, const char *rightEQ);

		void CalcFullCLAPartialInternalRateHet(CondLikeArray *destCLA, const CondLikeArray *LCLA, const FLOAT_TYPE *pr1, CondLikeArray *partialCLA);
		void CalcFullCLAPartialTerminalRateHet(CondLikeArray *destCLA, const CondLikeArray *partialCLA, const FLOAT_TYPE *Lpr, char *Ldata);

		void CalcFullCLAInternalInternalNState(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr);
		void CalcFullCLAInternalTerminalNState(CondLikeArray *destCLA, const CondLikeArray *LCLA, const FLOAT_TYPE *pr1, const FLOAT_TYPE *pr2, char *data2);
		void CalcFullCLATerminalTerminalNState(CondLikeArray *destCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, const char *Ldata, const char *Rdata);
		
		FLOAT_TYPE OptimizeBranchLength(FLOAT_TYPE optPrecision, TreeNode *nd, bool goodGuess);
		FLOAT_TYPE OptimizeAllBranches(FLOAT_TYPE optPrecision);
		void OptimizeBranchesAroundNode(TreeNode *nd, FLOAT_TYPE optPrecision, int subtreeNode);
		void OptimizeBranchesWithinRadius(TreeNode *nd, FLOAT_TYPE optPrecision, int subtreeNode, TreeNode *prune);
		void OptimizeBranchesInArray(int *nodes, int numNodes, FLOAT_TYPE optPrecision);
		FLOAT_TYPE RecursivelyOptimizeBranches(TreeNode *nd, FLOAT_TYPE optPrecision, int subtreeNode, int radius, bool dontGoNext, FLOAT_TYPE scoreIncrease, bool ignoreDelta=false);
		FLOAT_TYPE RecursivelyOptimizeBranchesDown(TreeNode *nd, TreeNode *calledFrom, FLOAT_TYPE optPrecision, int subtreeNode, int radius, FLOAT_TYPE scoreIncrease);
		FLOAT_TYPE BrentOptimizeBranchLength(FLOAT_TYPE accuracy_cutoff, TreeNode *here, bool firstPass);
		FLOAT_TYPE BranchLike(TreeNode *optNode);
		FLOAT_TYPE OptimizeAlpha(FLOAT_TYPE);
		FLOAT_TYPE OptimizeOmegaParameters(FLOAT_TYPE prec);
		FLOAT_TYPE OptimizeFlexRates(FLOAT_TYPE prec);
		FLOAT_TYPE OptimizeBoundedParameter(FLOAT_TYPE optPrecision, FLOAT_TYPE prevVal, int which, FLOAT_TYPE lowBound, FLOAT_TYPE highBound, void (Model::*SetParam)(int, FLOAT_TYPE));
		FLOAT_TYPE OptimizeTreeScale(FLOAT_TYPE);
		FLOAT_TYPE OptimizePinv();
		void SetNodesUnoptimized();
		void RescaleRateHet(CondLikeArray *destCLA);
		void RescaleRateHetNState(CondLikeArray *destCLA);

		pair<FLOAT_TYPE, FLOAT_TYPE> OptimizeSingleSiteTreeScale(FLOAT_TYPE optPrecision);

		//functions for dealing with conditional likelihood arrays
		void MarkUpwardClasToReclaim(int subtreeNode);
		void MarkDownwardClasToReclaim(int subtreeNode);
		void MarkClasNearTipsToReclaim(int subtreeNode);
		void ProtectClas();
		void UnprotectClas();
		inline CondLikeArray *GetClaDown(TreeNode *nd, bool calc=true);
		inline CondLikeArray *GetClaUpLeft(TreeNode *nd, bool calc=true);
		inline CondLikeArray *GetClaUpRight(TreeNode *nd, bool calc=true);
		void OutputValidClaIndeces();
		void OutputNthClaAcrossTree(ofstream &deb, TreeNode *nd, int site);
		void ClaReport(ofstream &cla);
		FLOAT_TYPE CountClasInUse();
		void OutputSiteLikelihoods(vector<double> &likes, const int *under1, const int *under2, ofstream &ordered, ofstream &packed);
		void OutputSiteDerivatives(vector<double> &likes, vector<double> &d1s, vector<double> &d2s, const int *under1, const int *under2, ofstream &ordered, ofstream &packed);
		void CountNumReservedClas(int &, int &, int&);
		void CheckClaAssignments(TreeNode *nd);
		void RemoveTempClaReservations();
		void SetupClasForSubtreeMode(int subtreeNode);
		void DirtyNodesOutsideOfSubtree(TreeNode *nd, int subtreeNodeAnc);
		void CopyClaIndecesInSubtree(const TreeNode *from, bool remove);
		void DirtyNodesInSubtree(TreeNode *nd);
		void ReclaimUniqueClas();
		void RemoveTreeFromAllClas();
		void TraceDirtynessToRoot(TreeNode *nd);
		void TraceDirtynessToNode(TreeNode *nd, int tonode);
		void SweepDirtynessOverTree(TreeNode *nd, TreeNode *from=NULL);
		void MakeNodeDirty(TreeNode *nd);
		void MakeAllNodesDirty();

		//accessor funcs
		bool IsGood() const {return root->IsGood();}
		int getNumTipsTotal() const {return numTipsTotal;}
		int getNumNodesTotal() const {return numNodesTotal;}
		int GetRandomInternalNode() const {return numTipsTotal+rnd.random_int(numTipsTotal-3)+1;}
		int GetRandomTerminalNode() const {return rnd.random_int(numTipsTotal)+1;}
		int GetRandomNonRootNode() const {return rnd.random_int(numNodesTotal-1)+1;}

		//odds and ends
		void PerturbAllBranches();
		void RandomizeBranchLengths(FLOAT_TYPE lowLimit, FLOAT_TYPE highLimit);
		void RandomizeBranchLengthsExponential(FLOAT_TYPE lambda);
		int NodeToNodeDistance(int num1, int num2);
		int NodesToRoot(TreeNode *nd);
		void SampleBlenCurve(TreeNode *nd, ofstream &out);
		void CalcEmpiricalDerivatives(TreeNode *nd, FLOAT_TYPE &D1, FLOAT_TYPE &D2);
		void SetDistanceBasedBranchLengthsAroundNode(TreeNode *nd);
		void FindNearestTerminalUp(TreeNode *start, TreeNode *&, FLOAT_TYPE &dist);
		void FindNearestTerminalsDown(TreeNode *start, TreeNode *from, TreeNode *&term1, TreeNode *&term2, FLOAT_TYPE &dist1, FLOAT_TYPE &dist2);
		void OutputTreeStructure(TreeNode *);
		void GetInternalStateString(char *string, int nodeNum);
		void RecursivelyCalculateInternalStateProbs(TreeNode *nd, ofstream &out);	
		void InferAllInternalStateProbs(const char *ofprefix);

		static void SetTreeStatics(ClaManager *, const SequenceData *, const GeneralGamlConfig *);
		};


inline void Tree::CopyBranchLens(const Tree *s){
	for(int i=1;i<numNodesTotal;i++)
		allNodes[i]->dlen=s->allNodes[i]->dlen;
	}

inline void Tree::MakeAllNodesDirty(){
	root->claIndexDown=claMan->SetDirty(root->claIndexDown);
	root->claIndexUL=claMan->SetDirty(root->claIndexUL);
	root->claIndexUR=claMan->SetDirty(root->claIndexUR);
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		allNodes[i]->claIndexDown=claMan->SetDirty(allNodes[i]->claIndexDown);
		allNodes[i]->claIndexUL=claMan->SetDirty(allNodes[i]->claIndexUL);
		allNodes[i]->claIndexUR=claMan->SetDirty(allNodes[i]->claIndexUR);
		}
	lnL=-ONE_POINT_ZERO;
	}
	
inline int Tree::FindUnusedNode(int start){
	for(int i=start;i<numNodesTotal;i++)
		if(!(allNodes[i]->attached))
			{allNodes[i]->left=allNodes[i]->right=NULL;
			return i;
			}
	assert(0);
	return -1;
	}	

inline void Tree::AssignCLAsFromMaster(){
	//remember that the root's down cla is actually the one that goes up 
	//the middle des
	assert(allNodes[0]->claIndexDown==-1);
	allNodes[0]->claIndexDown=claMan->AssignClaHolder();
	allNodes[0]->claIndexUL=claMan->AssignClaHolder();
	allNodes[0]->claIndexUR=claMan->AssignClaHolder();
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		assert(allNodes[i]->claIndexDown==-1);
		allNodes[i]->claIndexDown=claMan->AssignClaHolder();
		allNodes[i]->claIndexUL=claMan->AssignClaHolder();
		allNodes[i]->claIndexUR=claMan->AssignClaHolder();
		}
	}
	
inline void Tree::CopyClaIndeces(const Tree *from, bool remove){
	//the bool argument "remove" designates whether the tree currently has cla arrays
	//assigned to it or not (if not, it must have come from the unused tree vector)

	//do the clas down
	if(remove) claMan->DecrementCla(allNodes[0]->claIndexDown);
	allNodes[0]->claIndexDown=from->allNodes[0]->claIndexDown;
	if(allNodes[0]->claIndexDown != -1) claMan->IncrementCla(allNodes[0]->claIndexDown);
	
#ifdef EQUIV_CALCS
	if(from->dirtyEQ == false){
		memcpy(allNodes[0]->tipData, from->allNodes[0]->tipData, data->NChar()*sizeof(char));
		for(int i=numTipsTotal+1;i<numNodesTotal;i++)
			memcpy(allNodes[i]->tipData, from->allNodes[i]->tipData, data->NChar()*sizeof(char));
		dirtyEQ = false;
		}
	else dirtyEQ = true;
#endif

	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		if(remove) claMan->DecrementCla(allNodes[i]->claIndexDown);
		allNodes[i]->claIndexDown=from->allNodes[i]->claIndexDown;
		if(allNodes[i]->claIndexDown != -1) claMan->IncrementCla(allNodes[i]->claIndexDown);
		}
		
	//do the clas up left
	if(remove) claMan->DecrementCla(allNodes[0]->claIndexUL);
	allNodes[0]->claIndexUL=from->allNodes[0]->claIndexUL;
	if(allNodes[0]->claIndexUL != -1) claMan->IncrementCla(allNodes[0]->claIndexUL);
	
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		if(remove) claMan->DecrementCla(allNodes[i]->claIndexUL);
		allNodes[i]->claIndexUL=from->allNodes[i]->claIndexUL;
		if(allNodes[i]->claIndexUL != -1) claMan->IncrementCla(allNodes[i]->claIndexUL);
		}
	
	//do the clas up right
	if(remove) claMan->DecrementCla(allNodes[0]->claIndexUR);
	allNodes[0]->claIndexUR=from->allNodes[0]->claIndexUR;
	if(allNodes[0]->claIndexUR != -1) claMan->IncrementCla(allNodes[0]->claIndexUR);
		
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		if(remove) claMan->DecrementCla(allNodes[i]->claIndexUR);
		allNodes[i]->claIndexUR=from->allNodes[i]->claIndexUR;
		if(allNodes[i]->claIndexUR != -1) claMan->IncrementCla(allNodes[i]->claIndexUR);
		}
	}

inline void Tree::RemoveTreeFromAllClas(){
	if(root->claIndexDown != -1){
		claMan->DecrementCla(root->claIndexDown);
		root->claIndexDown=-1;
		}
	if(root->claIndexUL != -1){
		claMan->DecrementCla(root->claIndexUL);
		root->claIndexUL=-1;
		}
	if(root->claIndexUR != -1){	
		claMan->DecrementCla(root->claIndexUR);
		root->claIndexUR=-1;
		}
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		if(allNodes[i]->claIndexDown != -1){
			claMan->DecrementCla(allNodes[i]->claIndexDown);
			allNodes[i]->claIndexDown=-1;
			}
		if(allNodes[i]->claIndexUL != -1){
			claMan->DecrementCla(allNodes[i]->claIndexUL);
			allNodes[i]->claIndexUL=-1;
			}
		if(allNodes[i]->claIndexUR != -1){
			claMan->DecrementCla(allNodes[i]->claIndexUR);
			allNodes[i]->claIndexUR=-1;
			}
		}
	}
	
inline void Tree::SetBranchLength(TreeNode *nd, FLOAT_TYPE len){
	assert(!(len < min_brlen) && !(len > max_brlen));
	nd->dlen=len;
	SweepDirtynessOverTree(nd);
	}

inline CondLikeArray *Tree::GetClaDown(TreeNode *nd, bool calc/*=true*/){
	if(claMan->IsDirty(nd->claIndexDown)){
		if(calc==true){
			ConditionalLikelihoodRateHet(DOWN, nd);
			}
		else claMan->FillHolder(nd->claIndexDown, 1);
		}
	if(memLevel > 1) claMan->ReserveCla(nd->claIndexDown);
	return claMan->GetCla(nd->claIndexDown);
	}
	
inline CondLikeArray *Tree::GetClaUpLeft(TreeNode *nd, bool calc/*=true*/){
	if(claMan->IsDirty(nd->claIndexUL)){
		if(calc==true){
			ConditionalLikelihoodRateHet(UPLEFT, nd);
			}
		else claMan->FillHolder(nd->claIndexUL, 2);
		}
	if(memLevel > 0) claMan->ReserveCla(nd->claIndexUL);
	return claMan->GetCla(nd->claIndexUL);
	}
	
inline CondLikeArray *Tree::GetClaUpRight(TreeNode *nd, bool calc/*=true*/){
	if(claMan->IsDirty(nd->claIndexUR)){
		if(calc==true){
			ConditionalLikelihoodRateHet(UPRIGHT, nd);
			}
		else claMan->FillHolder(nd->claIndexUR, 2);
		}
	if(memLevel > 0) claMan->ReserveCla(nd->claIndexUR);
	return claMan->GetCla(nd->claIndexUR);
	}

inline void Tree::ProtectClas(){
	if(memLevel != 3){
		for(int i=numTipsTotal+1;i<numNodesTotal;i++){
			claMan->ReserveCla(allNodes[i]->claIndexDown, false);
			}
		}
	else{
		for(int i=numTipsTotal+1;i<numNodesTotal;i++){
			if(allNodes[i]->left->IsInternal() && allNodes[i]->right->IsInternal()) claMan->ReserveCla(allNodes[i]->claIndexDown, false);
			}
		}
	}

inline void Tree::UnprotectClas(){
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		if(allNodes[i]->claIndexDown > -1)
			claMan->UnreserveCla(allNodes[i]->claIndexDown);
		}
	}

inline int Tree::NodeToNodeDistance(int num1, int num2){
	TreeNode *nd1=allNodes[num1];
	TreeNode *nd2=allNodes[num2];
	int dist=0;
	
	int height1=NodesToRoot(nd1);
	int height2=NodesToRoot(nd2);
	
	while(height1 > height2){
		nd1=nd1->anc;
		dist++;
		height1--;
		}
	while(height2 > height1){
		nd2=nd2->anc;
		dist++;
		height2--;
		}
	
	while(nd1 != nd2){
		nd1=nd1->anc;
		nd2=nd2->anc;
		dist += 2;
		}	
	
	return dist;
	}

inline int Tree::NodesToRoot(TreeNode *nd){
	int i=0;
	while(nd->anc){
		nd=nd->anc;
		i++;
		}
	return i;
	}

#endif

