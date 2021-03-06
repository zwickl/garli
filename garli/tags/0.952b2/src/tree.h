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
#include "mlhky.h"
#include "reconnode.h"


#undef BRENT

class DNAData;
class HKYData;
class ClaManager;
class GeneralGamlConfig;
class Model;
extern rng rnd;

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

		//a bunch of statics
		static FLOAT_TYPE meanBrlenMuts;
		static FLOAT_TYPE alpha; //alpha shape of blen mutation, not gamma rate het
		static FLOAT_TYPE min_brlen;
		static FLOAT_TYPE max_brlen;
		static FLOAT_TYPE exp_starting_brlen;
		static ClaManager *claMan;
		static const HKYData *data;
		static FLOAT_TYPE treeRejectionThreshold;
		static vector<Constraint> constraints;
		static AttemptedSwapList attemptedSwaps;
		static FLOAT_TYPE uniqueSwapBias;
		static FLOAT_TYPE distanceSwapBias;
		static unsigned rescaleEvery;
		static list<TreeNode *> nodeOptVector;
		
		static float uniqueSwapPrecalc[500];
		static float distanceSwapPrecalc[1000];

		int calcs;

	enum{//the directions for sweeping of CLAs
		DOWN = 1,
		UPLEFT = 2,
		UPRIGHT = 3,
		ROOT = 4
		};		
		
	public: 
		//construction and allocation functions
		Tree();
		Tree(HKYData*,CondLikeArray **sharedcl);
		Tree(const char*, bool numericalTaxa, bool allowPolytomies=false);
		void AllocateTree();
		void AssignCLAsFromMaster();
	
		//destructor
		~Tree();

		//functions for manipulating and making trees
		void AddRandomNode(int nodenum, int & );
		void AddRandomNodeWithConstraints(int nodenum, int &placeInAllNodes, Bipartition &mask);
		void MakeTrifurcatingRoot(bool reducenodes, bool clasAssigned);
		void SortAllNodesArray();
		void EliminateNode(int nn);
		int FindUnusedNode(int start);
		inline void SetBranchLength(TreeNode *nd, FLOAT_TYPE len);
		bool IdenticalSubtreeTopology(const TreeNode *other);
		bool IdenticalTopology(const TreeNode *other);
		void RotateNodesAtRoot(TreeNode *newroot);
		void RerootHere(int newroot);
		void SwapNodeDataForReroot(TreeNode *nroot);
		void CheckBalance();
		void SwapAndFreeNodes(TreeNode *cop);

		//functions for copying trees
		void MimicTopologyButNotInternNodeNums(TreeNode *copySource,TreeNode *replicate,int &placeInAllNodes);
		void MimicTopo(TreeNode *, bool firstNode, bool sameModel);	
     	void MimicTopo(const Tree *source);
        void CopyBranchLens(const Tree *s);
		void CopyClaIndeces(const Tree *source, bool remove);

		// mutation functions
		int TopologyMutator(FLOAT_TYPE optPrecision, int range, int subtreeNode);
		void GatherValidReconnectionNodes(int maxRange, TreeNode *cut, const TreeNode *subtreeNode);
		void AssignWeightsToSwaps(TreeNode *cut);
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
		void LoadConstraints(ifstream &con, int nTaxa);
		bool AllowedByConstraint(Constraint *constr, TreeNode *cut, ReconNode *broken, Bipartition &proposed) ;
		bool AllowedByPositiveConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *cut, TreeNode *broken);
		bool AllowedByPositiveBackboneConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *cut, TreeNode *broken);
		bool AllowedByNegativeConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *cut, TreeNode *broken);
		bool RecursiveAllowedByPositiveConstraintWithMask(Constraint *constr, const Bipartition *mask, TreeNode *nd);
		bool RecursiveAllowedByNegativeConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *nd);
		void CalcBipartitions();
		void OutputBipartitions();
		TreeNode *ContainsBipartition(const Bipartition *bip);
		TreeNode *ContainsBipartitionOrComplement(const Bipartition *bip);
		TreeNode *ContainsMaskedBipartitionOrComplement(const Bipartition *bip, const Bipartition *mask);

		// functions for computing likelihood
		bool ConditionalLikelihood(int direction, TreeNode* nd);	
		int ConditionalLikelihoodRateHet(int direction, TreeNode* nd, bool fillFinalCLA=false);
		FLOAT_TYPE GetScorePartialTerminalRateHet(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const char *Ldata);
		FLOAT_TYPE GetScorePartialInternalRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat);
		int Score(int rootNodeNum =0);
	
		//functions to optimize blens and params
		pair<FLOAT_TYPE, FLOAT_TYPE> CalcDerivativesRateHet(TreeNode *nd1, TreeNode *nd2);
		FLOAT_TYPE NewtonRaphsonOptimizeBranchLength(FLOAT_TYPE precision1, TreeNode *nd, bool goodGuess);
#ifdef OPT_DEBUG
		FLOAT_TYPE NewtonRaphsonSpoof(FLOAT_TYPE precision1, TreeNode *nd, bool goodGuess);
#endif
		void GetDerivsPartialTerminal(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, const char *Ldata, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot, const unsigned *ambigMap =NULL);
		void GetDerivsPartialInternal(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, FLOAT_TYPE &d1, FLOAT_TYPE &d2);
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
		FLOAT_TYPE OptimizeTreeScale(FLOAT_TYPE);
		FLOAT_TYPE OptimizePinv();
		void SetNodesUnoptimized();

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
		void OutputFirstClaAcrossTree(ofstream &deb, TreeNode *nd);
		void ClaReport(ofstream &cla);
		FLOAT_TYPE CountClasInUse();
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
		int NodeToNodeDistance(int num1, int num2);
		int NodesToRoot(TreeNode *nd);
		void SampleBlenCurve(TreeNode *nd, ofstream &out);
		void SetDistanceBasedBranchLengthsAroundNode(TreeNode *nd);
		void FindNearestTerminalUp(TreeNode *start, TreeNode *&, FLOAT_TYPE &dist);
		void FindNearestTerminalsDown(TreeNode *start, TreeNode *from, TreeNode *&term1, TreeNode *&term2, FLOAT_TYPE &dist1, FLOAT_TYPE &dist2);
		void OutputTreeStructure(TreeNode *);
		void GetInternalStateString(char *string, int nodeNum);
		void RecursivelyCalculateInternalStateProbs(TreeNode *nd, ofstream &out);	
		void InferAllInternalStateProbs(const char *ofprefix);

		static void SetTreeStatics(ClaManager *, const HKYData *, const GeneralGamlConfig *);
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
	lnL=-1;
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

