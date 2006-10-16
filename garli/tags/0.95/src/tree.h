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
		double lnL;	// holds likelihood score
		Model *mod;
		TreeNode *root;
		TreeNode **allNodes;
		ReconList sprRang;

		//a bunch of statics
		static double meanBrlenMuts;
		static double alpha; //alpha shape of blen mutation, not gamma rate het
		static double min_brlen;
		static double max_brlen;
		static double exp_starting_brlen;
		static ClaManager *claMan;
		static HKYData *data;
		static double treeRejectionThreshold;
		static vector<Constraint> constraints;
		static AttemptedSwapList attemptedSwaps;
		static double uniqueSwapBias;
		static double distanceSwapBias;
		static int rescaleEvery;
		static list<TreeNode *> nodeOptVector;

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
		inline void SetBranchLength(TreeNode *nd, double len);
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
		int TopologyMutator(double optPrecision, int range, int subtreeNode);
		void GatherValidReconnectionNodes(int maxRange, TreeNode *cut, const TreeNode *subtreeNode);
		void AssignWeightsToSwaps(TreeNode *cut);
		int SPRMutate(int cutnum, ReconNode *broke, double optPrecision, int subtreeNode);
		void ReorientSubtreeSPRMutate(int oldRoot, ReconNode *newRoot, double optPrecision);
		int BrlenMutate();
		int BrlenMutateSubset(const vector<int> &subtreeList);
		void ScaleWholeTree(double factor=-1.0);
		//deprecated mutation functions
		int VariableSPRMutate(int range, double optPrecision);
		void SPRMutate(int cutnum, int broknum, double optPrecision, const vector<int> &nonSubNodes);
		void NNIMutate(int node, int branch, double optPrecision, int subtreeNode);
		void VariableNNIMutate(int node, int branch, double optPrecision, int subtreeNode);
		void LocalMove();

		//recombination
		int BipartitionBasedRecombination( Tree *t, bool sameModel, double optPrecision);
		int SubtreeBasedRecombination( Tree *t, int recomNodeNum, bool sameModel, double optPrecision);
		void RecombineWith( Tree *t, bool sameModel , double optPrecision );

		//functions for dealing with constraints and bipartitions
		void LoadConstraints(ifstream &con);
		bool AllowedByConstraint(Constraint *constr, TreeNode *cut, ReconNode *broken, Bipartition &proposed) const;
		bool AllowedByPositiveConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *cut, TreeNode *broken);
		bool AllowedByNegativeConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *cut, TreeNode *broken);
		bool RecursiveAllowedByPositiveConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *nd);
		bool RecursiveAllowedByNegativeConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *nd);
		void CalcBipartitions();
		void OutputBipartitions();
		TreeNode *ContainsBipartition(const Bipartition *bip);
		TreeNode *ContainsBipartitionOrComplement(const Bipartition *bip);

		// functions for computing likelihood
		bool ConditionalLikelihood(int direction, TreeNode* nd);	
		int ConditionalLikelihoodRateHet(int direction, TreeNode* nd, bool fillFinalCLA=false);
		double GetScorePartialTerminalRateHet(const CondLikeArray *partialCLA, const double *prmat, const char *Ldata);
		double GetScorePartialInternalRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const double *prmat);
		int Score(int rootNodeNum =0);
	
		//functions to optimize blens and params
		pair<double, double> CalcDerivativesRateHet(TreeNode *nd1, TreeNode *nd2);
		double NewtonRaphsonOptimizeBranchLength(double precision1, TreeNode *nd, bool goodGuess);
#ifdef OPT_DEBUG
		double NewtonRaphsonSpoof(double precision1, TreeNode *nd, bool goodGuess);
#endif
		void GetDerivsPartialTerminal(const CondLikeArray *partialCLA, const double *prmat, const double *d1mat, const double *d2mat, const char *Ldata, double &d1Tot, double &d2Tot, const unsigned *ambigMap =NULL);
		void GetDerivsPartialInternal(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const double *prmat, const double *d1mat, const double *d2mat, double &d1, double &d2);
		double OptimizeBranchLength(double optPrecision, TreeNode *nd, bool goodGuess);
		double OptimizeAllBranches(double optPrecision);
		void OptimizeBranchesAroundNode(TreeNode *nd, double optPrecision, int subtreeNode);
		void OptimizeBranchesWithinRadius(TreeNode *nd, double optPrecision, int subtreeNode, TreeNode *prune);
		void OptimizeBranchesInArray(int *nodes, int numNodes, double optPrecision);
		double RecursivelyOptimizeBranches(TreeNode *nd, double optPrecision, int subtreeNode, int radius, bool dontGoNext, double scoreIncrease, bool ignoreDelta=false);
		double RecursivelyOptimizeBranchesDown(TreeNode *nd, TreeNode *calledFrom, double optPrecision, int subtreeNode, int radius, double scoreIncrease);
		double BrentOptimizeBranchLength(double accuracy_cutoff, TreeNode *here, bool firstPass);
		double BranchLike(TreeNode *optNode);
		double OptimizeAlpha();
		double OptimizeTreeScale(double);
		double OptimizePinv();
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
		double CountClasInUse();
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
		void FindNearestTerminalUp(TreeNode *start, TreeNode *&, double &dist);
		void FindNearestTerminalsDown(TreeNode *start, TreeNode *from, TreeNode *&term1, TreeNode *&term2, double &dist1, double &dist2);
		void OutputTreeStructure(TreeNode *);
		void GetInternalStateString(char *string, int nodeNum);
		void RecursivelyCalculateInternalStateProbs(TreeNode *nd, ofstream &out);	
		void InferAllInternalStateProbs(const char *ofprefix);

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
	
inline void Tree::SetBranchLength(TreeNode *nd, double len){
	assert(!(len < DEF_MIN_BRLEN) && !(len > DEF_MAX_BRLEN));
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
			if(allNodes[i]->left->left!=NULL && allNodes[i]->right->left!=NULL) claMan->ReserveCla(allNodes[i]->claIndexDown, false);
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

