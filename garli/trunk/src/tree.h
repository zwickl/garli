// GARLI version 0.93 source code
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


#ifndef __TREE_H
#define __TREE_H

#include <vector>
#include <list>

using namespace std;

#include "defs.h"
#include "rng.h"
#include "treenode.h"
#include "clamanager.h"
#include "subset.h"
#include "hashdefines.h"
#include "model.h"
#include "mlhky.h"
#include "reconnode.h"

#undef OPT_DEBUG
#undef BRENT

class DNAData;
class HKYData;
class Model;
extern rng rnd;

class Tree{

#ifdef GANESH
/* TODO explain these private variables */
/* TODO variables like A are declared as private, yet are passed around in
 * function calls unnecessarily. Identify and eliminate such redundancies.
 * */
    private:
        int *iedges;    
        int running_index;
        int pedges_index; 
        int **treecatalan;
        int n_p_subtree_tips;
        TreeNode *p_subtree_root;
        int *postorder;
        int *postorder_reverse;
        int postorder_index;

        int *masked_postorder;
        int *masked_postorder_reverse;
        int masked_postorder_index;

        TreeNode **p_subtree_leaves;
        int *pedges;
        int *A;
#endif

	protected:
		int numTipsTotal;
		int numTipsAdded;
		int numNodesAdded;
		int numBranchesAdded;
		int numNodesTotal;
		int *taxtags;//int[ntax+1] used in tagging terminals in recombine function
			//allocated in SharedTreeConstruction, deleted in dest

	public:
		TreeNode *root;
		TreeNode **allNodes;
		ReconList sprRang;
//		subset sprRange;
		static double meanBrlenMuts;
		static double alpha; //alpha shape of blen mutation, not gamma rate het
		static double min_brlen;
		static double max_brlen;
		static double exp_starting_brlen;
		static ClaManager *claMan;
		static HKYData *data;
		static double treeRejectionThreshold;
		static vector<Constraint> constraints;

		int calcs;
		
		void TopologyMutator(double optPrecision, int range, int subtreeNode);
		void GatherValidReconnectionNodes(int maxRange, TreeNode *cut, const TreeNode *subtreeNode);
		
	enum{
		DOWN = 1,
		UPLEFT = 2,
		UPRIGHT = 3,
		ROOT = 4
		};		
		
#ifdef GANESH
/* the value of p in pECR */
        static int p_value;
        static bool random_p;
        static int *C;
        static int **realcat_intervals;
        static int **inv_realcat_intervals;
#endif
			
		double lnL;	// holds likelihood score
		Model *mod;
		
		static int rescaleEvery;
		static list<TreeNode *> nodeOptVector;
		
	public: 
		Tree(HKYData*,CondLikeArray **sharedcl);
		Tree(const char*, bool numericalTaxa, bool allowPolytomies=false);
		~Tree();
		int FindUnusedNode(int start);
		void MimicTopologyButNotInternNodeNums(TreeNode *copySource,TreeNode *replicate,int &placeInAllNodes);
		void MimicTopo(TreeNode *, bool firstNode, bool sameModel);	
     	void MimicTopo(const Tree *source);
       
		// mutation and recombination functions
		Tree();
		void AllocateTree();
		void SharedTreeConstruction(CondLikeArray **sharedcl);
		void RecombineWith( Tree *t, bool sameModel , double optPrecision );
		void AddRandomNode(int nodenum, int & );
		void AddRandomNodeWithConstraints(int nodenum, int &placeInAllNodes, Bipartition &mask);
#ifdef GANESH
		void PECRMutate(rng& rnd, double optPrecision);
        int IdentifyPSubtree(TreeNode **p_subtree_leaves, int *pedges, int p, rng& rnd);
        void GenerateMatching(int *A, int matching_number,
                              int num_matchings, int num_p_subtree_vertices);
        void RootWithMiddle(int node_number);
        int  AlternateChoosePSubtreeGivenRoot(TreeNode **p_subtree_leaves,
                                              int *pedges,
                                              int p, TreeNode *subtree_root,  
                                               rng& rnd);
        int  ChoosePSubtreeGivenRoot(TreeNode **p_subtree_leaves,
                                     int *pedges,
                                     int p, TreeNode *subtree_root,  
                                     rng& rnd);
        static void ComputeRealCatalan();
        void ComputeTreeCatalan(int p);
        void ComputeNumInternalEdges(TreeNode *node, int *iedges);
        bool IsALeaf(TreeNode *node);
        bool IsALeaf(int nodeNum);
        void SortMatching(int *A, int max_index);
        void MatchingToTree(int *A, int max_index, TreeNode **p_subtree_leaves, 
                            int *pedges);
        void QSort(int *A, int left, int right);
        int Partition(int *A, int left, int right, int pivot_index);
        int Key(int *A, int i);
        int MaxLabelPair(int *A, int i);
        int Swap(int *A, int index1, int index2);
        int CompareLabelPairs(int *A, int index1, int index2);
        void PostOrderTraverse(TreeNode *node);
        void PECRCleanUp(int p, int code);
#ifdef PECR_SET_PARSIMONY_BRLEN
        void MaskedPostOrderTraverse(TreeNode *node);
        void ComputeParsimonyBrLen(int *pedges, TreeNode **p_subtree_leaves, int W[][4]);
        bool IsLeafMaskOn(TreeNode *node);
        bool IsLeafMaskOn(int nodeNum);
        void TranslateTipString(const char *tip, char *string, int nsites);
#endif
#endif
		inline void SetBranchLength(TreeNode *nd, double len);
		int VariableSPRMutate(int range, double optPrecision);
		int SPRMutate(int range, double optPrecision);
		int ConstrainedSPRMutate(int range, double optPrecision);
		int SPRMutate(int cutnum, int broknum, double optPrecision, int subtreeNode, int range);
		void SPRMutate(int cutnum, int broknum, double optPrecision, const vector<int> &nonSubNodes);
		void NNIMutate(int node, int branch, double optPrecision, int subtreeNode);
		void VariableNNIMutate(int node, int branch, double optPrecision, int subtreeNode);
		void TaxonSwap(int tax1, int tax2, double optPrecision);
		void TaxonSwap(int range, double optPrecision);
		void MakeTrifurcatingRoot(bool reducenodes, bool clasAssigned);
		int BrlenMutate();
		int BrlenMutateSubset(const vector<int> &subtreeList);
		void LocalMove();
		void LoadConstraints(ifstream &con);
		bool AllowedByConstraint(Constraint *constr, TreeNode *cut, TreeNode *broken) const;
		bool AllowedByPositiveConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *cut, TreeNode *broken);
		bool AllowedByNegativeConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *cut, TreeNode *broken);
		bool RecursiveAllowedByPositiveConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *nd);
		bool RecursiveAllowedByNegativeConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *nd);
		// functions for computing likelihood
		bool ConditionalLikelihood(int direction, TreeNode* nd);	
		int ConditionalLikelihoodRateHet(int direction, TreeNode* nd, bool fillFinalCLA=false);

		int Score(int rootNodeNum =0);
		double SumSiteLikes(const double *cla, const int *underflow_mult);
		void FillSiteLikes(const CondLikeArray *fullCla, double *dest, bool addPinv);
		int ScoreRateHet( const HKYData* data);
		double ScoreAtSubtree(int subtreeNode);
		bool IsGood();
		// operators
		friend ostream& operator <<( ostream&, Tree& );
		void CopyBranchLens(const Tree *s);
		void SortAllNodesArray();
		void EliminateNode(int nn);
		void AssignCLAs(CondLikeArray **sharedcl);
		void AllocateCLAs(int ntax, int nchar);
		void TraceDirtynessToRoot(TreeNode *nd);
		void TraceDirtynessToNode(TreeNode *nd, int tonode);
		void SweepDirtynessOverTree(TreeNode *nd, TreeNode *from=NULL);
		void MakeNodeDirty(TreeNode *nd);
		void MakeAllNodesDirty();
		void CopyCLAs(Tree *source, int nchar);
		void AssignCLAsFromMaster();
		void CopyClaIndeces(const Tree *source, bool remove);
		void RotateNodesAtRoot(TreeNode *newroot);
		void RerootHere(int newroot);
		void SwapNodeDataForReroot(TreeNode *nroot);
		void CheckBalance();
		void RemoveTreeFromAllClas();
		void CalcPmat(double *prmat, double t, double kappa, double beta, double piAG, double piCT, bool flip = 0);
		void SwapAndFreeNodes(TreeNode *cop);
		
		//functions to partially optimize blens
		pair<double, double> CalcDerivativesRateHet(TreeNode *nd1, TreeNode *nd2);
		double NewtonRaphsonOptimizeBranchLength(double precision1, TreeNode *nd, bool goodGuess);
		void GetDerivsPartialTerminalRateHet(const CondLikeArray *partialCLA, const double *prmat, const double *d1mat, const double *d2mat, const char *Ldata, double &d1, double &d2);
		void GetDerivsPartialTerminalFlexRates(const CondLikeArray *partialCLA, const double *prmat, const double *d1mat, const double *d2mat, const char *Ldata, double &d1Tot, double &d2Tot);
		void GetDerivsPartialInternalRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const double *prmat, const double *d1mat, const double *d2mat, double &d1, double &d2);
		void GetDerivsPartialInternalFlexRates(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const double *prmat, const double *d1mat, const double *d2mat, double &d1, double &d2);
		double OptimizeBranchLength(double optPrecision, TreeNode *nd, bool goodGuess);
		double OptimizeAlpha();
		double GetScorePartialTerminalRateHet(const CondLikeArray *partialCLA, const double *prmat, const char *Ldata);
		double GetScorePartialInternalRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const double *prmat);
		
		void PerturbAllBranches();
		void ScaleWholeTree(double factor=-1.0);
		double OptimizeTreeScale();
		double OptimizePinv();
		double OptimizeAllBranches(double optPrecision);
		void SetDistanceBasedBranchLengthsAroundNode(TreeNode *nd);
		void FindNearestTerminalUp(TreeNode *start, TreeNode *&, double &dist);
		void FindNearestTerminalsDown(TreeNode *start, TreeNode *from, TreeNode *&term1, TreeNode *&term2, double &dist1, double &dist2);
		void TraceLikeUpFromRoot(TreeNode *here, TreeNode *calledfrom, bool firstNode);
		void TraceLikeUpFromRootRateHet(TreeNode *here, TreeNode *calledfrom, bool firstNode);
		double LocalLike(TreeNode *here);
		double LocalLikeRateHet(TreeNode *here);
		double BranchLike(TreeNode *optNode);
		void OptimizeBranchesAroundNode(TreeNode *nd, double optPrecision, int subtreeNode);
		void OptimizeBranchesWithinRadius(TreeNode *nd, double optPrecision, int subtreeNode, TreeNode *prune=NULL);
		void OptimizeBranchesInArray(int *nodes, int numNodes, double optPrecision);
		double RecursivelyOptimizeBranches(TreeNode *nd, double optPrecision, int subtreeNode, int radius, bool dontGoNext, double scoreIncrease, bool ignoreDelta=false);
		double RecursivelyOptimizeBranchesDown(TreeNode *nd, TreeNode *calledFrom, double optPrecision, int subtreeNode, int radius, double scoreIncrease);
		double BrentOptimizeBranchLength(double accuracy_cutoff, TreeNode *here, bool firstPass);
		void CalcPartialCla(TreeNode *optNode);
		void CalcPartialClaRateHet(TreeNode *optNode);
		double SubTreeScore( TreeNode *nd);
		double SubTreeScoreRateHet( TreeNode *nd);
		int BipartitionBasedRecombination( Tree *t, bool sameModel, double optPrecision);
		int SubtreeBasedRecombination( Tree *t, int recomNodeNum, bool sameModel, double optPrecision);
		void CalcBipartitions();
		void OutputBipartitions();
		TreeNode *ContainsBipartition(const Bipartition *bip);
		TreeNode *ContainsBipartitionOrComplement(const Bipartition *bip);
		void SetAllTempClasDirty();
		void SetAppropriateTempClasDirty(int first, int sec);
		void SetTempClasDirtyWithinSubtree(int subtreeNode);
		void SetSpecifiedTempClasDirty(int *list);
		void OutputTempClaStatus(bool test=false);
		void OutputTreeStructure(TreeNode *);
		int getNumTipsTotal();
		int getNumNodesTotal() {return numNodesTotal;}
		int GetRandomInternalNode() const;
		int GetRandomTerminalNode() const;
		int GetRandomNonRootNode() const;
		void GetInternalStateString(char *string, int nodeNum);
		void RecursivelyCalculateInternalStateProbs(TreeNode *nd, ofstream &out);	
		void InferAllInternalStateProbs(char *ofprefix);
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
		bool IdenticalSubtreeTopology(const TreeNode *other);
		bool IdenticalTopology(const TreeNode *other);
		int NodeToNodeDistance(int num1, int num2);
		int NodesToRoot(TreeNode *nd);
		int GatherNodesInRadius(int range, TreeNode *sib, TreeNode *subtreeNode);
		int NumNodesTotal() {return numNodesTotal;}
		void SampleBlenCurve(TreeNode *nd, ofstream &out);
		};

inline int Tree::getNumTipsTotal()
{
	return numTipsTotal;
}

inline bool Tree::IsGood(){
	return root->IsGood();
	}
	
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
//	bool poo=true;
//	while(poo);

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
	
inline int Tree::GetRandomInternalNode() const{
	return numTipsTotal+rnd.random_int(numTipsTotal-3)+1;
	}	

inline int Tree::GetRandomTerminalNode() const{
	return rnd.random_int(numTipsTotal)+1;
	}	
	
inline int Tree::GetRandomNonRootNode() const{
	return rnd.random_int(numNodesTotal-1)+1;
	}
	
inline void Tree::SetBranchLength(TreeNode *nd, double len){
	assert(!(len < DEF_MIN_BRLEN) && !(len > DEF_MAX_BRLEN));
	nd->dlen=len;
	SweepDirtynessOverTree(nd);
	}

inline CondLikeArray *Tree::GetClaDown(TreeNode *nd, bool calc/*=true*/){
	if(claMan->IsDirty(nd->claIndexDown)){
		if(calc==true){
			if(mod->NRateCats()>1) ConditionalLikelihoodRateHet(DOWN, nd);
			else ConditionalLikelihood(DOWN, nd);
			}
		else claMan->FillHolder(nd->claIndexDown, 1);
		}
	if(memLevel > 1) claMan->ReserveCla(nd->claIndexDown);
	return claMan->GetCla(nd->claIndexDown);
	}
	
inline CondLikeArray *Tree::GetClaUpLeft(TreeNode *nd, bool calc/*=true*/){
	if(claMan->IsDirty(nd->claIndexUL)){
		if(calc==true){
			if(mod->NRateCats()>1) ConditionalLikelihoodRateHet(UPLEFT, nd);
			else ConditionalLikelihood(UPLEFT, nd);
			}
		else claMan->FillHolder(nd->claIndexUL, 2);
		}
	if(memLevel > 0) claMan->ReserveCla(nd->claIndexUL);
	return claMan->GetCla(nd->claIndexUL);
	}
	
inline CondLikeArray *Tree::GetClaUpRight(TreeNode *nd, bool calc/*=true*/){
	if(claMan->IsDirty(nd->claIndexUR)){
		if(calc==true){
			if(mod->NRateCats()>1) ConditionalLikelihoodRateHet(UPRIGHT, nd);
			else ConditionalLikelihood(UPRIGHT, nd);
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

