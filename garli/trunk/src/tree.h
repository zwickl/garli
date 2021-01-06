// GARLI version 2.1 source code
// Copyright 2005-2014 Derrick J. Zwickl
// email: garli.support@gmail.com
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
#include "model.h"
#include "sequencedata.h"
#include "reconnode.h"


#undef BRENT

class SequenceData;
class NucleotideData;
class ClaManager;
class GeneralGamlConfig;
class ModelPartition;
class Individual;
extern rng rnd;

#define RESCALE_ARRAY_LENGTH 90

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
		ModelPartition *modPart;
		TreeNode *root;
		TreeNode *dummyRoot;//when we are dummy rootinging this will just alias allNodes[numTipsTotal]
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
		static const DataPartition *dataPart;
		static FLOAT_TYPE treeRejectionThreshold;
		static vector<Constraint> constraints;
		static AttemptedSwapList attemptedSwaps;
		static FLOAT_TYPE uniqueSwapBias;
		static FLOAT_TYPE distanceSwapBias;
		static unsigned rescaleEvery;
		static FLOAT_TYPE rescaleBelow;
		static FLOAT_TYPE reduceRescaleBelow;
		static FLOAT_TYPE bailOutBelow;
		static list<TreeNode *> nodeOptVector;
		
		static bool useOptBoundedForBlen;
		static bool rootWithDummy;
		static bool dummyRootBranchMidpoint;
		static bool someOrientedGap;
		
		static FLOAT_TYPE uniqueSwapPrecalc[500];
		static FLOAT_TYPE distanceSwapPrecalc[1000];
		static FLOAT_TYPE expectedPrecision;

		static FLOAT_TYPE rescalePrecalcThresh[RESCALE_ARRAY_LENGTH];
		static FLOAT_TYPE rescalePrecalcMult[RESCALE_ARRAY_LENGTH];
		static int rescalePrecalcIncr[RESCALE_ARRAY_LENGTH];

		static Bipartition *outgroup;

		static int siteToScore;

		int calcs;

		//this controls the amount of site likelihood output. It is easier to just set it for the whole
		//tree instead of passing it around a lot.  0 = no sitelikes, 1 = user level sitelikes, 2 = debugging
		//it is NOT PERSISTENT, so after OutputSitelikes is called it is reset to 0
		int sitelikeLevel;
		string ofprefix;

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
		void AllocateTree(bool withExtraNode);
		void AssignDataToTips();
		void AssignCLAsFromMaster();
	
		//destructor
		~Tree();

		//functions for manipulating and making trees
		void RandomlyAttachTip(int nodenum, int & );
		void RandomlyAttachTipWithConstraints(int nodenum, int &placeInAllNodes, Bipartition *mask);
		void MakeTrifurcatingRoot(bool reducenodes, bool clasAssigned);
		bool ArbitrarilyBifurcate();
		void SortAllNodesArray();
		void EliminateNode(int nn);
		int FindUnusedNode(int start);
		inline void SetBranchLength(TreeNode *nd, FLOAT_TYPE len, bool dummyRootDontRecurse=false);
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
		void GenerateTopologiesAtSprDistance(Individual *source, double optPrecision, int range);
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
		FLOAT_TYPE Treelength();
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
		FLOAT_TYPE GetScorePartialTerminalRateHet(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const char *Ldata, int modIndex, int dataIndex);
		FLOAT_TYPE GetScorePartialTerminalNState(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const char *Ldata, int modIndex, int dataIndex);
		FLOAT_TYPE GetScorePartialInternalRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, int modIndex, int dataIndex);
		FLOAT_TYPE GetScorePartialInternalNState(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, int modIndex, int dataIndex);
		int Score(int rootNodeNum =0);

		FLOAT_TYPE GetScorePartialTerminalOrientedGap(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const char *Ldata, int modIndex, int dataIndex);
	
		//functions to optimize blens and params
		pair<FLOAT_TYPE, FLOAT_TYPE> CalcDerivativesRateHet(TreeNode *nd1, TreeNode *nd2);
		FLOAT_TYPE NewtonRaphsonOptimizeBranchLength(FLOAT_TYPE precision1, TreeNode *nd, bool goodGuess);
#ifdef OPT_DEBUG
		FLOAT_TYPE NewtonRaphsonSpoof(FLOAT_TYPE precision1, TreeNode *nd, bool goodGuess);
#endif
		void GetDerivsPartialTerminal(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, const char *Ldata, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot, int modIndex, int dataIndex, const unsigned *ambigMap =NULL);
		void GetDerivsPartialTerminalNState(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, const char *Ldata, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot, int modIndex, int dataIndex);
		void GetDerivsPartialTerminalNStateRateHet(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, const char *Ldata, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot, int modIndex, int dataIndex);
		void GetDerivsPartialInternal(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, FLOAT_TYPE &d1, FLOAT_TYPE &d2, int modIndex, int dataIndex);
		void GetDerivsPartialInternalNState(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, FLOAT_TYPE &d1, FLOAT_TYPE &d2, int modIndex, int dataIndex);
		void GetDerivsPartialInternalNStateRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot, int modIndex, int dataIndex);
		void GetDerivsPartialInternalEQUIV(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, FLOAT_TYPE &d1, FLOAT_TYPE &d2, char *equiv, int modIndex, int dataIndex);
		void CalcFullCLAInternalInternal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, int modIndex, int dataIndex);
		void CalcFullCLATerminalTerminal(CondLikeArray *destCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, const char *Ldata, const char *Rdata, int modIndex, int dataIndex);
		void CalcFullCLAInternalTerminal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const FLOAT_TYPE *pr1, const FLOAT_TYPE *pr2, char *data2, const unsigned *ambigMap, int modIndex, int dataIndex);
		void CalcFullCLAInternalTerminalOpenMP(CondLikeArray *destCLA, const CondLikeArray *LCLA, const FLOAT_TYPE *pr1, const FLOAT_TYPE *pr2, char *data2, const unsigned *ambigMap, int modIndex, int dataIndex);

		void CalcFullCLAInternalInternalEQUIV(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr,const char *leftEQ, const char *rightEQ, int modIndex, int dataIndex);

		void CalcFullCLAPartialInternalRateHet(CondLikeArray *destCLA, const CondLikeArray *LCLA, const FLOAT_TYPE *pr1, CondLikeArray *partialCLA, int modIndex, int dataIndex);
		void CalcFullCLAPartialTerminalRateHet(CondLikeArray *destCLA, const CondLikeArray *partialCLA, const FLOAT_TYPE *Lpr, char *Ldata, int modIndex, int dataIndex);

		void CalcFullCLAInternalInternalNState(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, int modIndex, int dataIndex);
		void CalcFullCLAInternalTerminalNState(CondLikeArray *destCLA, const CondLikeArray *LCLA, const FLOAT_TYPE *pr1, const FLOAT_TYPE *pr2, char *data2, int modIndex, int dataIndex);
		void CalcFullCLATerminalTerminalNState(CondLikeArray *destCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, const char *Ldata, const char *Rdata, int modIndex, int dataIndex);

		//for all internal state recon
		void GetStatewiseUnscaledPosteriorsPartialInternalNState(CondLikeArray *destCLA, const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, int modIndex, int dataIndex);
		void GetStatewiseUnscaledPosteriorsPartialTerminalNState(CondLikeArray *destCLA, const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const char *Ldata, int modIndex, int dataIndex);
		int FillStatewiseUnscaledPosteriors(CondLikeArraySet *partialCLAset, CondLikeArraySet *childCLAset, TreeNode *child, FLOAT_TYPE blen1);
/*		
		void CalcFullCLATerminalTerminalOrientedGap(CondLikeArray *destCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, const char *Ldata, const char *Rdata, int modIndex, int dataIndex);
		void CalcFullCLAInternalTerminalOrientedGap(CondLikeArray *destCLA, const CondLikeArray *LCLA, const FLOAT_TYPE *pr1, const FLOAT_TYPE *pr2, char *data2, int modIndex, int dataIndex);
		void CalcFullCLAInternalInternalOrientedGap(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, int modIndex, int dataIndex);
*/
		//new version that should work for any node - for each child data or a CLA will be passed in
		void CalcFullCLAOrientedGap(CondLikeArray *destCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const char *Ldata, const char *Rdata, int modIndex, int dataIndex);

		void UpdateCLAs(CondLikeArraySet *destCLA, CondLikeArraySet *firstCLA, CondLikeArraySet *secCLA, TreeNode *firstChild, TreeNode *secChild, FLOAT_TYPE blen1, FLOAT_TYPE blen2);
		void GetTotalScore(CondLikeArraySet *partialCLA, CondLikeArraySet *childCLA, TreeNode *child, FLOAT_TYPE blen1);

		FLOAT_TYPE OptimizeBranchLength(FLOAT_TYPE optPrecision, TreeNode *nd, bool goodGuess);
		FLOAT_TYPE OptimizeAllBranches(FLOAT_TYPE optPrecision);
		int PushBranchlengthsToMin();
		void OptimizeBranchesAroundNode(TreeNode *nd, FLOAT_TYPE optPrecision, int subtreeNode);
		void OptimizeBranchesWithinRadius(TreeNode *nd, FLOAT_TYPE optPrecision, int subtreeNode, TreeNode *prune);
		void OptimizeBranchesInArray(int *nodes, int numNodes, FLOAT_TYPE optPrecision);
		FLOAT_TYPE RecursivelyOptimizeBranches(TreeNode *nd, FLOAT_TYPE optPrecision, int subtreeNode, int radius, bool dontGoNext, FLOAT_TYPE scoreIncrease, bool ignoreDelta=false);
		FLOAT_TYPE RecursivelyOptimizeBranchesDown(TreeNode *nd, TreeNode *calledFrom, FLOAT_TYPE optPrecision, int subtreeNode, int radius, FLOAT_TYPE scoreIncrease);
		FLOAT_TYPE BrentOptimizeBranchLength(FLOAT_TYPE accuracy_cutoff, TreeNode *here, bool firstPass);
		FLOAT_TYPE BranchLike(TreeNode *optNode);
		FLOAT_TYPE OptimizeAlpha(FLOAT_TYPE, int modnum);
		FLOAT_TYPE OptimizeOmegaParameters(FLOAT_TYPE prec, int modnum);
		FLOAT_TYPE OptimizeFlexRates(FLOAT_TYPE prec, int modnum);
		FLOAT_TYPE OptimizeEquilibriumFreqs(FLOAT_TYPE prec, int modnum);
		FLOAT_TYPE OptimizeRelativeNucRates(FLOAT_TYPE prec, int modnum);
		FLOAT_TYPE OptimizeInsertDeleteRates(FLOAT_TYPE prec, int modnum);
		FLOAT_TYPE OptimizeSubsetRates(FLOAT_TYPE prec);
//the new versions from the trunk
#ifdef SINGLE_PRECISION_FLOATS
		FLOAT_TYPE OptimizeBoundedParameter(int modnum, FLOAT_TYPE optPrecision, FLOAT_TYPE initialVal, int which, FLOAT_TYPE lowBound, FLOAT_TYPE highBound, void (Model::*SetParam)(int, FLOAT_TYPE), FLOAT_TYPE targetScoreDigits = 5.0f);
#else
		FLOAT_TYPE OptimizeBoundedParameter(int modnum, FLOAT_TYPE optPrecision, FLOAT_TYPE initialVal, int which, FLOAT_TYPE lowBound, FLOAT_TYPE highBound, void (Model::*SetParam)(int, FLOAT_TYPE), FLOAT_TYPE targetScoreDigits = 9.0);
#endif
		FLOAT_TYPE OptimizeBoundedParameter(FLOAT_TYPE optPrecision, FLOAT_TYPE prevVal, int which, FLOAT_TYPE lowBound, FLOAT_TYPE highBound, int modnum, void (Model::*SetParam)(int, FLOAT_TYPE));
		template<class T> FLOAT_TYPE OptimizeBoundedParameter(FLOAT_TYPE optPrecision, FLOAT_TYPE prevVal, int which, FLOAT_TYPE lowBound, FLOAT_TYPE highBound, T *obj, void (T::*SetParam)(int, FLOAT_TYPE));
		template<class T> void TraceParameterLikelihood(ofstream &out, int which, FLOAT_TYPE prevVal, FLOAT_TYPE startVal, FLOAT_TYPE endVal, FLOAT_TYPE incr, T *obj, void (T::*SetParam)(int, FLOAT_TYPE));

		void TraceLikelihoodForParameter(int modnum, int which, FLOAT_TYPE init, FLOAT_TYPE min, FLOAT_TYPE max, FLOAT_TYPE interval, void (Model::*SetParam)(int, FLOAT_TYPE), bool append);
		//FLOAT_TYPE OptimizeBoundedParameter(FLOAT_TYPE optPrecision, FLOAT_TYPE prevVal, int which, FLOAT_TYPE lowBound, FLOAT_TYPE highBound, void (Model::*SetParam)(int, FLOAT_TYPE));
		FLOAT_TYPE SetAndEvaluateParameter(int modnum, int which, FLOAT_TYPE val, FLOAT_TYPE &bestKnownScore, FLOAT_TYPE &bestKnownVal, void (Model::*SetParam)(int, FLOAT_TYPE));
		bool CheckScoreAndRestore(int modnum, int which, void (Model::*SetParam)(int, FLOAT_TYPE), FLOAT_TYPE otherScore, FLOAT_TYPE otherVal, FLOAT_TYPE bestScore, FLOAT_TYPE bestVal, FLOAT_TYPE tolerance);

		FLOAT_TYPE OptimizeTreeScale(FLOAT_TYPE);
		FLOAT_TYPE OptimizePinv();
		void SetNodesUnoptimized();
		void RescaleRateHet(CondLikeArray *destCLA, int dataIndex);
		void RescaleRateHetNState(CondLikeArray *destCLA, int dataIndex);

		void StoreBranchlengths(vector<FLOAT_TYPE> &blens){
			for(int n=1;n<numNodesTotal;n++)
				blens.push_back(allNodes[n]->dlen);
			assert(blens.size() == numNodesTotal - 1);
			}
		void RestoreBranchlengths(vector<FLOAT_TYPE> &blens){
			for(int n=1;n<numNodesTotal;n++)
				SetBranchLength(allNodes[n], blens[n-1]);
			MakeAllNodesDirty();
			}

		void VerifyScore(double tol){
			if(!FloatingPointEquals(-1.0, lnL, 1e-8)){
				double bef = lnL;
				MakeAllNodesDirty();
				Score();
				outman.UserMessage("Verify score: %.6f -> %.6f = %.6f", bef, lnL, bef - lnL);
				assert(FloatingPointEquals(bef, lnL, tol));
				}
			return;
			}

		pair<FLOAT_TYPE, FLOAT_TYPE> OptimizeSingleSiteTreeScale(FLOAT_TYPE optPrecision);

		//functions for dealing with conditional likelihood arrays
		void MarkUpwardClasToReclaim(int subtreeNode);
		void MarkDownwardClasToReclaim(int subtreeNode);
		void MarkClasNearTipsToReclaim(int subtreeNode);
		void ProtectClas();
		void UnprotectClas();
		inline CondLikeArraySet *GetClaDown(TreeNode *nd, bool calc=true);
		inline CondLikeArraySet *GetClaUpLeft(TreeNode *nd, bool calc=true);
		inline CondLikeArraySet *GetClaUpRight(TreeNode *nd, bool calc=true);
		void OutputValidClaIndeces();
		void OutputNthClaAcrossTree(ofstream &deb, TreeNode *nd, int site, int modIndex);
		void ClaReport(ofstream &cla);
		FLOAT_TYPE CountClasInUse();
		void OutputSiteLikelihoods(int partnum, vector<FLOAT_TYPE> &likes, const int *under1, const int *under2);
		void OutputSiteDerivatives(int partNum, vector<double> &likes, vector<double> &d1s, vector<double> &d2s, const int *under1, const int *under2, ofstream &ordered, ofstream &packed);
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
		int NTax() const {return numTipsTotal;}
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
		void MoveDummyRootToBranchMidpoint();

		static void SetTreeStatics(ClaManager *, const DataPartition *, const GeneralGamlConfig *);

		void C4(const FLOAT_TYPE *a);
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
	if(claMan == NULL)
		return;
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
	
inline void Tree::SetBranchLength(TreeNode *nd, FLOAT_TYPE len, bool dummyRootDontRecurse /*=false*/){
	assert(!(len < min_brlen) && !(len > max_brlen));
	nd->dlen=len;

	//the dontRecurse bit here just keeps it from bouncing back and forth setting the 
	//lengths of the two "root" branches, since changing one triggers a change to the other
	//There are are few posibilities:
	//1. nd->anc is the root and dummyRoot->anc is the root
	//		The other branch to adjust is the descendent of the root that is not the dummyRoot nor nd
	//2. nd->anc is not the root
	//	2a. nd is the branch "above" where dummyRoot attaches
	//		In this case dummyRoot is the next or prev of nd, and the other branch to adjust is nd->anc
	//	2b. nd is the branch "below" where dummyRoot attaches
	//		In this case dummyRoot is left or right of nd.  The branch to adjust is the other descendent of nd
	//	3.	nd is not related to the dummyRooted branch
	if(rootWithDummy && dummyRootBranchMidpoint && dummyRootDontRecurse == false){
		TreeNode *otherNode = NULL;
		if(nd->anc == root && dummyRoot->anc == root){
			otherNode = root->left;
			do{
				if(otherNode != dummyRoot && otherNode != nd){
					break;
					}
				else
					otherNode = otherNode->next;
				}while(otherNode);
			}
		else{
			if(nd->prev == dummyRoot || nd->next == dummyRoot)
				otherNode = nd->anc;
			else{
				if(nd->left == dummyRoot)
					otherNode = nd->right;
				else if(nd->right == dummyRoot)
					otherNode = nd->left;
				}
			}
		if(otherNode)
			SetBranchLength(otherNode, len, true);
		}

	SweepDirtynessOverTree(nd);
	}

inline void Tree::MoveDummyRootToBranchMidpoint(){
	TreeNode *branch1, *branch2;
	if(dummyRoot->anc == root){
		if(root->left != dummyRoot){
			branch1 = root->left;
			if(root->left->next != dummyRoot){
				branch2 = root->left->next;
				assert(root->right == dummyRoot);
				}
			else{
				branch2 = root->right;
				}
			}
		else{
			branch1 = root->left->next;
			branch2 = root->right;
			}
		}
	else{
		branch1 = dummyRoot->anc;
		branch2 = (dummyRoot->next ? dummyRoot->next : dummyRoot->prev);
		}
	double sum = branch1->dlen + branch2->dlen;
	//this should automatically adjust the length of branch2 because of code in SetBranchLength
	SetBranchLength(branch1, sum / 2.0);
	}

inline CondLikeArraySet *Tree::GetClaDown(TreeNode *nd, bool calc/*=true*/){
	if(claMan->IsDirty(nd->claIndexDown)){
		if(calc==true){
			ConditionalLikelihoodRateHet(DOWN, nd);
			}
		else claMan->FillHolder(nd->claIndexDown, 1);
		}
	if(memLevel > 1) claMan->ReserveCla(nd->claIndexDown);
	return claMan->GetCla(nd->claIndexDown);
	}
	
inline CondLikeArraySet *Tree::GetClaUpLeft(TreeNode *nd, bool calc/*=true*/){
	if(claMan->IsDirty(nd->claIndexUL)){
		if(calc==true){
			ConditionalLikelihoodRateHet(UPLEFT, nd);
			}
		else claMan->FillHolder(nd->claIndexUL, 2);
		}
	if(memLevel > 0) claMan->ReserveCla(nd->claIndexUL);
	return claMan->GetCla(nd->claIndexUL);
	}
	
inline CondLikeArraySet *Tree::GetClaUpRight(TreeNode *nd, bool calc/*=true*/){
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

template<class T>
void Tree::TraceParameterLikelihood(ofstream &out, int which, FLOAT_TYPE prevVal, FLOAT_TYPE startVal, FLOAT_TYPE endVal, FLOAT_TYPE incr, T *obj, void (T::*SetParam)(int, FLOAT_TYPE)){
	for(FLOAT_TYPE val = startVal;val <=endVal;val += incr){
		CALL_SET_PARAM_FUNCTION(*obj, SetParam)(which, val);
		MakeAllNodesDirty();
		Score();
		out << val << "\t" << lnL << endl;
		}
	CALL_SET_PARAM_FUNCTION(*obj, SetParam)(which, prevVal);
	MakeAllNodesDirty();
	Score();
	out << prevVal << "\t" << lnL << endl;
	}

//a templated version
template<class T>
FLOAT_TYPE Tree::OptimizeBoundedParameter(FLOAT_TYPE optPrecision, FLOAT_TYPE prevVal, int which, FLOAT_TYPE lowBound, FLOAT_TYPE highBound, T *obj, void (T::*SetParam)(int, FLOAT_TYPE)){

	FLOAT_TYPE epsilon = min(optPrecision, (FLOAT_TYPE) 1.0e-5);

	assert(prevVal > lowBound - epsilon && prevVal < highBound + epsilon);

	//if the initial value is already very near or equal to a bound, bump it off a tad so the emirical derivs below work right.
	if(prevVal - lowBound < epsilon){
		prevVal = lowBound + epsilon;
		CALL_SET_PARAM_FUNCTION(*obj, SetParam)(which, prevVal);
		MakeAllNodesDirty();
		}
	else if(highBound - prevVal < epsilon){
		prevVal = highBound - epsilon;
		CALL_SET_PARAM_FUNCTION(*obj, SetParam)(which, prevVal);
		MakeAllNodesDirty();
		}

	if(FloatingPointEquals(lnL, -ONE_POINT_ZERO, 1e-8)) Score();
	FLOAT_TYPE start, prev, cur;
	prev = start = cur = lnL;
	FLOAT_TYPE lastChange=(FLOAT_TYPE)9999.9;
	FLOAT_TYPE upperBracket = highBound;   //the smallest value we know of with a negative d1, or the minimum allowed value
	FLOAT_TYPE lowerBracket = lowBound;   //the largest value we know of with a positive d1 , or the maximum allowed value
	FLOAT_TYPE incr;
	int lowBoundOvershoot = 0;
	int upperBoundOvershoot = 0;
	int positiveD2Num = 0;
	int pass = 0;

/*	ofstream curves("lcurve.log");
	curves.precision(8);
	curves << endl;
	for(double c = max(prevVal - 0.01, lowBound); c < min(prevVal + 0.01, highBound) ; c += 0.001){
		CALL_SET_PARAM_FUNCTION(*obj, SetParam)(which, c);
		MakeAllNodesDirty();
		Score();
		curves << c << "\t" << lnL << endl;;
		}
	curves.close();

	CALL_SET_PARAM_FUNCTION(*obj, SetParam)(which, prevVal);
	MakeAllNodesDirty();
	Score();
*/
	while(1){
#ifdef SINGLE_PRECISION_FLOATS
		incr=0.005f;
#else
		incr=min(0.0001*optPrecision, min(prevVal - lowerBracket, upperBracket - prevVal));
		//incr=min(0.0001, min(prevVal - lowerBracket, upperBracket - prevVal));
#endif
		CALL_SET_PARAM_FUNCTION(*obj, SetParam)(which, prevVal+incr);
		MakeAllNodesDirty();
		Score();
		cur=lnL;
		FLOAT_TYPE d11=(cur-prev)/incr;
		
		CALL_SET_PARAM_FUNCTION(*obj, SetParam)(which, prevVal-incr);
		MakeAllNodesDirty();
		Score();
		cur=lnL;
		FLOAT_TYPE d12=(cur-prev)/-incr;
		
		FLOAT_TYPE d1=(d11+d12)*ZERO_POINT_FIVE;
		//if the evaluation points straddle the optimum, leave now
		if((d11 - d12) == ZERO_POINT_ZERO || (d11 > ZERO_POINT_ZERO && d12 < ZERO_POINT_ZERO) || (d11 < ZERO_POINT_ZERO && d12 > ZERO_POINT_ZERO)){
			CALL_SET_PARAM_FUNCTION(*obj, SetParam)(which, prevVal);;
			MakeAllNodesDirty();
			lnL = prev;
			return prev-start;
			}
		
		FLOAT_TYPE d2=(d11-d12)/incr;
		
		FLOAT_TYPE est=-d1/d2;
		
		FLOAT_TYPE proposed = prevVal + est;

	//	outman.UserMessage("%f\t%f\t%f\t%f\t%f", d1, d2, prevVal, est, proposed);

		if(d2 > ZERO_POINT_ZERO){
			positiveD2Num++;
			if(d1 > ZERO_POINT_ZERO) proposed=prevVal*(FLOAT_TYPE)(ONE_POINT_ZERO+0.01*positiveD2Num);
			else proposed=prevVal*(FLOAT_TYPE)(ONE_POINT_ZERO-0.01*positiveD2Num);
			}
		if(d1 < ZERO_POINT_ZERO && proposed < (lowerBracket + epsilon)){
			if(prevVal - lowerBracket - epsilon < epsilon * ZERO_POINT_FIVE){
				CALL_SET_PARAM_FUNCTION(*obj, SetParam)(which, prevVal);;
				MakeAllNodesDirty();			
				lnL = prev;
				return prev-start;
				}
			lowBoundOvershoot++;
			if(lowBoundOvershoot > 1)
				proposed = lowerBracket + epsilon;
			else
				proposed = (prevVal + lowerBracket) * ZERO_POINT_FIVE;
			}
		else if(d1 > ZERO_POINT_ZERO && proposed > upperBracket - epsilon){
			if(upperBracket - epsilon - prevVal < epsilon * ZERO_POINT_FIVE){
				CALL_SET_PARAM_FUNCTION(*obj, SetParam)(which, prevVal);;
				MakeAllNodesDirty();
				lnL = prev;
				return prev-start;
				}
			upperBoundOvershoot++;
			if(upperBoundOvershoot > 1)
				proposed = upperBracket - epsilon;
			else
				proposed = (prevVal + upperBracket) * ZERO_POINT_FIVE;
			}

		FLOAT_TYPE estImprove;
		if(d2 < ZERO_POINT_ZERO) estImprove = d1*(proposed - prevVal) + (d2 * (proposed - prevVal) * (proposed - prevVal)) * ZERO_POINT_FIVE;
		else estImprove = 9999.9;

		//require that we didn't significantly worsen the likelihood 
		if(estImprove < optPrecision && prev >= start - 1.0e-6){
			CALL_SET_PARAM_FUNCTION(*obj, SetParam)(which, prevVal);;
			MakeAllNodesDirty();			
			lnL = prev;
			return prev-start;
			}

		//don't allow infinite looping if something goes wrong
		if(pass > 1000){
			throw ErrorException("too many passes in OptimizeBoundedParameter");
			}
		
		//update the brackets
		if(d1 <= ZERO_POINT_ZERO && prevVal < upperBracket)
			upperBracket = prevVal;
		else if(d1 > ZERO_POINT_ZERO && prevVal > lowerBracket)
			lowerBracket = prevVal;

		CALL_SET_PARAM_FUNCTION(*obj, SetParam)(which, proposed);;
		MakeAllNodesDirty();			
		Score();
		lastChange = lnL - prev;
		prev=lnL;
		prevVal=proposed;
		pass++;
		}
	return -1;
	}

#endif

