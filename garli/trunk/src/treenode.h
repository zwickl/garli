// GARLI version 0.96b4 source code
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


#ifndef __TREE_NODE_H
#define __TREE_NODE_H
#include <iostream>
#include <vector>
#include <cassert>
using namespace std;

#include "condlike.h"
#include "clamanager.h"
#include "bipartition.h"

/* GANESH included the following #includes */
#include "hashdefines.h"

class NucleotideData;
class MFILE;

class TreeNode{
	public:
		TreeNode* left,* right,* next,* prev,* anc;
 		int nodeNum;
 		int claIndexDown;
 		int claIndexUL;
 		int claIndexUR;
 		FLOAT_TYPE dlen;
		bool attached;
		bool alreadyOptimized;
		Bipartition *bipart;
		char *tipData;
#ifdef OPEN_MP
		unsigned *ambigMap;
#endif

		TreeNode( const int i = -1 );
		~TreeNode();
		
		//functions for manipulating nodes within a tree
		TreeNode * AddDes(TreeNode *);//returns argument
		void RemoveDes(TreeNode *d);
		void MoveDesToAnc(TreeNode *d);
		void SubstituteNodeWithRespectToAnc(TreeNode *subs);
		int CountBranches(int s);
		int CountTerminals(int s);
		int CountTerminalsDown(int s, TreeNode *calledFrom);
		void CountSubtreeBranchesAndDepth(int &branches, int &sum, int depth, bool first) const;
		void MarkTerminals(int *taxtags);
		void Prune();
        TreeNode* FindNode( int &n);
        TreeNode* FindNode( int &n,TreeNode *tempno);
		void CountNumberofNodes(int &nnodes);
		void MarkUnattached(bool includenode);
		void RotateDescendents();
		void AdjustClasForReroot(int dir);
		void AddNodesToList(vector<int> &list);
		void FlipBlensToRoot(TreeNode *from);
		void FlipBlensToNode(TreeNode *from, TreeNode *stopNode);
		void RecursivelyAddOrRemoveSubtreeFromBipartitions(Bipartition *);

		//misc functions
		char *MakeNewick(char *s, bool internalNodes, bool branchLengths, bool highPrec=false) const;
		void MakeNewickForSubtree(char *s) const;
		bool IsGood();
		bool IsTerminal() const{
			return left == NULL;
			}
		bool IsInternal() const{
			return left != NULL;
			}
		bool IsRoot() const{
			return anc==NULL;
			}
		bool IsNotRoot() const{
			return anc!=NULL;
			}
		void CalcDepth(int &dep);
		void CopyOneClaIndex(const TreeNode *from, ClaManager *claMan, int dir);
		Bipartition* CalcBipartition();
		void StandardizeBipartition();
		void GatherConstrainedBiparitions(vector<Bipartition> &biparts);
		void OutputBipartition(ostream &out);
		void PrintSubtreeMembers(ofstream &out);
		void SetUnoptimized(){
			alreadyOptimized=false;
			if(left) left->SetUnoptimized();
			if(right) right->SetUnoptimized();
			}
		void SetEquivalentConditionalVectors(const SequenceData *data);
		void OutputBinaryNodeInfo(OUTPUT_CLASS &out) const;
		     
		//debugging functions for checking tree formation
		void CheckforPolytomies();
		void CheckforLeftandRight();
		void FindCrazyLongBranches();
		void FindCrazyShortBranches();
		void CheckTreeFormation();
};

inline void TreeNode::CopyOneClaIndex(const TreeNode *from, ClaManager *claMan, int dir){

	const int *indexF;
	int *indexT;
	if(dir==1){
		indexF=&from->claIndexDown;
		indexT=&claIndexDown;
		}
	else if(dir==2){
		indexF=&from->claIndexUL;
		indexT=&claIndexUL;
		}
	else if(dir==3){
		indexF=&from->claIndexUR;
		indexT=&claIndexUR;
		}
	else assert(0);
	
	claMan->DecrementCla(*indexT);
	claMan->IncrementCla(*indexF);
	*indexT=*indexF;
	}



#endif








