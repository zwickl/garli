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


#ifndef __TREE_NODE_H
#define __TREE_NODE_H
#include <iostream>
#include <vector>
#include <cassert>
using namespace std;
//#include "gaml.h"
#include "condlike.h"
#include "clamanager.h"
#include "bipartition.h"

/* GANESH included the following #includes */
#include "hashdefines.h"

class subset;
class HKYdata;

class TreeNode
{
	public:
		TreeNode* left,* right,* next,* prev,* anc;
 		int nodeNum;
 		int claIndexDown;
 		int claIndexUL;
 		int claIndexUR;
 		double dlen;
 /*		short rescaleRankDown;
 		short rescaleRankUL;
 		short rescaleRankUR;
 */		bool attached;
		Bipartition *bipart;
		//unsigned char *tipData;
		char *tipData;

#ifdef PECR_SET_PARSIMONY_BRLEN
        bool leaf_mask;		
#endif

		TreeNode( const int i = -1 );
		~TreeNode();
		
		char *MakeNewick(char *s, bool internalNodes) const;
		TreeNode * AddDes(TreeNode *);//returns argument
		void SubstituteNodeWithRespectToAnc(TreeNode *subs);
		int CountBranches(int s);
		int CountTerminals(int s);
		int CountTerminalsDown(int s, TreeNode *calledFrom);
		void CountSubtreeBranchesAndDepth(int &branches, int &sum, int depth, bool first) const;
		void MarkTerminals(int *taxtags);
		void Prune();
        TreeNode* FindNode( int &n);
        TreeNode* FindNode( int &n,TreeNode *tempno);
		bool IsGood();
		void CheckforPolytomies();
		void CountNumberofNodes(int &nnodes);
		void CheckforLeftandRight();
		void FindCrazyLongBranches();
		void FindCrazyShortBranches();
		void CheckTreeFormation();
		void CalcDepth(int &dep);
		void CopyOneClaIndex(const TreeNode *from, ClaManager *claMan, int dir);
		void AllocateMultipliers(int nchar);
		Bipartition* CalcBipartition();
		void OutputBipartition(ostream &out);
		void MarkUnattached(bool includenode);
		void RotateDescendents();
		void AdjustClasForReroot(int dir);
		void PrintSubtreeMembers(ofstream &out);
		
		void RemoveSubTreeFromSubset(subset &sub, bool startnode);
		void AddNodesToList(vector<int> &list);
		void FlipBlensToRoot(TreeNode *from);
		
		void getNNIList(vector<int> &NNIList);
		void getTaxonSwapList(vector<int> &TaxonSwapList);
		void getSPRList(int cut, vector<int> &SPRList);
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








