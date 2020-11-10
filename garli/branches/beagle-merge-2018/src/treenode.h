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


#ifndef __TREE_NODE_H
#define __TREE_NODE_H
#include <iostream>
#include <vector>
#include <cassert>
using namespace std;

#include "condlike.h"
#include "clamanager.h"
#include "bipartition.h"
#ifdef USE_BEAGLE
#include "calculationmanager.h"
#endif

class NucleotideData;
class MFILE;

class TreeNode{
	public:
		TreeNode* left,* right,* next,* prev,* anc;
 		int nodeNum;
#ifdef USE_BEAGLE
		NodeClaManager myMan;
#endif
 		int claIndexDown;
 		int claIndexUL;
 		int claIndexUR;
 		FLOAT_TYPE dlen;
		bool attached;
		bool alreadyOptimized;
		Bipartition *bipart;
		vector<char *> tipData;
#ifdef OPEN_MP
		//unsigned *ambigMap;
		vector<unsigned *> ambigMap;
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
		void RecursivelyAddOrRemoveSubtreeFromBipartitions(const Bipartition &subtree);
		void CollapseMinLengthBranches(int &);
		//misc functions
		char *MakeNewick(char *s, bool internalNodes, bool branchLengths, bool highPrec=false) const;
		void MakeNewick(string &outStr, const DataPartition *data, bool internalNodes, bool branchLengths, bool taxonNames = false, bool highPrec = false) const;
		void MakeNewickForSubtree(char *s) const;
		void MakeNewickForSubtree(string &s, const DataPartition *data, bool internalNodes, bool branchLengths, bool taxonNames = false, bool highPrec = false) const;
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
		Bipartition* CalcBipartition(bool standardize);
		Bipartition* VerifyBipartition(bool standardize);
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
		void OutputNodeConnections();
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








