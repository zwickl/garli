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

#include <cassert>
#include <iostream>
#include <fstream>

using namespace std;

#include "treenode.h"
#include "clamanager.h"
#include "bipartition.h"
#include "subset.h"

#undef DEBUG_RECOMBINEWITH

TreeNode::TreeNode( const int no )
	: left(0), right(0), next(0), prev(0), anc(0), tipData(0L), bipart(0L)
{	attached =false;
	
	claIndexDown=-1;
	claIndexUL=-1;
	claIndexUR=-1;
		
	nodeNum  = no;
	dlen	= 0.0;
/* GANESH added this */
#ifdef PECR_SET_PARSIMONY_BRLEN
    /* every node is a leaf (no descendants) when created */
    leaf_mask = true;
#endif    
}

TreeNode::~TreeNode(){
	if(bipart!=NULL) delete bipart;
}

TreeNode* TreeNode::AddDes(TreeNode *d){
	d->anc=this;
	d->next=NULL;
	if(left)
		{
		if(right){
			right->next=d;
			d->prev=right;
			right=d;
			}
		else{
			right=d;
			left->next=d;
			d->prev=left;
			}
		}
	else
		{left=d;
		d->prev=NULL;
		}
//	right=d;
	d->attached=true;

/* GANESH added this */
#ifdef PECR_SET_PARSIMONY_BRLEN
    /* not a leaf any more once we add descendants */
    leaf_mask = false;
#endif
	return d;

}
char *TreeNode::MakeNewick(char *s, bool internalNodes, bool branchLengths /*=true*/) const{
	if(left){
		if(internalNodes==true && nodeNum!=0){
			sprintf(s, "%d", nodeNum);
			while(*s)s++;
			}
		*s++='(';
		s=left->MakeNewick(s, internalNodes, branchLengths);
		if(anc){
			if(branchLengths==true){
				*s++=':';
				sprintf(s, "%.8lf", dlen);
	//			strcpy(s, ".01");
				while(*s)s++;
				}
			}
		else
			{*s='\0';
			return s;
			}
		}
	else {
		sprintf(s, "%d", nodeNum);
		while(*s)s++;
		if(branchLengths==true){
			*s++=':';
			sprintf(s, "%.8lf", dlen);
	//		strcpy(s, ".01");
			while(*s)s++;
			}
		}
		
	if(next){
		*s++=',';
		s=next->MakeNewick(s, internalNodes, branchLengths);
		}
	else {
		if(anc){
			*s++=')';
			}
		}
	return s;
	}

void TreeNode::MakeNewickForSubtree(char *s) const{
	assert(left);
	*s++='(';
	s=left->MakeNewick(s, false, false);
	*s++=';';
	*s++='\0';
	}

//MTH
void TreeNode::Prune()
{	//DZ 7-6 removing adjustments to branch lengths when pruning, which just result in adding a whole bunch of length
	//to the whole tree, making the dlens get longer and longer and longer as the run progresses
	assert(anc);//never call with this=root
	attached=false;
	if(anc->anc)
		{//not connected to the root
		if(anc->left->next==anc->right)
			{TreeNode *sis;
			if(anc->left==this)
				sis=anc->right;
			else
				sis=anc->left;
//			sis->dlen+=anc->dlen;
			anc->SubstituteNodeWithRespectToAnc(sis);
			anc->attached=false;
			}
		else
			{
			anc=anc;
			assert(0);//internal polytomy
			}
		}
	else
		{
		//these assume a trifurcating root
		if(anc->left==this){
			anc->left=next;
			anc->left->prev=NULL;
			anc->left->next=anc->right;
			}
		else
			if(anc->right==this){
				anc->right=prev;
				anc->right->next=NULL;
				anc->left->next=anc->right;
				}
			else 
				{assert(anc->left==prev && anc->right==next);
		next->prev=prev;
		prev->next=next;
		assert(anc->left &&  anc->right);
		TreeNode *temp;
				if(anc->left->left){
					//anc->right->dlen+=anc->left->dlen;
			temp=anc->left;
			temp->SubstituteNodeWithRespectToAnc(temp->left);
			anc->AddDes(temp->right);
			temp->attached=false;
			}
		else
					if(anc->right->left){
						//anc->left->dlen+=anc->right->dlen;
				temp=anc->right;
				temp->SubstituteNodeWithRespectToAnc(temp->left);
				anc->AddDes(temp->right);
				temp->attached=false;
						}
				}
		
		}
}

//MTH	
void TreeNode::SubstituteNodeWithRespectToAnc(TreeNode *subs)//note THIS DOESN't do anything with numBranchesAdded or any other tree fields that describe the tree
{		//this function moves subs into the place in the tree that had been occupied by this
		// it is called in swapping and addRandomNode and can't be called with the root as this
		//  Nothing is done with branch lengths OR LEFT OR RIGHT (this and subs still keep their descendants)
		subs->anc=anc;
		subs->prev=prev;
		subs->next=next;
		assert(anc);
		if(anc->left==this)
			anc->left=subs;
		if(anc->right==this)
			anc->right=subs;
		if(next)
			next->prev=subs;
		if(prev)
			prev->next=subs;
		subs->attached=true;
		attached=false;
		next=prev=anc=NULL;
}

int TreeNode::CountBranches(int s)
{
	if(left)
		s=left->CountBranches(++s);
	if(nodeNum==0)
	  s=left->next->CountBranches(++s);
	if(right)
		s=right->CountBranches(++s);
	return s;
}

int TreeNode::CountTerminals(int s){
	if(left)
		s=left->CountTerminals(s);
	else s++;
	if(nodeNum==0)
		s=left->next->CountTerminals(s);
	if(right)
		s=right->CountTerminals(s);
	return s;
}

int TreeNode::CountTerminalsDown(int s, TreeNode *calledFrom){
	
	TreeNode *sib;
	
	if(nodeNum!=0){		
		if(left==calledFrom) sib=right;
		else sib=left;

		if(sib)
			s=sib->CountTerminals(s);
		else s++;

		s=anc->CountTerminalsDown(s, this);
		}
	else {
		if(left!=calledFrom) s=left->CountTerminals(s);
		if(left->next!=calledFrom) s=left->next->CountTerminals(s);
		if(right!=calledFrom) s=right->CountTerminals(s);
		}
		
		
	return s;
}

void TreeNode::CountSubtreeBranchesAndDepth(int &branches, int &sum, int depth, bool first) const
//this is the version to use if you want to be
//sure not to just to another subtree (ie, don't 
//go ->next from the calling node)
{
	if(left){
		sum+=depth;
		left->CountSubtreeBranchesAndDepth(++branches, sum, depth+1, false);
		}
	if(next&&!first){
		sum+=depth-1;
		next->CountSubtreeBranchesAndDepth(++branches, sum, depth, false);
		}
}
	

void TreeNode::CalcDepth(int &dep){
	dep++;
	int l=0, r=0;
	if(left){
		left->CalcDepth(l);
		}
	if(right){
		right->CalcDepth(r);
		}
	dep += (r > l ? r : l);
	}

void TreeNode::MarkTerminals(int *taxtags)
{	if(left)
		left->MarkTerminals(taxtags);
	else
		taxtags[nodeNum]=1;
	if(next)
		next->MarkTerminals(taxtags);
}

void TreeNode::MarkUnattached(bool includenode){
	attached=false;
	if(left)
		left->MarkUnattached(false);
	if(next&&includenode==false)
		next->MarkUnattached(false);	
	}
	
TreeNode* TreeNode::FindNode( int &n, TreeNode *tempno){		
		//my version DZ.  It returns nodeNum n
		if(left&&tempno!=NULL){
			tempno=left->FindNode(n, tempno);
			}
		if(next&&tempno!=NULL){
			tempno=left->FindNode(n, tempno);
			}
		if(nodeNum==n){
			tempno=this;
			}			
		return tempno;
	}
	
//MTH
TreeNode* TreeNode::FindNode( int &n){
		//note that this function does NOT look for the node with nodeNum n, but rather
		//counts nodes and returns the nth one that it finds
		n--;
		if(n<0)
			return this;
		if(left)
			{TreeNode* tempno;
			tempno=left->FindNode(n);
			if(tempno)
				return tempno;
			}
		if(next)
			return next->FindNode(n);
		return NULL;
}

bool TreeNode::IsGood()
{	if(attached || !anc)
		{if(!left && right)
			return false;
		if(!anc)
			{if(nodeNum!=0 || next || prev)
				return false;
			}
		else
			{TreeNode *tempno;
			tempno=anc->left;
			bool found=false;
			int nsibs=0;
			while(tempno)
				{if(tempno->anc!=anc)
					return false;
				if(tempno==this)
					found=true;
				tempno=tempno->next;
				nsibs++;
				if(nsibs>3)
					return false;
				}
			if(!found)
				return false;
			}
		if(left){
			
			if(!left->IsGood())
				return false;
			}
		if(next)
			return next->IsGood();
		return true;
		}
	else
		return false;
}

void TreeNode::CountNumberofNodes(int &nnodes){
	if(left!=NULL){
		left->CountNumberofNodes(nnodes);
		}
	if(next!=NULL){
		next->CountNumberofNodes(nnodes);
		}
	nnodes++;
	}
	
void TreeNode::CheckforLeftandRight(){
	if(left!=NULL){
		left->CheckforLeftandRight();
		}
	
	if(next!=NULL){
		next->CheckforLeftandRight();
		}
	
	if((left&&!right)||(right&&!left)){
		assert(0);
		}
	}
	
void TreeNode::FindCrazyLongBranches(){
	if(left!=NULL){
		left->FindCrazyLongBranches();
		}
	
	if(next!=NULL){
		next->FindCrazyLongBranches();
		}
	if(dlen>1.0){
		cout << "WTF?\n";
		}
	}
			
void TreeNode::FindCrazyShortBranches(){
	if(left!=NULL){
		left->FindCrazyShortBranches();
		}
	
	if(next!=NULL){
		next->FindCrazyShortBranches();
}
	if(anc&&dlen<.0001){
		cout << "WTF?\n";
		}
	}
	
void TreeNode::CheckTreeFormation()	{
	//make sure that nodes that this node points to also point back (ie this->ldes->anc=this)
	if(left){
		assert(left->anc==this);
		left->CheckTreeFormation();
		}
	if(right){
		assert(right->anc==this);
		}
	if(next){
		assert(next->prev==this);
		next->CheckTreeFormation();
		}
	if(prev){
		assert(prev->next==this);
		}
	assert(!anc||dlen>0.0);
	}

void TreeNode::CheckforPolytomies(){
	
	if(left!=NULL){
		left->CheckforPolytomies();
		}

	if(next!=NULL){
		next->CheckforPolytomies();
		}
		
	if(anc!=NULL){
		if(left!=NULL){
			if(left->next!=right){
				assert(0);
				}
			}
		}
	}
/*	
void TreeNode::AllocateMultipliers(int nchar){	
	underflow_mult=new double[nchar];
	for(int i=0;i<nchar;i++) underflow_mult[i]=0.0;
	}
*/

Bipartition* TreeNode::CalcBipartition(){	
	if(left&&anc){
		*bipart=left->CalcBipartition();
		*bipart+=left->next->CalcBipartition();
		return bipart;
		}
	else if(anc){
		return bipart->TerminalBipart(nodeNum);	
		}
	else{
		left->CalcBipartition();
		left->next->CalcBipartition();
		left->next->next->CalcBipartition();
		}
	return NULL;
	}

void TreeNode::OutputBipartition(ostream &out){	
	if(left&&anc){
		left->OutputBipartition(out);
		left->next->OutputBipartition(out);
		out << bipart->Output() << endl;
		}
	else if(!anc){
		left->OutputBipartition(out);
		left->next->OutputBipartition(out);
		left->next->next->OutputBipartition(out);
		}
	}

void TreeNode::RotateDescendents(){
	//don't call this with the root!
	assert(anc);
	TreeNode* tmp=right;
	right=left;
	left=tmp;
	left->prev=NULL;	
	left->next=right;
	right->next=NULL;
	}
	
void TreeNode::RemoveSubTreeFromSubset(subset &sub, bool startnode){
	sub.elementremove(nodeNum);
	if(left != NULL) left->RemoveSubTreeFromSubset(sub, false);
	if(startnode==false && next != NULL) next->RemoveSubTreeFromSubset(sub, false);
	}

void TreeNode::AddNodesToList(vector<int> &list){
	list.push_back(nodeNum);
	if(left!=NULL) left->AddNodesToList(list);
	if(next!=NULL) next->AddNodesToList(list);
	}
	
void TreeNode::FlipBlensToRoot(TreeNode *from){
	if(anc!=NULL) anc->FlipBlensToRoot(this);
	if(from==NULL) dlen=-1;
	else dlen=from->dlen;
	}

void TreeNode::PrintSubtreeMembers(ofstream &out){
	if(left==NULL) out << nodeNum << "\t"; 
	else left->PrintSubtreeMembers(out);
	if(next!=NULL) next->PrintSubtreeMembers(out);
	}
	
void TreeNode::AdjustClasForReroot(int dir){
	int tmp=claIndexDown;
	if(dir==2){//the ancestor and left des have been swapped
		claIndexDown=claIndexUL;
		claIndexUL=tmp;
		}
	else if(dir==3){//the ancestor and right des have been swapped
		claIndexDown=claIndexUR;
		claIndexUR=tmp;		
		}
	else assert(0);
	}	

void TreeNode::getTaxonSwapList(vector<int> &TaxonSwapList){
  	if(left == NULL){
		TaxonSwapList.push_back(nodeNum); 
		}
	else {
	  	left->getTaxonSwapList(TaxonSwapList);
	    right->getTaxonSwapList(TaxonSwapList);
		}
	}

void TreeNode::getNNIList(vector<int> &NNIList){
	if(left!=NULL){
		NNIList.push_back(nodeNum);
    	left->getNNIList(NNIList);
    	right->getNNIList(NNIList);
	    }
	else
		return;
	}

void TreeNode::getSPRList(int cut, vector<int> &SPRList){ 
	if(nodeNum==cut) return;

	if(left){
		//if this is an internal node, only add it if it's not cut's ancestor
		if((right->nodeNum!=cut) && (left->nodeNum!=cut)){
			SPRList.push_back(nodeNum);
			}
		left->getSPRList(cut,SPRList);
		right->getSPRList(cut,SPRList);
		}
	else{//if its a terminal node
		SPRList.push_back(nodeNum);		
		}
	}
