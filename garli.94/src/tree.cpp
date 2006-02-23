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

#include <algorithm>
#include <vector>
#include <list>
#ifdef UNIX
	#include <sys/mman.h>
#endif

using namespace std;

#include "mlhky.h"
#include "parameters.h"
#include "clamanager.h"
#include "funcs.h"
#include "stopwatch.h"
#include "model.h"
#include "defs.h"
#include "tree.h"
#include "datamatr.h"

#include "memchk.h"

extern rng rnd;

//external global variables
extern int calcCount;
extern int optCalcs;
extern ofstream opt;
extern ofstream optsum;
extern int memLevel;

//Tree static definitions
double Tree::meanBrlenMuts;
double Tree::alpha;
double Tree::min_brlen	   = DEF_MIN_BRLEN;  // branch lengths never below this value
double Tree::max_brlen	   = DEF_MAX_BRLEN;
double Tree::exp_starting_brlen    = DEF_STARTING_BRLEN; // expected starting branch length
ClaManager *Tree::claMan;
list<TreeNode *> Tree::nodeOptVector;
HKYData *Tree::data;
int Tree::rescaleEvery;
double Tree::treeRejectionThreshold;

void InferStatesFromCla(char *states, double *cla, int nchar);
double CalculateHammingDistance(const char *str1, const char *str2, int nchar);
void SampleBranchLengthCurve(double (*func)(TreeNode*, Tree*, double, bool), TreeNode *thisnode, Tree *thistree);
double CalculatePDistance(const char *str1, const char *str2, int nchar);
inline double CallBranchLike(TreeNode *thisnode, Tree *thistree, double blen, bool brak);
bool Rescale(CondLikeArray *destCLA, int nsites);
bool RescaleRateHet(CondLikeArray *destCLA, int nsites);

#ifdef GANESH //some stuff added by Ganesh for pecr move
/* the value of p in p_value */
int    Tree::p_value = 3; 
int   *Tree::C = 0;
int   **Tree::realcat_intervals = 0;
/* inv_realcat_intervals[0...C[3]-1] */
int   **Tree::inv_realcat_intervals = 0;
bool   Tree::random_p = false;

void Tree::ComputeRealCatalan() {

    int p = p_value;

    C = new int[p+1];

    C[0] = 1;
    C[1] = 2;

    for (int i = 2; i <= p; i++) {
        C[i] = C[i-1];
        for (int j = 1; j <= i-1; j++) {
            C[i] = C[i] + C[j-1]*C[i-1-j];
        }
        C[i] = C[i] + C[i-1];
    }

    realcat_intervals = new int*[p+1];
    inv_realcat_intervals = new int*[p+1];

    for (int a = p; a >= 0; a--) {
        realcat_intervals[a] = new int[a+1]; 
        inv_realcat_intervals[a] = new int[C[a]];

        realcat_intervals[a][0] = C[a-1];

        for (int j = 1; j <= a-1; j++) {
            realcat_intervals[a][j] = realcat_intervals[a][j-1] + (C[j-1] * C[a-j-1]);
        }

        realcat_intervals[a][a] = realcat_intervals[a][a-1] + C[a-1];

        int index = 0;
        for (int j = 0; j < C[a]; j++) {
#if 1
            if (j < realcat_intervals[a][index]) {
                inv_realcat_intervals[a][j] = index;
            }
            else {
                index++;
                inv_realcat_intervals[a][j] = index;
            }
#else
            for (int b = 0; b <= a; b++) {
                if (j < realcat_intervals[a][b]) {
                    inv_realcat_intervalsa[a][j] = b;
                    break;
                }
            }
#endif
            assert (realcat_intervals[a][inv_realcat_intervals[a][j]] > j);
        }
        assert (index == a);
    }    
    return;
}
#endif

int Tree::BrlenMutate(){
	//random_binomial is now called with the mean number of blen muts, which is easiser to specify across datasets
	//than is a per branch probability
	int numBrlenMuts;
	do{
		numBrlenMuts=rnd.random_binomial(numNodesTotal-1, meanBrlenMuts);
		}while(numBrlenMuts==0);
	for(int i=0;i<numBrlenMuts;i++){
		int branch=GetRandomNonRootNode();
		allNodes[branch]->dlen*=rnd.gamma( Tree::alpha );
		SweepDirtynessOverTree(allNodes[branch]);
		allNodes[branch]->dlen = (allNodes[branch]->dlen > min_brlen ? (allNodes[branch]->dlen < max_brlen ? allNodes[branch]->dlen : max_brlen) : min_brlen);
		}	
	return numBrlenMuts;
	}

void Tree::PerturbAllBranches(){
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		allNodes[i]->dlen*=rnd.gamma(100);
		}
	MakeAllNodesDirty();
	}

void Tree::ScaleWholeTree(double factor/*=-1.0*/){
	if(factor==-1.0) factor = rnd.gamma( Tree::alpha );
	for(int i=numTipsTotal;i<numNodesTotal;i++){
		allNodes[i]->dlen*=factor;
		allNodes[i]->dlen = (allNodes[i]->dlen > min_brlen ? (allNodes[i]->dlen < max_brlen ? allNodes[i]->dlen : max_brlen) : min_brlen);
		assert(!(allNodes[i]->dlen < min_brlen));
		}
	MakeAllNodesDirty();
	lnL=-1;
	}


int Tree::BrlenMutateSubset(vector<int> const &subtreeMemberNodes){
	int numBrlenMuts;
	do{
		numBrlenMuts=rnd.random_binomial(subtreeMemberNodes.size(), meanBrlenMuts);
		}while(numBrlenMuts==0);
	for(int i=0;i<numBrlenMuts;i++){
		int branch=subtreeMemberNodes[(int)(rnd.uniform()*subtreeMemberNodes.size())];//can't mutate the root
		allNodes[branch]->dlen*=rnd.gamma( Tree::alpha );
		SweepDirtynessOverTree(allNodes[branch]);
		allNodes[branch]->dlen = (allNodes[branch]->dlen > min_brlen ? (allNodes[branch]->dlen < max_brlen ? allNodes[branch]->dlen : max_brlen) : min_brlen);
		}	
	return numBrlenMuts;
	}

//DZ 10-31-02
//separating general tree construction stuff from CLA assignment/allocation
//which will happend differently if the CLAs are shared or not
Tree::Tree(){
	allNodes=new TreeNode*[2*data->NTax()-2];
	for(int i=0;i<2*data->NTax()-2;i++){
		allNodes[i]=new TreeNode(i);
		allNodes[i]->bipart=new Bipartition();
		}
	root=allNodes[0];
	root->attached=true;
	//assign data to tips
	for(int i=1;i<=data->NTax();i++){
		allNodes[i]->tipData=data->GetAmbigString(i-1);
		}
	
	numTipsAdded=0;
	numNodesAdded=1;//root
	numTipsTotal=data->NTax();
	numNodesTotal=2*data->NTax()-2;
	lnL=0.0;
	calcs=0;
	numBranchesAdded=0;
	taxtags=new int[numTipsTotal+1];
	}

void Tree::AllocateTree(){
	calcs=0;
	allNodes=new TreeNode*[2*data->NTax()-1];
	for(int i=0;i<2*data->NTax()-1;i++){
		allNodes[i]=new TreeNode(i);
		allNodes[i]->bipart=new Bipartition();
		}
	root=allNodes[0];
	//assign data to tips
	for(int i=1;i<=data->NTax();i++){
		allNodes[i]->tipData=data->GetAmbigString(i-1);
		}
	
	numTipsAdded=0;
	numNodesAdded=1;//root
	numTipsTotal=data->NTax();
	numNodesTotal=2*data->NTax()-1;
	lnL=0.0;
	numBranchesAdded=0;
	taxtags=new int[numTipsTotal+1];
	}

//DJZ 4-28-04
//adding the ability to read in treestrings in which the internal node numbers are specified.  I'd like to make the
//internal numbers be specified the way that internal node labels are according to the newick format, ie directly after
//the closing paren that represents the internal node.  But, that makes going from string -> tree annoying
//because by the time the internal node number would be read the treeNode structure would have already been created.
//So, the internal node numbers will go just BEFORE the opening paren that represents that node
//Example:  50(1:.05, 2:.02):.1 signifies a node numbered 50 that is ancestral to 1 and 2.
Tree::Tree(const char* s, bool numericalTaxa /*=true*/){
	//if we are using this constructor, we can't guarantee that the tree will be specified unrooted (with 
	//a trifurcating root), so use an allocation function that is guaranteed to have enough room and then
	//trifurcate and delete if necessary
	
	AllocateTree();
	TreeNode *temp=root;
	root->attached=true;
	int current=numTipsTotal+1;
	bool cont=false;
//	numBranchesAdded=-1;//if reading in a tree start at -1 so opening ( doesn't add a branch
	while(*s)
		{
		cont = false;

		if(*s == ';')
			break;  // ignore semicolons
		else if(*s == ' ')
			break;  // ignore spaces
		else if(*s == ')'){
			assert(temp->anc);
			temp=temp->anc;
			s++;
			while(*s && !isgraph(*s))
				s++;
			if(*s==':')
				{NxsString info = "";
				while( *(s+1)!=')'&& *(s+1)!=','){
						info+=*(s+1);
						s++;
					}
				temp->dlen = atof( info.c_str() );
				temp->dlen = (temp->dlen > min_brlen ? temp->dlen : min_brlen);
				s++;
				}
			else if(*s==','||*s==')')
				temp->dlen=.05;
			else
				{
				if(*s==';'){
					s++;
					}
				assert(!*s  || *s==';');
				}
			}
		else if(*s == ','){
			assert(temp->anc);
			temp=temp->anc;
			numBranchesAdded--;//sloppy way to avoid over incrementing numBranchesAdded
			if(*(s+1)!='(') s++;
			cont = true;
			}
		if(*s == '(' || isdigit(*s) || cont==true){
			//here we're about to add a node of some sort
			numBranchesAdded++;
			if(*(s+1)=='('){//add an internal node
				temp=temp->AddDes(allNodes[current++]);
				numNodesAdded++;
				s++;
				}
			else{
				//this gets ugly.  At this point we could be adding an internal node with the internal node
				//num specifed, or a terminal node.  Either way the next characters in the string will be
				//digits.  We'll have to look ahead to see what the next non-digit character is.  If it's
				//a '(', we know we are adding a prenumbered internal
				if(*s=='(') s++;
				int i=0;
				bool term=true;
				while(isdigit(*(s+i))) i++;
				if(*(s+i) == '(') term=false;
				
				if(term == false){//add an internal node with the nodenum specified in the string
					NxsString num;
					num = *s;
					while(isdigit(*(s+1))){
						assert(*s);
						num += *++s;
						}
					int internalnodeNum = atoi( num.c_str() );
	                temp=temp->AddDes(allNodes[internalnodeNum]);
	                numNodesAdded++;
	                s++;							
					}
				else{//add a terminal node
					// read taxon name
	                 NxsString name;
	                 name = *s;
	                 int taxonnodeNum;
					if(numericalTaxa==true){
						while(isdigit(*(s+1))){
							assert(*s);
							name += *++s;
							}						
						taxonnodeNum = atoi( name.c_str() );
						}
					else{
						while(*(s+1) != ':' && *(s+1) != ',' && *(s+1) != ')'){
							assert(*s);
							name += *++s;
							}
						taxonnodeNum = data->TaxonNameToNumber(name);
						if(taxonnodeNum < 0) throw ErrorException("Unknown taxon \"%s\" encountered in tree description!", name.c_str());
						}
	                temp=temp->AddDes(allNodes[taxonnodeNum]);
	                numNodesAdded++;
	                s++;
			
					if(*s!=':' && *s!=',' && *s!=')'){
						s--;	
						ofstream str("treestring.log", ios::app);
						str << s << endl;
						str.close();
						assert(0);
						}
						
	                if(*s==':'){
						NxsString info = "";
						while( *(s+1)!=')'&& *(s+1)!=','){
							info+=*(s+1);
							s++;
							}
						s++;
						temp->dlen = atof( info.c_str() );
						temp->dlen = (temp->dlen > min_brlen ? temp->dlen : min_brlen);
						}
					else temp->dlen= Tree::exp_starting_brlen;// * rnd.gamma( alpha );
					}
				}
			}
		}
	if(root->left->next==root->right){
		MakeTrifurcatingRoot(true, false);	
		}
	else	{
		EliminateNode(2*data->NTax()-2);
		numBranchesAdded--;
		}
	assert(root->left->next!=root->right);

	root->CheckforLeftandRight();
	root->CheckforPolytomies();
	root->CheckTreeFormation();
	} 

Tree::~Tree(){
	if(taxtags!=NULL) delete []taxtags;
	if(allNodes!=NULL){
		for(int x=0; x<numNodesTotal; x++){
			delete *(allNodes+x);
			}
		delete []allNodes;
		}
	}

void Tree::MakeTrifurcatingRoot(bool reducenodes, bool clasAssigned ){
	//reducenodes should only =1 if this function is called after generating a random tree
	//or after reading in a tree with a bifurcating root.  DO NOT call with reducenodes=1 if
	//this is being used after one of the initial root branches was pruned off

	//clasAssigned should be true if the clas have been assigned to the nodes by the claManager.
	//(ie, not right after tree creation)
	TreeNode *t1, *t2, *removedNode;
	double l;
	assert(root->left->next==root->right);
	if(root->left->left!=NULL){
		removedNode=root->left;
		t1=root->left->left;
		t2=root->left->right;
		l=root->left->dlen;
		root->right->dlen+=l;
		root->left->attached=false;
		if(clasAssigned){
			root->left->claIndexDown=claMan->SetDirty(root->left->claIndexDown);
			root->left->claIndexUL=claMan->SetDirty(root->left->claIndexUL);
			root->left->claIndexUR=claMan->SetDirty(root->left->claIndexUR);			
			}
		root->left=t1;
		t1->next=t2;
		t2->prev=t1;
		t2->next=root->right;
		root->right->prev=t2;
		t1->anc=root;
		t2->anc=root;
		}
	else {
		removedNode=root->right;
		t1=root->right->left;
		t2=root->right->right;
		l=root->right->dlen;
		root->left->dlen+=l;
	 	root->right->attached=false;
		if(clasAssigned){
			root->right->claIndexDown=claMan->SetDirty(root->right->claIndexDown);
			root->right->claIndexUL=claMan->SetDirty(root->right->claIndexUL);
			root->right->claIndexUR=claMan->SetDirty(root->right->claIndexUR);
			}			
	 	root->left->next=t1;
	 	t1->prev=root->left;
		t1->next=t2;
		t2->prev=t1;
		t2->next=NULL;
		t1->anc=root;
		t2->anc=root;
	 	root->right=t2;
		}
	if(reducenodes==1){
		//we need to permanently get rid of the node that was removed and decrement the nodeNums of those greater 
		//than it.
		SortAllNodesArray();
		EliminateNode(removedNode->nodeNum);
		numBranchesAdded--;
		}
	}
		

void Tree::AddRandomNode(int nodenum , int &placeInAllNodes){
	
	assert(nodenum>0 && nodenum<=numTipsTotal);  //should be adding a terminal
	TreeNode* nd=allNodes[nodenum];
	nd->dlen = Tree::exp_starting_brlen * rnd.gamma( 100 );
	nd->dlen = (nd->dlen > min_brlen ? nd->dlen : min_brlen);
	
	nd->next=nd->prev=NULL;//in case this node was connected in some other tree
	
	// Pick a random tip node to add new node to
	//DZ 11-5-02 Will this work to make a trifurcating root?
	if(numBranchesAdded<3)
		{root->AddDes(nd);
		}
	else
		{// If we're not adding directly to the root node, then we will need
		// a connector node and make the new terminal its left des
		TreeNode* connector=allNodes[placeInAllNodes++];
		numNodesAdded++;
		connector->dlen = Tree::exp_starting_brlen * rnd.gamma( 100 );
		nd->dlen = (nd->dlen > min_brlen ? nd->dlen : min_brlen);
		connector->left=connector->right=NULL;
		connector->AddDes(nd);
		
		//select a branch to break with the connector
		int k = rnd.random_int( numBranchesAdded );
		k++;
		TreeNode* otherDes = root->FindNode( k );
		
		// replace put connection in the tree where otherDes had been
		otherDes->SubstituteNodeWithRespectToAnc(connector);
		
		//add otherDes back to the tree as the sister to the new tip
		connector->AddDes(otherDes);
		numBranchesAdded++;//numBranchesAdded needs to be incremented twice because a total of two branches have been added
		}
	numBranchesAdded++;
	numNodesAdded++;
	numTipsAdded++;
	}

void Tree::MimicTopologyButNotInternNodeNums(TreeNode *copySource,TreeNode *replicate,int &placeInAllNodes){
	//used in recombine so internal node nodeNums don't have to match
	TreeNode *tempno=copySource->left;
	assert(copySource->left);
	while(tempno)
		{if(tempno->left)
			{//tempno isn't a terminal
			placeInAllNodes=FindUnusedNode(placeInAllNodes);
			allNodes[placeInAllNodes]->dlen=tempno->dlen;
//			allNodes[placeInAllNodes]->CopyOneClaIndex(copySource, claMan);
			MimicTopologyButNotInternNodeNums(tempno,replicate->AddDes(allNodes[placeInAllNodes]),placeInAllNodes);
			}
		else
			{allNodes[tempno->nodeNum]->dlen=tempno->dlen;
			replicate->AddDes(allNodes[tempno->nodeNum]);
			}
		tempno=tempno->next;
		}	
	}
	
void Tree::RecombineWith( Tree *t, bool sameModel, double optPrecision ){
	//note that this function will loop infinately right now if the tree is too small
	//(ie, there are no suitable nodes to choose to recombine with)

	//mark all of the tags as present in this;
	for(int i=1;i<=numTipsTotal;i++)
		taxtags[i]=0;

	// Pick a random internal node that is the source of the subtree that will be copied into both trees
	int k;
	TreeNode* cop;
	bool sfound=false;
	
	while(!sfound){//find a non trivial clade to add to this
		//k = rnd.random_int( t->numBranchesAdded-1);
		//cop = t->root->FindNode( ++k);
		//don't bother picking terminal nodes
		k=t->GetRandomInternalNode();
		cop=t->allNodes[k];
		if(cop->left->left || cop->right->left){ // cop isn't a two node sub tree
			if(cop->anc) //check to make sure there are at least 2 nodes "below"cop on the source tree
					{if(cop->anc->anc)
						sfound=true;
					else
						{if(t->root->left!=cop)
							{if(t->root->left->left)
								sfound=true;
							}
						if(!sfound && t->root->left->next!=cop)
							{if(t->root->left->next->left)
								sfound=true;
							}
						if(!sfound && t->root->right!=cop)
							{if(t->root->right->left)
								sfound=true;
							}
						}
					}
				}
			}
	//Prune terminals off of this to prepare for attachement of a copy of cop
	cop->left->MarkTerminals(taxtags);
	for(int i=1;i<=numTipsTotal;i++){
		if(taxtags[i]){
			//before removing the tip, trace dirtyness from its anc to the root
			//make sure to set any
			TraceDirtynessToRoot(allNodes[i]->anc);
//			TraceDirtynessToRoot(allNodes[i]);
			allNodes[i]->Prune();
			if(root->left->next==root->right) MakeTrifurcatingRoot(false, true);
			}
		}
	int numAttachedToRoot=root->CountBranches(0);
	
	//what we'd like to do now is make the nodeNums of the subtree that will be attached to this
	//the same as they were in the source tree.  This will require swapping some nodes in the allNodes array,
	//but will simplify other things, and allow us not to recalc some clas.  This is a bit dangerous though, as 
	//the nodeNums in this that correspond to those in the cop subtree are now technically free, but are still
	//marked as attached.  There should still be one node in this marked as unattached that will be used for
	//the connector
	SwapAndFreeNodes(cop);
		
	// Pick a random node whose branch we will bisected by the new subtree.
	int n = rnd.random_int( numAttachedToRoot );
	TreeNode* broken = root->FindNode( ++n );
	assert(broken->anc);//broken can't be the root;
	
	//DZ 7-6 rewritting this so that broken keeps it's original dlen and connector has a new one
	//generated.  Exactly how this would be best done is not clear.  For now picking uniform[0.05,0.2]
	TreeNode *connector;
	int nextUnconnectedNode=FindUnusedNode(numTipsTotal+1);
	connector=allNodes[nextUnconnectedNode];
	connector->left=connector->right=NULL;
	broken->SubstituteNodeWithRespectToAnc(connector);
	connector->AddDes(broken);
	connector->AddDes(allNodes[cop->nodeNum]);
	MimicTopo(cop, 1, sameModel);

	//place connector midway along the broken branch
	connector->dlen=broken->dlen*.5;
	broken->dlen-=connector->dlen;

	TraceDirtynessToRoot(connector);
	OptimizeBranchesAroundNode(connector, optPrecision, 0);
	}

TreeNode *Tree::ContainsBipartition(Bipartition *bip){
	//note that this doesn't work for terminals (but there's no reason to call for them anyway)
	//find a taxon that appears "on" in the bipartition
	int tax=bip->FirstPresentTaxon();
	
	//now start moving down the tree from that taxon until a bipart that
	//conflicts or a match is found
	TreeNode *nd=allNodes[tax]->anc;
	while(nd->anc){
		if(nd->bipart->IsIncompatible(bip))	return NULL;
		else if(nd->bipart->EqualsEquals(bip)) return nd;
		else nd=nd->anc;
		}
	return NULL;
	}

TreeNode *Tree::ContainsBipartitionOrComplement(Bipartition *bip){
	//this version will detect if the same bipartition exists in the trees, even
	//if it is in different orientation, which could happen due to rooting
	//differences

	//find a taxon that appears "on" in the bipartition
	int tax=bip->FirstPresentTaxon();
	
	//now start moving down the tree from that taxon until a bipart that
	//conflicts or a match is found
	TreeNode *nd=allNodes[tax]->anc;
	while(nd->anc){
		if(nd->bipart->IsIncompatible(bip)) break;
		else if(nd->bipart->EqualsEquals(bip)) return nd;
		else nd=nd->anc;
		}
		
	//find a taxon that is NOT "on" in the bipartition
	tax=bip->FirstNonPresentTaxon();
	
	//now start moving down the tree from that taxon until a bipart that
	//conflicts or a match is found
	nd=allNodes[tax]->anc;
	while(nd->anc){
		if(nd->bipart->IsIncompatibleComplement(bip)){
			return NULL;
			}
		else if(nd->bipart->EqualsEqualsComplement(bip)) return nd;
		else nd=nd->anc;
		}
		
	return NULL;
	}

int Tree::SubtreeBasedRecombination( Tree *t, int recomNodeNum, bool sameModel, double optPrecision){
	//this will work more or less like the normal bipartition based recombination, except
	//that the node at which the recombination will occur will be passed in from the population
	//which knows what subtree each remote is working on
	
	//we are assuming that the recomNodeNum represents the same bipartition (subtree) in each tree
	
	TreeNode *tonode=allNodes[recomNodeNum];
	TreeNode *fromnode=t->allNodes[recomNodeNum];
		
	tonode->MarkUnattached(true);
	SwapAndFreeNodes(fromnode);
	//manually set up the base of the subtree in the totree and point tonode to it
	TreeNode *tempanc=tonode->anc;
	TreeNode *tempnext=tonode->next;
	TreeNode *tempprev=tonode->prev;
	if(tempanc->left==tonode){
		tempanc->left=allNodes[fromnode->nodeNum];
		tonode=tempanc->left;
		}
	else if(tempanc->right==tonode){
		tempanc->right=allNodes[fromnode->nodeNum];
		tonode=tempanc->right;			
		}
	else{
		tempanc->left->next=allNodes[fromnode->nodeNum];
		tonode=tempanc->left->next;				
		}
	tonode->anc=tempanc;
	tonode->next=tempnext;
	tonode->prev=tempprev;
	if(tempnext) tempnext->prev=tonode;
	if(tempprev) tempprev->next=tonode;
	MimicTopo(fromnode, 1, sameModel);
	if(sameModel==true) CopyClaIndecesInSubtree(fromnode, true);
	else DirtyNodesInSubtree(tonode);
	
	SweepDirtynessOverTree(tonode);
	
	//try branch length optimization of tonode's branch, to make sure it fits in it's new tree background 
	OptimizeBranchLength(optPrecision, tonode, true);
	return 1;
	}


bool Tree::IdenticalSubtreeTopology(const TreeNode *other){
	//This should not be called with the root, and only detects identical subtrees
	//in the same orientation (ie rooting can fool it)
	assert(other->anc != NULL);
	bool identical;
	
	if(other->anc!=NULL){
		if(other->left==NULL) return true;
		identical=ContainsBipartition(other->bipart);
		if(identical==true){
			identical=IdenticalSubtreeTopology(other->left);
			if(identical==true)
				identical=IdenticalSubtreeTopology(other->right);
			}
		}
	
	return identical;
	}

bool Tree::IdenticalTopology(const TreeNode *other){
	//this is intitially called with the root, it will detect any difference in the 
	//overall topology
	bool identical;
	
	if(other->anc!=NULL){
		if(other->left==NULL) return true;
		identical=ContainsBipartitionOrComplement(other->bipart);
		if(identical==true){
			identical=IdenticalTopology(other->left);
			if(identical==true)
				identical=IdenticalTopology(other->right);
			}
		}
	else{
		TreeNode *nd=other->left;
		while(nd != NULL){
			identical=IdenticalTopology(nd);
			if(identical == false){
				return identical;
				}
			nd=nd->next;
			}
		}
	
	return identical;

	}

int Tree::BipartitionBasedRecombination( Tree *t, bool sameModel, double optPrecision){
	//find a bipartition that is shared between the trees
	TreeNode *tonode, *fromnode;
	bool found=false;
	int tries=0;
	while(!found && (++tries<50)){
		int i;
		do{
			i=t->GetRandomInternalNode();
			}while((t->allNodes[i]->left->left==NULL && t->allNodes[i]->right->left==NULL));
		fromnode=t->ContainsBipartition(allNodes[i]->bipart);
		if(fromnode != NULL){
			//OK the biparts match, but see if they share the same clas!!!!
			//Not much point in scoring them then.
			tonode=allNodes[i];
			if(!((tonode->nodeNum == fromnode->nodeNum) && (tonode->claIndexDown == fromnode->claIndexDown))){
				if(IdenticalSubtreeTopology(fromnode->left)==false) found=true;
				if(found==false) if(IdenticalSubtreeTopology(fromnode->right)==false) found=true;
				}
			}
		}
		//sum the two subtrees as if they were the root to see which is better in score
/*		if(found==true){
			double toscore, fromscore;
			toscore=SubTreeScore(tonode);
			fromscore=t->SubTreeScore(fromnode);
		
			if(fromscore > (toscore + .1)){
				found=true;
				break;
				}
			else found=false;
			}
*/	
	if(found==true){
		tonode->MarkUnattached(true);
		SwapAndFreeNodes(fromnode);
		//manually set up the base of the subtree in the totree and point tonode to it
		TreeNode *tempanc=tonode->anc;
		TreeNode *tempnext=tonode->next;
		TreeNode *tempprev=tonode->prev;
		if(tempanc->left==tonode){
			tempanc->left=allNodes[fromnode->nodeNum];
			tonode=tempanc->left;
			}
		else if(tempanc->right==tonode){
			tempanc->right=allNodes[fromnode->nodeNum];
			tonode=tempanc->right;			
			}
		else{
			tempanc->left->next=allNodes[fromnode->nodeNum];
			tonode=tempanc->left->next;				
			}
		tonode->anc=tempanc;
		tonode->next=tempnext;
		tonode->prev=tempprev;
		if(tempnext) tempnext->prev=tonode;
		if(tempprev) tempprev->next=tonode;
		MimicTopo(fromnode, 1, sameModel);
		if(sameModel==true) CopyClaIndecesInSubtree(fromnode, true);
		else DirtyNodesInSubtree(tonode);
		
		//try branch length optimization of tonode's branch, to make sure it fits in it's new tree background 
		SweepDirtynessOverTree(tonode);
		OptimizeBranchLength(optPrecision, tonode, true);
		Score(tonode->nodeNum);
		}
	else return -1;
	return 1;
	}



void Tree::LocalMove(){
	//this will all assume that there are no polytomies besides the root
	TreeNode *a, *b, *c, *d;
	int cPosition;

	//pick a random TreeNode and set up the rest of the nodes in relation to it//
	int tmp=numTipsTotal+rnd.random_int(numTipsTotal-3)+1;
	TreeNode *u=allNodes[tmp];
	//set up the vicinity of u
	TreeNode *v=u->anc;

	//STANDARDIZE by making v->left=u
	if(u!=v->left){
		if(u==v->left->next){
			if(v->anc==NULL){
				TreeNode *tempnode=v->left;
				TreeNode *tempnode2;
				if(v->left->next==u) tempnode2=u->next;
				else tempnode2=v->left->next;
				v->left=u;
				u->next=tempnode;
				tempnode->next=tempnode2;
				tempnode2->next=NULL;
				}
			else{
				v->RotateDescendents();
			/*	TreeNode *tempnode=v->left;
				v->left=u;
				u->next=tempnode;
				tempnode->next=NULL;
			*/	}
			}
		else{
			//v must be the root, and u must be v->left->next->next
			TreeNode *tempnode=v->left;
			TreeNode *tempnode2=v->left->next;
			v->left=u;
			u->next=tempnode;
			u->next->next=tempnode2;
			tempnode2->next=NULL;
			}
		}
	//determine a and b
	if(rnd.uniform()<.5){
		a=u->left;
		b=a->next;
		}
	else{
		b=u->left;
		a=b->next;
		}
	//STANDARDIZE by making u->left=a
	if(u->left!=a){
		u->RotateDescendents();
/*		u->left=a;
		u->left->next=b;
		b->next=NULL;
*/		}
	//set up the vicinity of v
	if(v->anc==NULL){
		//if v is the root
		if(rnd.uniform()<.05){
			c=u->next;
			d=c->next;
			//STANDARDIZE by making c=v->left->next->next
			u->next=d;
			d->next=c;
			c->next=NULL;
			cPosition=2;
			}
		else{
			d=u->next;
			c=d->next;//left->next->next
			if(c==NULL){
				c=c;
				}
			cPosition=2;
			}
		}
	else{
		//if v is not the root...
		if(rnd.uniform()<.5){
			c=u->next;
			cPosition=1;//left->next
			d=v->anc;
			//STANDARDIZE by making d->left=v if cPosition==1
			if(d->left!=v){
				if(d->anc!=NULL){
					d->RotateDescendents();
/*					TreeNode *tempnode=d->left;
					d->left=v;
					v->next=tempnode;
					tempnode->next=NULL;
*/					}
				else{
					TreeNode *tempnode=d->left;
					TreeNode *tempnode2;
					if(d->left->next==v) tempnode2=v->next;
					else tempnode2=d->left->next;
					d->left=v;
					v->next=tempnode;
					tempnode->next=tempnode2;
					tempnode2->next=NULL;				
					}
				}
			}
		else{
			d=u->next;
			c=v->anc;
			cPosition=3;//anc
			//STANDARDIZE by making c->left=v if cPosition==3
			if(c->left!=v){
				if(c->anc==NULL){
					TreeNode *tempnode=c->left;
					TreeNode *tempnode2;
					if(tempnode->next==v) tempnode2=v->next;
					else tempnode2=tempnode->next;
					c->left=v;
					v->next=tempnode;
					tempnode->next=tempnode2;
					tempnode2->next=NULL;
					}
				else{
					c->RotateDescendents();
/*					TreeNode *tempnode=c->left;
					c->left=v;
					v->next=tempnode;
					tempnode->next=NULL;
*/					}
				}
			}
		}
	/*Now that things are set up, we can count on the following being true:
		u->left=a;
		u->left->next=b;
		v->left=a;
		if(v->anc!=NULL){
			v->left->next=c && d->left=v (case 1)
			else c->left=v && v->left->next=d (case 2)
			}
		else{
			v->left->next->next=c && v->left->next=d (case 3)
			}
	*/	



	//Ok, the nodes are defined.
	//Calculate the backbone length and the new length
	double m;
	double changing_blens[3];
	double new_blens[3];
	changing_blens[0]=a->dlen;
	changing_blens[1]=u->dlen;
	if(cPosition==3){
		changing_blens[2]=v->dlen;
		}
	else {
		changing_blens[2]=c->dlen;
		}
	m=changing_blens[0]+changing_blens[1]+changing_blens[2];
	double r=rnd.uniform();



//	double tuning=.25;
//	double tuning=.1;
	double mprime=m*exp(.5*(rnd.uniform()-.5));
	double x, y;
	//choose whether to "detach" u or v.  Don't actually detach anything though
	if(rnd.uniform()<.5){ //detach u
		//calculate x and y
		x=rnd.uniform()*mprime;
		y=(a->dlen+u->dlen) * (mprime/m);
		
		if(x<y){//all cases
			//no topo change
			a->dlen=x;
			u->dlen=y-x;
			if(cPosition==3) v->dlen=mprime-y;
			else c->dlen = mprime-y;
			TraceDirtynessToRoot(a->anc);
//			tree->AdjustCLArrayFlagsBelow(a->anc, curMove);
			}
		else{
			//case 1
			if(cPosition==1){
				u->left=b;
				u->left->next=c;
				c->next=NULL;
				c->anc=u;
				v->left->next=a;
				a->anc=v;
				a->next=NULL;
				a->dlen=y;
				u->dlen=x-y;
				c->dlen=mprime-x;
				TraceDirtynessToRoot(c->anc);
				//tree->AdjustCLArrayFlagsBelow(c->anc, curMove);
				}
			//case 3
			else if(cPosition==2){
				u->left=b;
				u->left->next=c;
				c->next=NULL;
				c->anc=u;
				v->left->next->next=a;
				a->anc=v;
				a->next=NULL;
				a->dlen=y;
				u->dlen=x-y;
				c->dlen=mprime-x;
				TraceDirtynessToRoot(c->anc);
				//tree->AdjustCLArrayFlagsBelow(c->anc, curMove);
				}
			//case 2
			else{//u and v physically swap positions in this case
				v->left=a;
				a->anc=v;
				a->next=d;
				d->next=NULL;
				u->left=v;
				u->next=v->next;
				v->next=b;
				b->next=NULL;
				u->anc=c;
				v->anc=u;
				c->left=u;
				a->dlen=y;
				v->dlen=x-y;
				u->dlen=mprime-x;
				TraceDirtynessToRoot(a->anc);
				//tree->AdjustCLArrayFlagsBelow(a->anc, curMove);
				}
			}
		}



	else{
		//"detach" v
		x=a->dlen*(mprime/m);
		y=rnd.uniform() * mprime;
		if(x<y){
			//no topo change
			a->dlen=x;
			u->dlen=y-x;
			if(cPosition==3) v->dlen=mprime-y;
			else c->dlen=mprime-y;
			TraceDirtynessToRoot(a->anc);
//			tree->AdjustCLArrayFlagsBelow(a->anc, curMove);
			}			
		else{
			//case 1
			if(cPosition==1){
				u->left=b;
				u->left->next=c;
				c->next=NULL;
				c->anc=u;
				v->left->next=a;
				a->anc=v;
				a->next=NULL;
				a->dlen=y;
				u->dlen=x-y;
				c->dlen=mprime-x;
				TraceDirtynessToRoot(c->anc);
		//		tree->AdjustCLArrayFlagsBelow(c->anc, curMove);
				}
			//case 3
			else if(cPosition==2){
				u->left=b;
				u->left->next=c;
				c->next=NULL;
				c->anc=u;
				v->left->next->next=a;
				a->anc=v;
				a->next=NULL;
				a->dlen=y;
				u->dlen=x-y;
				c->dlen=mprime-x;
				TraceDirtynessToRoot(c->anc);
	//			tree->AdjustCLArrayFlagsBelow(c->anc, curMove);
				}
			//case 2
			else{//u and v physically swap positions in this case
				v->left=a;
				a->anc=v;
				a->next=d;
				d->next=NULL;
				u->left=v;
				u->next=v->next;
				v->next=b;
				b->next=NULL;
				u->anc=c;
				v->anc=u;
				c->left=u;
				a->dlen=y;
				v->dlen=x-y;
				u->dlen=mprime-x;
				TraceDirtynessToRoot(a->anc);
//				tree->AdjustCLArrayFlagsBelow(a->anc, curMove);
				}
			}
		}
	}	


void Tree::TaxonSwap(int tax1, int tax2, double optPrecision){
//	assert(0);
	TreeNode swap;
	TreeNode *t1, *t2;
	t1=allNodes[tax1];
	t2=allNodes[tax2];

	SweepDirtynessOverTree(t1);
	SweepDirtynessOverTree(t2);
	
//temporarily store t2's info	
	swap.anc=t2->anc;
	swap.prev=t2->prev;
	swap.next=t2->next;
	swap.dlen=t2->dlen;
//copy t1's info into t2
	t2->anc=t1->anc;
	t2->prev=t1->prev;
	t2->next=t1->next;
	t2->dlen=t1->dlen;	
	
	if(t1->anc->left == t1){
		t1->anc->left=t2;
		t1->anc->right->prev=t2;
		}
	else{
		t1->anc->right=t2;
		t1->anc->left->next=t2;
		}
	
//copy t2's info into t1
	t1->anc=swap.anc;
	t1->prev=swap.prev;
	t1->next=swap.next;
	t1->dlen=swap.dlen;	
	
	if(t1->anc->left == t2){
		t1->anc->left=t1;
		t1->anc->right->prev=t1;
		}
	else{
		t1->anc->right=t1;
		t1->anc->left->next=t1;
		}
	
	//Optionally, do some branchlength optimization
	OptimizeBranchesWithinRadius(t1->anc, optPrecision, 0, t2->anc);
/*	OptimizeBranchesAroundNode(t1->anc, optPrecision, 0);
	OptimizeBranchesAroundNode(t2->anc, optPrecision, 0);*/
	}

void Tree::TaxonSwap(int range, double optPrecision){
	assert(0);
	//swap two terminal nodes, within a radius
	int tax1, tax2;
	
	//choose the first taxon
	tax1=GetRandomTerminalNode();
	while(allNodes[tax1]->anc==root)
		tax1=GetRandomTerminalNode();
	sprRange.setseed(tax1);
	
	if(range == 0 ){
	    tax2 = GetRandomTerminalNode();
		}
	else  {
	    for(int i = 0;i<range;i++){
		      int j = sprRange.total;
		      for(int k=0; k < j; k++)
			if(sprRange.front[k]==i){
			  TreeNode *cur=allNodes[sprRange.element[k]];
			  if(cur->left!=NULL) 
			    sprRange.addelement(cur->left->nodeNum, i+1, sprRange.pathlength[k] + cur->left->dlen);
			  if(cur->right!=NULL) 
			    sprRange.addelement(cur->right->nodeNum, i+1, sprRange.pathlength[k] + cur->right->dlen);
			  if(allNodes[sprRange.element[k]]->anc!=NULL) 
			    sprRange.addelement(cur->anc->nodeNum, i+1, sprRange.pathlength[k] + cur->dlen);
				}// end of loop through element of current subset
	   		}// end of loop to findrange
	    int n = 0;
	    do{
			n = rnd.random_int(sprRange.total);
		    } while(allNodes[sprRange.element[n]]->left!=NULL || allNodes[sprRange.element[n]]->anc->nodeNum == 0 || allNodes[sprRange.element[n]]->next == allNodes[tax1] || allNodes[sprRange.element[n]]->prev == allNodes[tax1]);
		tax2 = sprRange.element[n];
		}// end of selection of broken node
	//mark the necessary nodes as dirty
	TraceDirtynessToRoot(allNodes[tax1]->anc);
	TraceDirtynessToRoot(allNodes[tax2]->anc);
	
	//now swap the nodes
	TreeNode swap;
	TreeNode *t1, *t2;
	t1=allNodes[tax1];
	t2=allNodes[tax2];
//temporarily store t2's info	
	swap.anc=t2->anc;
	swap.prev=t2->prev;
	swap.next=t2->next;
	swap.dlen=t2->dlen;
//copy t1's info into t2
	t2->anc=t1->anc;
	t2->prev=t1->prev;
	t2->next=t1->next;
	t2->dlen=t1->dlen;	
	
	if(t1->anc->left == t1){
		t1->anc->left=t2;
		t1->anc->right->prev=t2;
		}
	else{
		t1->anc->right=t2;
		t1->anc->left->next=t2;
		}
	
//copy c2's info into t1
	t1->anc=swap.anc;
	t1->prev=swap.prev;
	t1->next=swap.next;
	t1->dlen=swap.dlen;	
	
	if(t1->anc->left == t2){
		t1->anc->left=t1;
		t1->anc->right->prev=t1;
		}
	else{
		t1->anc->right=t1;
		t1->anc->left->next=t1;
		}
	
	//Optionally, do some branchlength optimization
	OptimizeBranchesAroundNode(t1->anc, optPrecision, 0);
	OptimizeBranchesAroundNode(t2->anc, optPrecision, 0);	
	}

void Tree::VariableNNIMutate(int node, int branch, double optPrecision, int subtreeNode){
	//this is just a spoof version of SPRMutate that will perform the same mutation
	//several times with different levels of optimiation, but will otherwise 
	//maintain exactly the same program flow because it resets the seed
	
	Individual tempIndiv;
	tempIndiv.treeStruct=new Tree();
	
	Individual sourceIndiv;
	sourceIndiv.treeStruct=this;
	sourceIndiv.mod->CopyModel(this->mod);
		
	int savedSeed;
	
	ofstream out("variable.log", ios::app);
	out.precision(9);
	out << "NNI\t" << lnL << "\t";
	
	
	tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, &sourceIndiv);
	
	double prec[5]={.5, .25, .1, .05, .01};
	for(int i=0;i<5;i++){
		savedSeed = rnd.seed();
		tempIndiv.treeStruct->NNIMutate(node, branch, prec[i], subtreeNode);
		out << tempIndiv.treeStruct->lnL << "\t";
		rnd.set_seed(savedSeed);
		tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, &sourceIndiv, true);
		}
	out << "\n";
	
	tempIndiv.treeStruct->RemoveTreeFromAllClas();
	delete tempIndiv.treeStruct;
	tempIndiv.treeStruct=NULL;
	sourceIndiv.treeStruct=NULL;

	NNIMutate(node, branch, optPrecision, subtreeNode);

	ofstream poo("3branchScores.log", ios::app);
	poo << endl;
	poo.close();
	}

void Tree::NNIMutate(int node, int branch, double optPrecision, int subtreeNode){

	TreeNode* connector=NULL;
	TreeNode* cut=NULL;
	TreeNode* broken=NULL;
	TreeNode* sib=NULL;

	assert(node<numNodesTotal);
	connector = allNodes[node];
	assert(connector->left!=NULL);
	
	if(branch==0){
		cut=connector->left;
		sib=connector->right;
		}
	else{
		cut=connector->right;
		sib=connector->left;
		}

	SweepDirtynessOverTree(cut);

	//cut will be attached to connector's next or prev
	if(connector->next!=NULL) broken=connector->next;
	else{
		if(connector->anc==root){
			//special case if connector's anc is root and connector is the rightmost decendent
			broken=connector->prev->prev;
			}
		else broken=connector->prev;
		}

	//take out connector and substitute cut's sib for it
   	connector->SubstituteNodeWithRespectToAnc(sib);

	//establish correct topology for connector and cut nodes
	connector->left=connector->right=cut;
	connector->next=connector->prev=connector->anc=cut->next=cut->prev=NULL;

	//assign branchlengths such that the previous blen of broken is divided between
	//broken and connector
	//cut will keep its original blen.  Connector's old blen will be added to sib
	sib->dlen+=connector->dlen;

	if(broken->dlen*.5 > DEF_MIN_BRLEN){
		connector->dlen=broken->dlen*.5;
		broken->dlen-=connector->dlen;
		}
	else connector->dlen=broken->dlen=DEF_MIN_BRLEN;

	//put everything in its place
	broken->SubstituteNodeWithRespectToAnc(connector);
	connector->AddDes(broken);
	
	//try some branch length optimization
	SweepDirtynessOverTree(connector, cut);
	MakeNodeDirty(connector);

#ifdef OPT_DEBUG
	opt << "NNI\n";
	optsum << "NNI\n";
#endif

	OptimizeBranchesWithinRadius(connector, optPrecision, subtreeNode);
	}
 
 
int Tree::VariableSPRMutate(int range, double optPrecision){
	//this is just a spoof version of SPRMutate that will perform the same mutation
	//several times with different levels of optimiation, but will otherwise 
	//maintain exactly the same program flow because it resets the seed
	
	Individual tempIndiv;
	tempIndiv.treeStruct=new Tree();
	
	Individual sourceIndiv;
	sourceIndiv.treeStruct=this;
	sourceIndiv.mod->CopyModel(this->mod);
		
	int savedSeed;
	
	ofstream out("variable.log", ios::app);
	out.precision(9);
	out << "SPR" << range << "\t" << lnL << "\t";
	
	
	tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, &sourceIndiv);
	
	double prec[5]={.5, .25, .1, .05, .01};
	for(int i=0;i<5;i++){
		savedSeed = rnd.seed();
		tempIndiv.treeStruct->SPRMutate(range, prec[i]);
		out << tempIndiv.treeStruct->lnL << "\t";
		rnd.set_seed(savedSeed);
		tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, &sourceIndiv, true);
		}
	out << "\n";

	tempIndiv.treeStruct->RemoveTreeFromAllClas();
	delete tempIndiv.treeStruct;
	tempIndiv.treeStruct=NULL;
	sourceIndiv.treeStruct=NULL;

	SPRMutate(range, optPrecision);
	ofstream poo("3branchScores.log", ios::app);
	poo << endl;
	poo.close();
	}
 
int Tree::SPRMutate(int range, double optPrecision){
	// Pick a random node whose branch we will cut.
	//the nodes range from 0 to numNodesTotal-1, but we don't want to choose 0
	TreeNode *cut=allNodes[GetRandomNonRootNode()];
	TreeNode *sib;
	//note that this assignment of the sib can be overridden below if cut is attached to the root
	if(cut->next!=NULL) sib=cut->next;
	else sib=cut->prev;

	TreeNode *connector=NULL;
	assert(cut!=NULL);
	
	//determine who the connector node will be.  It will be cut->anc unless that is the root
	//if cut->anc is the root, connector will be one of cut's siblings, which is freed when
	//the basal trichotomy is reestablished after removing cut.
	if(cut->anc->anc != NULL){
		connector=cut->anc;
		}
	else{
		if(root->left!=cut && root->left->left != NULL) connector = root->left;
		else if(root->left->next!=cut && root->left->next->left != NULL) connector = root->left->next;
		else if(root->right!=cut && root->right->left != NULL) connector = root->right;
		else{//this should be quite rare, and means that the three descendents of the root
			//are cut and two terminals, so no viable swap exists, just try again
			int numterms=0;
			if(root->left->left==NULL) numterms++;
			if(root->left->next->left==NULL) numterms++;
			if(root->right->left==NULL) numterms++;
			assert(numterms==2);
			int cutNodeNum=SPRMutate(range, optPrecision);
			return cutNodeNum;
			}
		}
	
	//all clas below cut will need to be recalced
	SweepDirtynessOverTree(cut);
	TreeNode *replaceForConn;
	if(cut->anc->anc){
		//cut is not connected to the root, so we can steal it's ancestor as the new connector
		//we've identified a connector node that we can recycle, now we need to disconnect it correctly
	   	if(cut==connector->left)
	   		{assert(cut->next==connector->right);
	   		replaceForConn=connector->right;
	   		}
	   	else
	   		{assert(cut==connector->right); 
	   		replaceForConn=connector->left;
	   		}
		replaceForConn->dlen+=connector->dlen;
	   	connector->SubstituteNodeWithRespectToAnc(replaceForConn);
	   	}
	else{//cut is connected to the root so we need to steal a non terminal sib node as the connector
		MakeNodeDirty(root);
		//Disconnect cut from the root
		if(cut==root->left){
			root->left=cut->next;
			cut->next->prev=NULL;
			}
		else if(cut==root->right){
			root->right=cut->prev;
			cut->prev->next=NULL;
			}
		else{
			assert(cut->prev==root->left && cut->next==root->right);//can only have a basal trifucation, or we're in trouble
			cut->prev->next=cut->next;
			cut->next->prev=cut->prev;
			}
		//root is now bifurcation
		//preserve branch length info
		if(root->right==connector){
			root->left->dlen+=	connector->dlen;
			sib=root->left;
			}
		else{
			root->right->dlen+=	connector->dlen;	
			sib=root->right;
			}
			
		//add the connectors two desccendants as descendants of the root
		assert(connector->right==connector->left->next);
		connector->SubstituteNodeWithRespectToAnc(connector->left);
		root->AddDes(connector->right);
		MakeNodeDirty(connector);
		}
		
	//establish correct topology for connector and cut nodes
	cut->anc=connector;
	connector->left=connector->right=cut;
	connector->next=connector->prev=connector->anc=cut->next=cut->prev=NULL;

	if(range!=0){
		range=GatherNodesInRadius(range, sib, NULL);
/*	
	//if we are reconnecting within some particular range of branches, gather all the nodes within that range
		sprRange.empty();
		sprRange.setseed(sib->nodeNum);//the subset will be centered on cut's sib
	    for(int i = 0;i<range;i++){
	      int j = sprRange.total;
			for(int k=0; k < j; k++){
				if(sprRange.front[k]==i){
					TreeNode *cur=allNodes[sprRange.element[k]];
					if(cur->left!=NULL) 
					    sprRange.addelement(cur->left->nodeNum, i+1, sprRange.pathlength[k] + cur->left->dlen);
					if(cur->right!=NULL) 
				    	sprRange.addelement(cur->right->nodeNum, i+1, sprRange.pathlength[k] + cur->right->dlen);
					if(cur->next!=NULL){
					    sprRange.addelement(cur->next->nodeNum, i+1, sprRange.pathlength[k] + cur->next->dlen);
					    if(cur->next->next!=NULL){//if cur is the left descendent of the root
					    	sprRange.addelement(cur->next->next->nodeNum, i+1, sprRange.pathlength[k] + cur->next->next->dlen);
					    	}
					    }
					if(cur->prev!=NULL){
					    sprRange.addelement(cur->prev->nodeNum, i+1, sprRange.pathlength[k] + cur->prev->dlen);
					    if(cur->prev->prev!=NULL){//if cur is the right descendent of the root
					    	sprRange.addelement(cur->prev->prev->nodeNum, i+1, sprRange.pathlength[k] + cur->prev->prev->dlen);
					    	}
					    }		    	
					if(cur->anc!=NULL) 
					    sprRange.addelement(cur->anc->nodeNum, i+1, sprRange.pathlength[k] + cur->anc->dlen);
					else
						sprRange.addelement(cur->left->next->nodeNum, i+1, sprRange.pathlength[k] + cur->left->next->dlen);
					}// end of loop through element of current subset
			    }// end of loop to findrange
			}
	    //remove unwanted nodes from the subset
	    sprRange.elementremove(sib->nodeNum);
	    if(sprRange.total>3) sprRange.removennis();
	    sprRange.compact();		
	    if(sprRange.total==0){//depending on how the tree was oriented and what was cut
		//it's possible (rare) that the range is empty now.  In that case just abort on the
		//range mutation and do a random one
			range=0;
			}
*/		}

	//--------------------------------------------------------------------
	// Pick a random node whose branch we will bisected by the new subtree.
	// Limited SPR within rangen, where range defined as the distance between two nodes.
	TreeNode * broken = NULL;
	do{
		if(range == 0 ){
		  int numAttachedToRoot=root->CountBranches(0);
		  int n = rnd.random_int( numAttachedToRoot );
		  broken = root->FindNode( ++n );
		}
		else  {
		    int n = rnd.random_int(sprRange.total);
		    broken = allNodes[sprRange.element[n]];
		  }
		}while(broken==sib);//if we reattach to the sib, we'll get the same topo back
		

	broken->SubstituteNodeWithRespectToAnc(connector);
	connector->AddDes(broken);

	if(broken->dlen*.5 > DEF_MIN_BRLEN){
		connector->dlen=broken->dlen*.5;
		broken->dlen-=connector->dlen;
		}
	else connector->dlen=broken->dlen=DEF_MIN_BRLEN;

	SweepDirtynessOverTree(connector, cut);
	MakeNodeDirty(connector);

#ifdef OPT_DEBUG
	opt << "SPR, lim " << range << "\n";
	optsum << "SPR, lim " << range << "\n";
#endif	

	OptimizeBranchesWithinRadius(connector, optPrecision, 0, sib->anc);
	return cut->nodeNum;
}

int Tree::GatherNodesInRadius(int range, TreeNode *sib, TreeNode *subtreeNode){
	//if we are reconnecting within some particular range of branches, gather all the nodes within that range
	sprRange.clear();
	sprRange.setseed(sib->nodeNum);//the subset will be centered on cut's sib
    for(int i = 0;i<range;i++){
      int j = sprRange.total;
		for(int k=0; k < j; k++){
			if(sprRange.front[k]==i){
				TreeNode *cur=allNodes[sprRange.element[k]];
				if(cur->left!=NULL) 
				    sprRange.addelement(cur->left->nodeNum, i+1, sprRange.pathlength[k] + cur->left->dlen);
				if(cur->right!=NULL) 
			    	sprRange.addelement(cur->right->nodeNum, i+1, sprRange.pathlength[k] + cur->right->dlen);
				if(cur->next!=NULL){
				    sprRange.addelement(cur->next->nodeNum, i+1, sprRange.pathlength[k] + cur->next->dlen);
				    if(cur->next->next!=NULL){//if cur is the left descendent of the root
				    	sprRange.addelement(cur->next->next->nodeNum, i+1, sprRange.pathlength[k] + cur->next->next->dlen);
				    	}
				    }
				if(cur->prev!=NULL){
				    sprRange.addelement(cur->prev->nodeNum, i+1, sprRange.pathlength[k] + cur->prev->dlen);
				    if(cur->prev->prev!=NULL){//if cur is the right descendent of the root
				    	sprRange.addelement(cur->prev->prev->nodeNum, i+1, sprRange.pathlength[k] + cur->prev->prev->dlen);
				    	}
				    }		    	
				if(cur->anc!=NULL && cur->anc!=subtreeNode) 
				    sprRange.addelement(cur->anc->nodeNum, i+1, sprRange.pathlength[k] + cur->anc->dlen);
				else if(cur->anc==NULL)
					sprRange.addelement(cur->left->next->nodeNum, i+1, sprRange.pathlength[k] + cur->left->next->dlen);
				}// end of loop through element of current subset
		    }// end of loop to findrange
		}
    //remove unwanted nodes from the subset
    sprRange.elementremove(sib->nodeNum);
    if(sprRange.total>4) sprRange.removennis();
    sprRange.compact();		
    if(sprRange.total==0){//depending on how the tree was oriented and what was cut
	//it's possible (rare) that the range is empty now.  In that case just abort on the
	//range mutation and do a random one
		assert(subtreeNode==0);
		return 0;
		}
	return range;
	}


//Here's another version of SPR that has the cut and broken nodenums passed in
void Tree::SPRMutate(int cutnum, int broknum, double optPrecision, int subtreeNode, int range){
	TreeNode* cut = allNodes[cutnum];
	TreeNode* prunePoint=NULL;
	TreeNode *sib;
	//note that this assignment of the sib can be overridden below if cut is attached to the root or the subtreeNode
	if(cut->next!=NULL) sib=cut->next;
	else sib=cut->prev;

	TreeNode *connector=NULL;
	
	//determine who the connector node will be.  It will be cut->anc unless that is the root
	//if cut->anc is the root, connector will be one of cut's siblings, which is freed when
	//the basal trichotomy is reestablished after removing cut.
	if(cut->anc->anc != NULL){
		if(cut->anc->nodeNum != subtreeNode){
			connector=cut->anc;
			}
		else{
			//cut is attached to the subtreeNode, so we will have to use it's sib as the connector
			connector=sib;
			sib=connector->left;
			}
		}
	else{
		if(root->left!=cut && root->left->left != NULL) connector = root->left;
		else if(root->left->next!=cut && root->left->next->left != NULL) connector = root->left->next;
		else if(root->right!=cut && root->right->left != NULL) connector = root->right;
		else{//this should be quite rare, and means that the three descendents of the root
			//are cut and two terminals, so no viable swap exists, just try again
			assert(0);
			}
		}
	
	//all clas below cut will need to be recalced
	SweepDirtynessOverTree(cut);
	TreeNode *replaceForConn;
	if(cut->anc->anc){
		if(cut->anc->nodeNum != subtreeNode){
			//cut is not connected to the root, so we can steal it's ancestor as the new connector
		   	if(cut==connector->left){
		   		assert(cut->next==connector->right);
				replaceForConn=connector->right;
		   		}
		   	else{
		   		assert(cut==connector->right); 
		   		replaceForConn=connector->left;
		   		}
			replaceForConn->dlen+=connector->dlen;
		   	connector->SubstituteNodeWithRespectToAnc(replaceForConn);
		   	}
		else{//cut is attached to the subtreeNode, so we will have to use it's sib as the connector
			//connector's two children become the subtreeNodes new children, and connector's dlen gets added to subtreeNodes
			TreeNode *subnode=allNodes[subtreeNode];
			subnode->dlen += connector->dlen;
			SweepDirtynessOverTree(connector);
			subnode->left=connector->left;
			subnode->right=connector->right;
			connector->left->anc=subnode;
			connector->right->anc=subnode;
			}
		}
	else{//cut is connected to the root so we need to steal a non terminal sib node as the connector
		MakeNodeDirty(root);
		//Disconnect cut from the root
		if(cut==root->left){
			root->left=cut->next;
			cut->next->prev=NULL;
			}
		else if(cut==root->right){
			root->right=cut->prev;
			cut->prev->next=NULL;
			}
		else{
			assert(cut->prev==root->left && cut->next==root->right);//can only have a basal trifucation, or we're in trouble
			cut->prev->next=cut->next;
			cut->next->prev=cut->prev;
			}
		//root is now bifurcation
		//preserve branch length info
		if(root->right==connector){
			root->left->dlen+=	connector->dlen;
			sib=root->left;
			}
		else{
			root->right->dlen+=	connector->dlen;	
			sib=root->right;
			}
			
		//add the connectors two desccendants as descendants of the root
		assert(connector->right==connector->left->next);
		connector->SubstituteNodeWithRespectToAnc(connector->left);
		root->AddDes(connector->right);
		MakeNodeDirty(connector);
		}
		
	//establish correct topology for connector and cut nodes
	cut->anc=connector;
	connector->left=connector->right=cut;
	connector->next=connector->prev=connector->anc=cut->next=cut->prev=NULL;

	TreeNode *broken;
	if(broknum > 0) broken=allNodes[broknum];
	else{
		GatherNodesInRadius(range, sib, allNodes[subtreeNode]);
	  	int n = rnd.random_int(sprRange.total);
		broken = allNodes[sprRange.element[n]];
		}

	broken->SubstituteNodeWithRespectToAnc(connector);
	connector->AddDes(broken);

	if(broken->dlen*.5 > DEF_MIN_BRLEN){
		connector->dlen=broken->dlen*.5;
		broken->dlen-=connector->dlen;
		}
	else connector->dlen=broken->dlen=DEF_MIN_BRLEN;

	SweepDirtynessOverTree(connector, cut);

	OptimizeBranchesWithinRadius(connector, optPrecision, subtreeNode, sib->anc);
#ifdef OPT_DEBUG
	opt << "SPR\n";
#endif	
}

//DJZ 8-11-04  This version is only for the master doing SPRs on nodes that aren't in a subtree when subtree
//mode is on.  Basically the only difference is that if the ancestor of the cut node is the root, we need to 
//choose one of the other nonSubtree nodes to make a connector to avoid screwing up the subtree partitioning
void Tree::SPRMutate(int cutnum, int broknum, double optPrecision, const vector<int> &nonSubNodes)
{	assert( numBranchesAdded > 3 );

	TreeNode* cut = allNodes[cutnum];
	assert(cut!=NULL);
	SweepDirtynessOverTree(cut->anc);
	TreeNode *connector;

	if(cut->anc->anc != NULL){
		connector=cut->anc;
		}
	else{
		bool foundAConn=false;
		connector=cut->prev;
		while(connector && !foundAConn)//try previous sibs
			{if(connector->left && find(nonSubNodes.begin(),nonSubNodes.end(),connector->nodeNum)!=nonSubNodes.end())//not a terminal
				foundAConn=true;
			else
				connector=connector->prev;
			}
		if(!foundAConn)
			{connector=cut->next;//that didn't work try the next sibs
			while(connector && !foundAConn)//try previous sibs
				{if(connector->left && find(nonSubNodes.begin(),nonSubNodes.end(),connector->nodeNum)!=nonSubNodes.end())//not a terminal
					foundAConn=true;
				else
					connector=connector->next;
				}
			}
		if(!foundAConn)
			return;//oops by chance we picked a trivial branch to cut, so it goes (if you want to call SPRMutate again that would make sure the tree always changes topo
		}
	
	SweepDirtynessOverTree(cut);
	TreeNode *replaceForConn;
	if(cut->anc->anc){
		//cut is not connected to the root, so we can steal it's ancestor as the new connector
	   	if(cut==connector->left){
	   		replaceForConn=connector->right;
	   		}
	   	else{
	   		replaceForConn=connector->left;
	   		}
		replaceForConn->dlen+=connector->dlen;
	   	connector->SubstituteNodeWithRespectToAnc(replaceForConn);
	   	}
	else{//cut is connected to the root so we need to steal a non terminal sib node as the connector
		//this makes the root totally dirty
		MakeNodeDirty(root);
		
		//Disconnect cut from the root
		if(cut==root->left){
			root->left=cut->next;
			cut->next->prev=NULL;
			}
		else if(cut==root->right){
			root->right=cut->prev;
			cut->prev->next=NULL;
			}
		else{
			assert(cut->prev==root->left && cut->next==root->right);//can only have a basal trifucation, or we're in trouble
			cut->prev->next=cut->next;
			cut->next->prev=cut->prev;
			}
		//root is now bifurcation
		//preserve branch length info
		if(root->right==connector)
			root->left->dlen+=	connector->dlen;
		else
			root->right->dlen+=	connector->dlen;
		//add the connectors two desccendants as descendants of the root
		assert(connector->right==connector->left->next);
		connector->SubstituteNodeWithRespectToAnc(connector->left);
		root->AddDes(connector->right);
		MakeNodeDirty(connector);
		}
		
	//establish correct topology for connector and cut nodes
	cut->anc=connector;
	connector->left=connector->right=cut;
	connector->next=connector->prev=connector->anc=cut->next=cut->prev=NULL;
		
	TreeNode *broken=allNodes[broknum];

	broken->SubstituteNodeWithRespectToAnc(connector);
	connector->AddDes(broken);

	if(broken->dlen*.5 > DEF_MIN_BRLEN){
		connector->dlen=broken->dlen*.5;
		broken->dlen-=connector->dlen;
		}
	else connector->dlen=broken->dlen=DEF_MIN_BRLEN;

	SweepDirtynessOverTree(connector, cut);
	MakeNodeDirty(connector);
	
#ifdef OPT_DEBUG
	opt << "SPR\n";
#endif
	OptimizeBranchesWithinRadius(connector, optPrecision, 0);
	}

void Tree::MimicTopo(const Tree *source){
//DZ 10-25-02 This should be much easier and faster using the allnodes array rather
//than being recursive.  Notice that even if the allNodes array of source is not
//ordered according to nodeNum, the new tree will be.
	TreeNode **allNs=source->allNodes;
	for(int i=0;i<source->numNodesTotal;i++){
		if(allNs[i]->anc!=NULL)
			allNodes[i]->anc=allNodes[allNs[i]->anc->nodeNum];
		else allNodes[i]->anc=NULL;
		if(allNs[i]->left!=NULL){
			allNodes[i]->left=allNodes[allNs[i]->left->nodeNum];
			allNodes[i]->right=allNodes[allNs[i]->right->nodeNum];
			}
		else{
			allNodes[i]->left=NULL;
			allNodes[i]->right=NULL;
			}
		if(allNs[i]->next!=NULL)
			allNodes[i]->next=allNodes[allNs[i]->next->nodeNum];
		else allNodes[i]->next=NULL;
		if(allNs[i]->prev!=NULL)
			allNodes[i]->prev=allNodes[allNs[i]->prev->nodeNum];
		else allNodes[i]->prev=NULL;
		allNodes[i]->dlen=allNs[i]->dlen;
		allNodes[i]->attached=true;
		}
	numNodesTotal=source->numNodesTotal;
	numNodesAdded=source->numNodesAdded;
	numTipsAdded=source->numTipsAdded;
	numBranchesAdded=source->numBranchesAdded;
	}

//this version is used for just copying a subtree,
//but assumes that the nodenums will match.  Automatically
//copys the cla indeces too
void Tree::MimicTopo(TreeNode *nd, bool firstNode, bool sameModel){
	//firstNode will be true if this is the base of the subtree to be copied.
	//if it is true, the anc, next and prev should not be copied for that node 
	//Above the firstNode, nodes will be assumed to be the same nodenum in both trees.  This 
	//allows replicating nodeNums from a certain subtree up, but not in the rest of the tree
	//The cla info will only be copied if the models are identical for the individuals (sameModel==true)
	//otherwise the replicated nodes will be marked as dirty
	TreeNode *mnd;
	mnd=allNodes[nd->nodeNum];
	mnd->attached=true;
	if(!firstNode){
		//stuff that should not be done for the root of the subtree
		if(nd->anc){
			mnd->anc=allNodes[nd->anc->nodeNum];
			}
		else{
			mnd->anc=NULL;
			}
		if(nd->next){
			mnd->next=allNodes[nd->next->nodeNum];
				MimicTopo(nd->next, false, sameModel);
			}
		else 
			mnd->next=NULL;
		if(nd->prev){
			mnd->prev=allNodes[nd->prev->nodeNum];
			}
		else 
			mnd->prev=NULL;
		}
	//this should apply to all nodes
	if(nd->left){ //if this is not a terminal
		mnd->left=allNodes[nd->left->nodeNum];
		mnd->right=allNodes[nd->right->nodeNum];
		MimicTopo(nd->left, false, sameModel);
		}
	else
		mnd->right=mnd->left=NULL;;

	//the clas are now taken care of back where this was called	
/*	if(nd->left){
		if(sameModel==true)
			mnd->CopyOneClaIndex(nd, claMan, DOWN);
		else mnd->claIndexDown=claMan->SetDirty(mnd->claIndexDown);
		}
*/
	mnd->dlen=nd->dlen;
}

void Tree::CopyClaIndecesInSubtree(const TreeNode *from, bool remove){
	//the bool argument "remove" designates whether the tree currently has cla arrays
	//assigned to it or not (if not, it must have come from the unused tree vector)
	//note that we assume that the node numbers and topologies match within the subtree
	assert(from->anc);

	//do the clas down
	if(remove) claMan->DecrementCla(allNodes[from->nodeNum]->claIndexDown);
	allNodes[from->nodeNum]->claIndexDown=from->claIndexDown;
	if(allNodes[from->nodeNum]->claIndexDown != -1) claMan->IncrementCla(allNodes[from->nodeNum]->claIndexDown);
	
	//do the clas up left
	if(remove) claMan->DecrementCla(allNodes[from->nodeNum]->claIndexUL);
	allNodes[from->nodeNum]->claIndexUL=from->claIndexUL;
	if(allNodes[from->nodeNum]->claIndexUL != -1) claMan->IncrementCla(allNodes[from->nodeNum]->claIndexUL);
	
	//do the clas up right
	if(remove) claMan->DecrementCla(allNodes[from->nodeNum]->claIndexUR);
	allNodes[from->nodeNum]->claIndexUR=from->claIndexUR;
	if(allNodes[from->nodeNum]->claIndexUR != -1) claMan->IncrementCla(allNodes[from->nodeNum]->claIndexUR);
	
	if(from->left->left != NULL) CopyClaIndecesInSubtree(from->left, remove);
	if(from->right->left != NULL) CopyClaIndecesInSubtree(from->right, remove);
	}

void Tree::DirtyNodesInSubtree(TreeNode *nd){
	
	MakeNodeDirty(nd);
	if(nd->left->left != NULL) DirtyNodesInSubtree(nd->left);
	if(nd->right->left != NULL) DirtyNodesInSubtree(nd->right);
	
	}

bool Rescale(CondLikeArray *destCLA, int nsites){
		double *destination=destCLA->arr;
		int *underflow_mult=destCLA->underflow_mult;
		destCLA->rescaleRank=0;
		//check if any clas are getting close to underflow
#ifdef UNIX
		madvise(destination, sizeof(double)*4*nsites, MADV_SEQUENTIAL);
		madvise(underflow_mult, sizeof(int)*nsites, MADV_SEQUENTIAL);
#endif
		bool reduceRescale=false;
		double smallest, largest;
		for(int i=0;i<nsites;i++){
			smallest=largest=*(destination);
			
			if(smallest>destination[1]) smallest = destination[1];
			else if(largest<destination[1]) largest = destination[1];
			if(smallest>destination[2]) smallest = destination[2];
			else if(largest<destination[2]) largest = destination[2];
			if(smallest>destination[3]) smallest = destination[3];
			else if(largest<destination[3]) largest = destination[3];
			if(smallest>destination[4]) smallest = destination[4];
			else if(largest<destination[4]) largest = destination[4];

			if(smallest<1e-100 && log(largest)-log(smallest)<1200){
				if(smallest < 1e-300 && largest < 1e-100){
					throw(true);
					}
/*				ofstream resc("rescale.log", ios::app);
				resc << smallest << "\t" << largest << "\n";
				resc.close();
*/				int incr=((int) -log(largest));
				underflow_mult[i]+=incr;
				double mult=exp((double)incr);
				for(int q=0;q<4;q++){
					destination[q]*=mult;
					assert(destination[q]<1e300);
					}
				}
			destination+=4;
			}
		return reduceRescale;
		}


bool RescaleRateHet(CondLikeArray *destCLA, int nsites){
		double *destination=destCLA->arr;
		int *underflow_mult=destCLA->underflow_mult;
		destCLA->rescaleRank=0;
		//check if any clas are getting close to underflow
#ifdef UNIX
		madvise(destination, sizeof(double)*16*nsites, MADV_SEQUENTIAL);
		madvise(underflow_mult, sizeof(int)*nsites, MADV_SEQUENTIAL);
#endif
		bool reduceRescale=false;
		double small1, large1, small2, large2;
		for(int i=0;i<nsites;i++){
#ifdef GCC295
			small1= min(destination[0] , destination[8]);
			large1= max(destination[0] , destination[8]);
			small2= min(destination[1] , destination[9]);
			large2= max(destination[1] , destination[9]);

			small1 = min(small1 , small2);
			large1= max(large1 , large2);
			small2= min(destination[2] , destination[10]);
			large2= max(destination[2] , destination[10]);
			
			small1 = min(small1 , small2);
			large1= max(large1 , large2);
			
			small2= min(destination[3] , destination[11]);
			large2= max(destination[3] , destination[11]);
			
			small1 = min(small1 , small2);
			large1= max(large1 , large2);
			
			small2= min(destination[4] , destination[12]);
			large2= max(destination[4] , destination[12]);

			small1 = min(small1 , small2);
			large1= max(large1 , large2);

			small2= min(destination[5] , destination[13]);
			large2= max(destination[5] , destination[13]);		
			
			small1 = min(small1 , small2);
			large1= max(large1 , large2);

			small2= min(destination[6] , destination[14]);
			large2= max(destination[6] , destination[14]);
			
			small1 = min(small1 , small2);
			large1= max(large1 , large2);

			small2= min(destination[7] , destination[15]);
			large2= max(destination[7] , destination[15]);
			
			small1 = min(small1 , small2);
			large1= max(large1 , large2);
			
#else
			small1= (destination[0] < destination[8] ? destination[0] : destination[8]);
			large1= (destination[0] > destination[8] ? destination[0] : destination[8]);
			small2= (destination[1] < destination[9] ? destination[1] : destination[9]);
			large2= (destination[1] > destination[9] ? destination[1] : destination[9]);			

			small1 = (small1 < small2 ? small1 : small2);
			large1= (large1 > large2 ? large1 : large2);
			small2= (destination[2] < destination[10] ? destination[2] : destination[10]);
			large2= (destination[2] > destination[10] ? destination[2] : destination[10]);				
			
			small1 = (small1 < small2 ? small1 : small2);
			large1= (large1 > large2 ? large1 : large2);
			
			small2= (destination[3] < destination[11] ? destination[3] : destination[11]);
			large2= (destination[3] > destination[11] ? destination[3] : destination[11]);
			
			small1 = (small1 < small2 ? small1 : small2);
			large1= (large1 > large2 ? large1 : large2);
			
			small2= (destination[4] < destination[12] ? destination[4] : destination[12]);
			large2= (destination[4] > destination[12] ? destination[4] : destination[12]);

			small1 = (small1 < small2 ? small1 : small2);
			large1= (large1 > large2 ? large1 : large2);

			small2= (destination[5] < destination[13] ? destination[5] : destination[13]);
			large2= (destination[5] > destination[13] ? destination[5] : destination[13]);				
			
			small1 = (small1 < small2 ? small1 : small2);
			large1= (large1 > large2 ? large1 : large2);

			small2= (destination[6] < destination[14] ? destination[6] : destination[14]);
			large2= (destination[6] > destination[14] ? destination[6] : destination[14]);
			
			small1 = (small1 < small2 ? small1 : small2);
			large1= (large1 > large2 ? large1 : large2);

			small2= (destination[7] < destination[15] ? destination[7] : destination[15]);
			large2= (destination[7] > destination[15] ? destination[7] : destination[15]);			
			
			small1 = (small1 < small2 ? small1 : small2);
			large1= (large1 > large2 ? large1 : large2);
#endif
			if(large1< 1e-5){
				if(large1 < 1e-150){
					throw(1);
					}
					
				int incr=((int) -log(large1))+5;

				underflow_mult[i]+=incr;
				double mult=exp((double)incr);
				for(int r=0;r<4;r++){
					for(int q=0;q<4;q++){
						destination[r*4 + q]*=mult;
						assert(destination[r*4 +q] == destination[r*4 +q]);
						assert(destination[r*4 +q] < 1e50);
						}
					}
				}
			destination+=16;
			}
		return reduceRescale;
		}

bool Tree::ConditionalLikelihood(int direction, TreeNode* nd){
	
	calcCount++;
//	if(calcCount % 1000 == 0) cout << calcCount << endl;
	
	int nsites=data->NChar();
	CondLikeArray *destCLA=NULL;
/*	if(direction==DOWN) destCLA=GetClaDown(nd, false);
	else if(direction==UPRIGHT) destCLA=GetClaUpRight(nd, false);
	else if(direction==UPLEFT) destCLA=GetClaUpLeft(nd, false);
	else destCLA=claMan->GetTempCla();
*/
	TreeNode* Lchild, *Rchild;
	CondLikeArray *LCLA=NULL, *RCLA=NULL, *partialCLA=NULL;
	double Lprmat[16], Rprmat[16];

	if(direction != ROOT){
		//the only complicated thing here will be to set up the two children depending on the direction
		//get all of the clas, underflow mults and pmat set up here, then the actual calc loops below 
		//won't depend on direction
		if(direction==DOWN){
			Lchild=nd->left;
			Rchild=nd->right;
			
			if(Lchild->left!=NULL)
				LCLA=GetClaDown(Lchild);
			if(Rchild->left!=NULL)
				RCLA=GetClaDown(Rchild);
			
			destCLA=GetClaDown(nd, false);
			mod->CalcPmat(Lchild->dlen, Lprmat, false);
			mod->CalcPmat(Rchild->dlen, Rprmat, false);
			}
		else if(direction==UPRIGHT || direction==UPLEFT){
			if(nd->anc){
				Lchild=nd->anc;
						
				if(nd->anc->left==nd)
					LCLA=GetClaUpLeft(Lchild);

				else if(nd->anc->right==nd)
					LCLA=GetClaUpRight(Lchild);

				else
					//watch out here.  This is the case in which we want the 
					//cla at the root including the left and right, but not the
					//middle.  We will confusingly store this in the root's DOWN
					//cla
					LCLA=GetClaDown(Lchild);

				mod->CalcPmat(nd->dlen, Lprmat, false);
			
				if(direction==UPRIGHT) Rchild=nd->left;
				else Rchild=nd->right;
				}
			else{
				if(direction==UPRIGHT){
					Lchild=nd->left;
					Rchild=nd->left->next;
					}
				else{
					Lchild=nd->left->next;
					Rchild=nd->right;
					}
				if(Lchild->left!=NULL)
					LCLA=GetClaDown(Lchild);

				mod->CalcPmat(Lchild->dlen, Lprmat, false);
				}
			
			if(Rchild->left!=NULL)
				RCLA=GetClaDown(Rchild);
			
			if(direction==UPRIGHT) destCLA=GetClaUpRight(nd, false);
			else if(direction==UPLEFT) destCLA=GetClaUpLeft(nd, false);
			
			mod->CalcPmat(Rchild->dlen, Rprmat, false);
			}
		
		if(LCLA!=NULL && RCLA!=NULL)
			//two internal children
			CalcFullCLAInternalInternal(destCLA, LCLA, RCLA, Lprmat, Rprmat, nsites);

		else if(LCLA==NULL && RCLA==NULL){
			//two terminal children
			char *Ldata=Lchild->tipData;
			char *Rdata=Rchild->tipData;
			CalcFullCLATerminalTerminal(destCLA, Lprmat, Rprmat, Ldata, Rdata, nsites);
			}

		else{
			//one terminal, one internal
			char *Ldata;
			if(LCLA==NULL){
				Ldata=Lchild->tipData;
				CalcFullCLAInternalTerminal(destCLA, RCLA, Rprmat, Lprmat, Ldata, nsites);
				}
			else {
				char *Rdata=Rchild->tipData;
				CalcFullCLAInternalTerminal(destCLA, LCLA, Lprmat, Rprmat, Rdata, nsites);
				}
			}
		}
	
	if(direction==ROOT){
		//at the root we need to include to contributions of 3 branches.  Check if we have a
		//valid CLA that already represents two of these three. If so we can save a bit of
		//computation.  This will mainly be the case during blen optimization, when when we 
		//only change one of the branches again and again.
		double prmat[64];
		TreeNode *child;
		CondLikeArray *childCLA=NULL;
		int *childUnderMult=NULL;
//		destCLA=claMan->GetTempCla();
		
		if(claMan->IsDirty(nd->claIndexUL) == false){
			partialCLA=GetClaUpLeft(nd, false);
			child=nd->left;
			if(child->left!=NULL){
				childCLA=GetClaDown(child, true);
				}
			mod->CalcPmat(child->dlen, prmat, false);
			}
		else if(claMan->IsDirty(nd->claIndexUR) == false){
			partialCLA=GetClaUpRight(nd, false);
			child=nd->right;
			if(child->left!=NULL){
				childCLA=GetClaDown(child, true);
				}
			mod->CalcPmat(child->dlen, prmat, false);
			}
		else{//both of the UP clas must be dirty.  We'll use the down one as the 
			//partial, and calc it now if necessary
			if(claMan->IsDirty(nd->claIndexDown) == true)
				partialCLA=GetClaDown(nd, true);
			else partialCLA=GetClaDown(nd, false);
			if(nd->anc!=NULL){
				child=nd->anc;
				if(child->left==nd){
					childCLA=GetClaUpLeft(child, true);							
					}
				else if(child->right==nd){
					childCLA=GetClaUpRight(child, true);
					}
				else{
					//the node down that we want to get must be the root, and this
					//node must be it's middle des.  Remember that the cla for that 
					//direction is stored as the root DOWN direction
					childCLA=GetClaDown(child);
					}
				mod->CalcPmat(nd->dlen, prmat, false);
				}
			else{
				child=nd->left->next;
				if(child->left!=NULL){
					childCLA=GetClaDown(child, true);
					}
				mod->CalcPmat(child->dlen, prmat, false);
				}
			}	
			
		if(childCLA!=NULL)//if child is internal
			CalcFullCLAPartialInternal(destCLA, childCLA, prmat, partialCLA, nsites);
		
		else{
			char *childData=child->tipData;
			CalcFullCLAPartialTerminal(destCLA, partialCLA, prmat, childData, nsites);
			}
		}

	if(destCLA->rescaleRank >= rescaleEvery && direction != ROOT)
		Rescale(destCLA, nsites);

/*	if(direction==DOWN) claMan->SetDirty(nd->claIndexDown, nd->nodeNum, false);
	else if(direction==UPRIGHT) claMan->SetDirty(nd->claIndexUR, nd->nodeNum, false);
	else if(direction==UPLEFT) claMan->SetDirty(nd->claIndexUL, nd->nodeNum, false);
*/	return true;
	}

int Tree::ConditionalLikelihoodRateHet(int direction, TreeNode* nd, bool fillFinalCLA /*=false*/){
	//note that fillFinalCLA just refers to whether we actually want to calc a CLA
	//representing the contribution of the entire tree vs just calcing the score
	//The only reason I can think of for doing that is to calc internal state probs
	//the fuction will then return a pointer to the CLA

	calcCount++;
	
	int nsites=data->NChar();
	CondLikeArray *destCLA=NULL;

	TreeNode* Lchild, *Rchild;
	CondLikeArray *LCLA=NULL, *RCLA=NULL, *partialCLA=NULL;
	double Lprmat[64], Rprmat[64];

	if(direction != ROOT){
		//the only complicated thing here will be to set up the two children depending on the direction
		//get all of the clas, underflow mults and pmat set up here, then the actual calc loops below 
		//won't depend on direction
		if(direction==DOWN){
			Lchild=nd->left;
			Rchild=nd->right;
			
			if(Lchild->left!=NULL)
				LCLA=GetClaDown(Lchild);
			if(Rchild->left!=NULL)
				RCLA=GetClaDown(Rchild);
			
			mod->CalcPmatRateHet(Lchild->dlen, Lprmat, false);
			mod->CalcPmatRateHet(Rchild->dlen, Rprmat, false);
			}
		else if(direction==UPRIGHT || direction==UPLEFT){
			if(nd->anc){
				Lchild=nd->anc;
						
				if(nd->anc->left==nd)
					LCLA=GetClaUpLeft(Lchild);

				else if(nd->anc->right==nd)
					LCLA=GetClaUpRight(Lchild);

				else//watch out here.  This is the case in which we want the cla at the root including the left
					//and right, but not the middle.  We will confusingly store this in the root's DOWN cla
					LCLA=GetClaDown(Lchild);

				mod->CalcPmatRateHet(nd->dlen, Lprmat, false);
			
				if(direction==UPRIGHT) Rchild=nd->left;
				else Rchild=nd->right;
				}
			else{
				if(direction==UPRIGHT){
					Lchild=nd->left;
					Rchild=nd->left->next;
					}
				else{
					Lchild=nd->left->next;
					Rchild=nd->right;
					}
				if(Lchild->left!=NULL)
					LCLA=GetClaDown(Lchild);

				mod->CalcPmatRateHet(Lchild->dlen, Lprmat, false);
				}
			
			if(Rchild->left!=NULL)
				RCLA=GetClaDown(Rchild);
			
			mod->CalcPmatRateHet(Rchild->dlen, Rprmat, false);
			}


		if(direction==DOWN) destCLA=GetClaDown(nd, false);
		else if(direction==UPRIGHT) destCLA=GetClaUpRight(nd, false);
		else if(direction==UPLEFT) destCLA=GetClaUpLeft(nd, false);
	
		if(LCLA!=NULL && RCLA!=NULL)
			//two internal children
			CalcFullCLAInternalInternalRateHet(destCLA, LCLA, RCLA, Lprmat, Rprmat, nsites);

		else if(LCLA==NULL && RCLA==NULL){
			//two terminal children
			CalcFullCLATerminalTerminalRateHet(destCLA, Lprmat, Rprmat, Lchild->tipData, Rchild->tipData, nsites);
			}

		else{
			//one terminal, one internal
			if(LCLA==NULL){
				CalcFullCLAInternalTerminalRateHet(destCLA, RCLA, Rprmat, Lprmat, Lchild->tipData, nsites);
				}
			else {
				CalcFullCLAInternalTerminalRateHet(destCLA, LCLA, Lprmat, Rprmat, Rchild->tipData, nsites);
				}
			}
		}
	
	if(direction==ROOT){
		//at the root we need to include the contributions of 3 branches.  Check if we have a
		//valid CLA that already represents two of these three. If so we can save a bit of
		//computation.  This will mainly be the case during blen optimization, when when we 
		//only change one of the branches again and again.
		double prmat[64];
		TreeNode *child;
		CondLikeArray *childCLA=NULL;
		int *childUnderMult=NULL;
		
		if(claMan->IsDirty(nd->claIndexUL) == false){
			partialCLA=GetClaUpLeft(nd, false);
			child=nd->left;
			if(child->left!=NULL){
				childCLA=GetClaDown(child, true);
				}
			mod->CalcPmatRateHet(child->dlen, prmat, false);
			}
		else if(claMan->IsDirty(nd->claIndexUR) == false){
			partialCLA=GetClaUpRight(nd, false);
			child=nd->right;
			if(child->left!=NULL){
				childCLA=GetClaDown(child, true);
				}
			mod->CalcPmatRateHet(child->dlen, prmat, false);
			}
		else{//both of the UP clas must be dirty.  We'll use the down one as the 
			//partial, and calc it now if necessary
			if(claMan->IsDirty(nd->claIndexDown) == true)
				partialCLA=GetClaDown(nd, true);
			else partialCLA=GetClaDown(nd, false);
			if(nd->anc!=NULL){
				child=nd->anc;
				if(child->left==nd){
					childCLA=GetClaUpLeft(child, true);							
					}
				else if(child->right==nd){
					childCLA=GetClaUpRight(child, true);
					}
				else{
					//the node down that we want to get must be the root, and this
					//node must be it's middle des.  Remember that the cla for that 
					//direction is stored as the root DOWN direction
					childCLA=GetClaDown(child);
					}
				mod->CalcPmatRateHet(nd->dlen, prmat, false);
				}
			else{
				child=nd->left->next;
				if(child->left!=NULL){
					childCLA=GetClaDown(child, true);
					}
				mod->CalcPmatRateHet(child->dlen, prmat, false);
				}
			}	
		
		if(fillFinalCLA==false){
			if(childCLA!=NULL)//if child is internal
				lnL = GetScorePartialInternalRateHet(partialCLA, childCLA, prmat);
			
			else
				lnL = GetScorePartialTerminalRateHet(partialCLA, prmat, child->tipData);
			}
		
		else{
			//this is only for inferring internal states
			//careful!  This will have to be returned manually!!
			int wholeTreeIndex=claMan->AssignClaHolder();
			claMan->FillHolder(wholeTreeIndex, ROOT);
			claMan->ReserveCla(wholeTreeIndex);
			if(childCLA!=NULL)//if child is internal
				CalcFullCLAPartialInternalRateHet(claMan->GetCla(wholeTreeIndex), childCLA, prmat, partialCLA, nsites);
			
			else
				CalcFullCLAPartialTerminalRateHet(claMan->GetCla(wholeTreeIndex), partialCLA, prmat, child->tipData, nsites);
			return wholeTreeIndex;
			}
		}

	if(direction != ROOT)
		if(destCLA->rescaleRank >= rescaleEvery)
			RescaleRateHet(destCLA, nsites);

	return -1;
	}

int Tree::Score(int rootNodeNum /*=0*/){

	TreeNode *rootNode=allNodes[rootNodeNum];

	bool scoreOK=true;
	do{
		try{
			scoreOK=true;
			if(mod->NRateCats()==1)
				ConditionalLikelihood( ROOT, rootNode);
			else
				ConditionalLikelihoodRateHet( ROOT, rootNode);
			}
			catch(int err){
				assert(err==1);
				scoreOK=false;
				MakeAllNodesDirty();
				rescaleEvery -= 2;
				ofstream resc("rescale.log", ios::app);
				resc << "rescale reduced to " << rescaleEvery << endl;
				resc.close();
				if(rescaleEvery<2) rescaleEvery=2;
				}
		}while(scoreOK==false);
/*
	double *cla=claMan->GetTempCla()->arr;
	int *underflow_mult=claMan->GetTempCla()->underflow_mult;

//	double templnL=lnL;
	lnL=SumSiteLikes(cla, underflow_mult);
//	assert(templnL == lnL);
*/
	return 1;
	}
/*
double Tree::SubTreeScore( TreeNode *nd){
	//calculates the likelihood of the tree above the node passed in
	double lnL = 0.0;
	int nSites = data->NChar();
	int ck;

	if(claMan->IsDirty(nd->claIndexDown)){
		if(mod->NRateCats()==1)
		ConditionalLikelihood( DOWN, nd);
		else
			ConditionalLikelihoodRateHet(DOWN, nd);
		}

	double *cla=claMan->GetCla(nd->claIndexDown)->arr;
	int *underflow_mult=claMan->GetCla(nd->claIndexDown)->underflow_mult;

	// loop over all patterns
	long double Lk;
	double siteL;
	int ufcount=0;
	const int *countit=data->GetCounts();
	if(mod->ProportionInvariant()==0.0){
		if(mod->NRateCats()==1){//no invariants or gamma
		for( int k = 0; k < nSites; k++ ){
				Lk =  mod->Pi(0) * cla[0] + mod->Pi(1) * cla[1] + mod->Pi(2) * cla[2] + mod->Pi(3) * cla[3];
				if(Lk<1e-300){
					printf("Underflow! site %d, multiplier %d\n", k, underflow_mult[k]);
					ufcount++;
					}
				cla+=4;
				siteL = (log( Lk ) - underflow_mult[k]);
				lnL += (  *countit++ *  siteL);
				}
			}
		else{//gamma, no invariants
			for( int k = 0; k < nSites; k++ ){
				Lk =  mod->Pi(0) * cla[0] + mod->Pi(1) * cla[1] + mod->Pi(2) * cla[2] + mod->Pi(3) * cla[3];
				Lk +=  mod->Pi(0) * cla[4] + mod->Pi(1) * cla[5] + mod->Pi(2) * cla[6] + mod->Pi(3) * cla[7];
				Lk +=  mod->Pi(0) * cla[8] + mod->Pi(1) * cla[9] + mod->Pi(2) * cla[10] + mod->Pi(3) * cla[11];
				Lk +=  mod->Pi(0) * cla[12] + mod->Pi(1) * cla[13] + mod->Pi(2) * cla[14] + mod->Pi(3) * cla[15];
				if(Lk<1e-300){
					printf("Underflow! site %d, multiplier %d\n", k, underflow_mult[k]);
					ufcount++;
					}
				cla+=16;
				//this is hard coded for 4 equal sized rate cats
				siteL = (log( Lk*.25 ) - underflow_mult[k]);
				lnL += (  *countit * siteL);
				countit++;
				}
			}		
		}
	else {
		double prI=mod->ProportionInvariant();
		int lastConst=data->LastConstant();
		const int *conBases=data->GetConstBases();
		
		if(mod->NRateCats()==1){//invariants without gamma
	for( int k = 0; k < nSites; k++ ){
		assert(0);
		//this isn't valid :mod->Pi(conBases[k]), because the con bases are coded as 1 2 4 8 for amiguity
		Lk =  mod->Pi(0) * cla[0] + mod->Pi(1) * cla[1] + mod->Pi(2) * cla[2] + mod->Pi(3) * cla[3];
		if(Lk<1e-300){
			printf("Underflow! site %d, multiplier %d\n", k, underflow_mult[k]);
			ufcount++;
			}
		cla+=4;
				if(k > lastConst){
					siteL = log( Lk * (1.0-prI)) - underflow_mult[k];
					lnL += ( *countit++ * siteL);
					}
				else{
					siteL = log( Lk * (1.0-prI) + (prI * mod->Pi(conBases[k])) * exp((double)underflow_mult[k]));
					lnL += ( *countit++ * (siteL + underflow_mult[k]));
					}
				}
			}
		else{//gamma and invariants
			double scaledGammaProp=0.25 * (1.0-prI);
			assert(0);
			//this isn't valid :mod->Pi(conBases[k]), because the con bases are coded as 1 2 4 8 for amiguity
			for( int k = 0; k < nSites; k++ ){
				Lk =  mod->Pi(0) * cla[0] + mod->Pi(1) * cla[1] + mod->Pi(2) * cla[2] + mod->Pi(3) * cla[3];
				Lk +=  mod->Pi(0) * cla[4] + mod->Pi(1) * cla[5] + mod->Pi(2) * cla[6] + mod->Pi(3) * cla[7];
				Lk +=  mod->Pi(0) * cla[8] + mod->Pi(1) * cla[9] + mod->Pi(2) * cla[10] + mod->Pi(3) * cla[11];
				Lk +=  mod->Pi(0) * cla[12] + mod->Pi(1) * cla[13] + mod->Pi(2) * cla[14] + mod->Pi(3) * cla[15];
				if(Lk<1e-300){
					printf("Underflow! site %d, multiplier %d\n", k, underflow_mult[k]);
					ufcount++;
					}
				cla+=16;
				if(k > lastConst){
					siteL = log( Lk * scaledGammaProp) - underflow_mult[k];
					lnL += ( *countit++ * siteL);
					}
				else{
					siteL = log( Lk * scaledGammaProp + (prI * mod->Pi(conBases[k])) * exp((double)underflow_mult[k]));
					lnL += ( *countit++ * (siteL + underflow_mult[k]));
					}
				}
			}
		}
	return lnL;
	}
*/
/*
double Tree::SubTreeScoreRateHet( TreeNode *nd){
	//calculates the likelihood of the tree above the node passed in
	double sublnL = 0.0;
	int nSites = data->NChar();
	int ck;

	if(claMan->IsDirty(nd->claIndexDown))
		ConditionalLikelihoodRateHet(DOWN, nd);


	double *cla=claMan->GetCla(nd->claIndexDown)->arr;
	int *underflow_mult=claMan->GetCla(nd->claIndexDown)->underflow_mult;

	// loop over all patterns
	long double Lk;
	int ufcount=0;
	const int *countit=data->GetCounts();
	for( int k = 0; k < nSites; k++ ){
		Lk =  mod->Pi(0) * cla[0] + mod->Pi(1) * cla[1] + mod->Pi(2) * cla[2] + mod->Pi(3) * cla[3];
		Lk +=  mod->Pi(0) * cla[4] + mod->Pi(1) * cla[5] + mod->Pi(2) * cla[6] + mod->Pi(3) * cla[7];
		Lk +=  mod->Pi(0) * cla[8] + mod->Pi(1) * cla[9] + mod->Pi(2) * cla[10] + mod->Pi(3) * cla[11];
		Lk +=  mod->Pi(0) * cla[12] + mod->Pi(1) * cla[13] + mod->Pi(2) * cla[14] + mod->Pi(3) * cla[15];
		if(Lk<1e-300){
			printf("Underflow! site %d, multiplier %d\n", k, underflow_mult[k]);
			ufcount++;
			}
		cla+=16;

		sublnL += (  *countit * (log( Lk*.25 ) - underflow_mult[k]) );
		countit++;
		}
	return sublnL;
}
*/
void Tree::TraceDirtynessToRoot(TreeNode *nd){
	SweepDirtynessOverTree(nd);

/*	
	while(nd){
		if(nd->nodeNum==0 || nd->nodeNum>numTipsTotal) nd->claIndexDown=claMan->SetDirty(nd->claIndexDown, true);
		nd=nd->anc;
		}
	
*/	}

void Tree::SweepDirtynessOverTree(TreeNode *nd, TreeNode *from/*=NULL*/){
	lnL=-1;

	if(from==NULL){
		//if this is the branch where the dirtyness starts
		if(nd->left!=NULL){
			nd->claIndexUL=claMan->SetDirty(nd->claIndexUL);
			nd->claIndexUR=claMan->SetDirty(nd->claIndexUR);
			if(nd->left->left!=NULL) SweepDirtynessOverTree(nd->left, nd);
			if(nd->right->left!=NULL) SweepDirtynessOverTree(nd->right, nd);
			}
		if(nd->anc!=NULL) SweepDirtynessOverTree(nd->anc, nd);	
		}
	else{
	//if the change was below, invalidating clas above, also if the change
	//was on the path connecting to the central des of the root
		if(from==nd->anc || (nd->anc==NULL && from==nd->left->next)){
			nd->claIndexUL=claMan->SetDirty(nd->claIndexUL);
			nd->claIndexUR=claMan->SetDirty(nd->claIndexUR);
			if(nd->left->left!=NULL) SweepDirtynessOverTree(nd->left, nd);
			if(nd->right->left!=NULL) SweepDirtynessOverTree(nd->right, nd);
			}
		else if(from==nd->left){
			nd->claIndexUR=claMan->SetDirty(nd->claIndexUR);
			nd->claIndexDown=claMan->SetDirty(nd->claIndexDown);
			if(nd->right->left!=NULL) SweepDirtynessOverTree(nd->right, nd);
			if(nd->anc!=NULL) SweepDirtynessOverTree(nd->anc, nd);		
			else if(nd->left->next->left!=NULL) SweepDirtynessOverTree(nd->left->next, nd);
			}
		else if(from==nd->right){
			nd->claIndexUL=claMan->SetDirty(nd->claIndexUL);
			nd->claIndexDown=claMan->SetDirty(nd->claIndexDown);
			if(nd->left->left!=NULL) SweepDirtynessOverTree(nd->left, nd);
			if(nd->anc!=NULL) SweepDirtynessOverTree(nd->anc, nd);		
			else if(nd->left->next->left!=NULL) SweepDirtynessOverTree(nd->left->next, nd);
			}
		}
	}

void Tree::TraceDirtynessToNode(TreeNode *nd, int tonode){
	if(nd->nodeNum==0 || nd->nodeNum>numTipsTotal) nd->claIndexDown=claMan->SetDirty(nd->claIndexDown);
	while(nd->nodeNum!=tonode){
		nd=nd->anc;
		if(nd->nodeNum==0 || nd->nodeNum>numTipsTotal) nd->claIndexDown=claMan->SetDirty(nd->claIndexDown);
		}
	}

void Tree::SortAllNodesArray(){
	//this function will simply sort the nodes in the allNodes **TreeNode array by their nodeNum
	//having the nodes always in order will make some other operations much simpler
	//the root(nodenum=0) and terminals(nodenums=1->Ntax) should already be in order, so just sort
	//starting at Ntax+1.  I'm making up a kind of wacky algorithm for this. DZ 10-30-02
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		if(allNodes[i]->nodeNum!=i){
			while(allNodes[i]->nodeNum!=i){
				TreeNode *toPlace=allNodes[i];
				int rightPlace=toPlace->nodeNum;
				TreeNode *temp=allNodes[rightPlace];//copy the node that is in toPlace's rightful place
				allNodes[rightPlace]=allNodes[i];   //put toPlace where it belongs
				allNodes[i]=temp;					//put the node that was moved in allNodes[i];
				}
			}
		}
	}

void Tree::EliminateNode(int nn){
	//DZ 10-30-02 this function will permenantly get rid of a node and correct all of the other nodeNums so that 
	//there isn't a hole in the middle.  I think this just needs to be called when an inital tree is made trifrcating
	//at the root.  This isn't the prettiest thing, but I can't think of an obvious way to make a tree that has 3 des
	//from the root in the first place
	delete allNodes[nn];
	for(int i=nn;i<numNodesTotal-1;i++){
		allNodes[i]=allNodes[i+1];
		allNodes[i]->nodeNum=i;
		}
	allNodes[numNodesTotal-1]=NULL;
	numNodesTotal--;
	numNodesAdded--;
	//now make a new allNodes array of the proper length
	TreeNode **newNodes=new TreeNode*[numNodesTotal];
	memcpy(newNodes, allNodes, sizeof(TreeNode*)*numNodesTotal);
	delete []allNodes;
	allNodes=newNodes;
	}
/*
void Tree::CopyCLAs(Tree *source, int nchar){
	//don't copy the root cla, since it is always being recaled if the tree is dirty
	for(int i=1;i<numTipsTotal-2;i++){
		memcpy(allCLAs[i]->arr, source->allCLAs[i]->arr, sizeof(double)*nchar*4.0);
		assert(!(source->allCLAs[i]->IsDirty()));
//		allCLAs[i]->SetDirty(source->allCLAs[i]->IsDirty());
		}
	}
*/	

void Tree::RotateNodesAtRoot(TreeNode *newroot){
	//DZ 11-3-02 This can be used to rebalance the tree
	//I'm assuming that this will be called with one of the des of the root;
	assert(newroot->anc==root);
	assert(newroot->left!=NULL);
	//detach the newroot from root, making it bifurcating
	if(newroot==root->left){
		root->left=newroot->next;
		root->left->prev=NULL;
		}
	else if(newroot==root->left->next){
		root->left->next=root->right;
		root->right->prev=root->left;
		}
	else{
		root->right=root->left->next;
		root->right->next=NULL;
		}
	//now make the root the middle des of newroot and correct the dlens
	root->anc=newroot;
	newroot->left->next=root;
	root->prev=newroot->left;
	root->next=newroot->right;
	newroot->right->prev=root;
	root->dlen=newroot->dlen;
	newroot->dlen=-1;
	newroot->anc=NULL;
	newroot->next=newroot->prev=NULL;
	 //now make the new root nodeNum 0 in the allNodes array
	TreeNode *tempnode=root;
	int tempindexdown=root->claIndexDown;
	root->claIndexDown=newroot->claIndexDown;
	newroot->claIndexDown=tempindexdown;
	int tempindexUL=root->claIndexUL;
	root->claIndexUL=newroot->claIndexUL;
	newroot->claIndexUL=tempindexUL;
	int tempindexUR=root->claIndexUR;
	root->claIndexUR=newroot->claIndexUR;
	newroot->claIndexUR=tempindexUR;
	root=newroot;
	allNodes[0]=newroot;
	tempnode->nodeNum=root->nodeNum;
	root->nodeNum=0;
	allNodes[tempnode->nodeNum]=tempnode;
	//this form of setdirty won't shift every copy to a new topo, but will set them to dirty
//	claMan->SetDirtyButDoNotMove(0, root->claIndex);
//	claMan->SetDirtyButDoNotMove(tempnode->nodeNum, tempnode->claIndex);
	root->claIndexDown=claMan->SetDirty(root->claIndexDown);
	tempnode->claIndexDown=claMan->SetDirty(tempnode->claIndexDown);
	}
	
void Tree::CheckBalance(){
	//DZ 11-=3=02
	//this function will keep the 3 subtrees decending from the root approximately the same
	//depth.  This insures that on average an optimally small number of clas will have to
	//be rescored due to blen or spr mutations.
	
/*	int l=0, m=0, r=0;
	root->left->CalcDepth(l);
	root->left->next->CalcDepth(m);
	root->right->CalcDepth(r);
	
	do{
		//its not exactly clear what a reasonable criterion for rebalancing the tree should be.
		//this seems to work pretty well
		if(l>m+1&&l>r+1) RotateNodesAtRoot(root->left);
		else if(m>l+1&&m>r+1) RotateNodesAtRoot(root->left->next);
		else if(r>m+1&&r>l+1) RotateNodesAtRoot(root->right);
		else return;
		
		l=m=r=0;
		root->left->CalcDepth(l);
		root->left->next->CalcDepth(m);
		root->right->CalcDepth(r);
		}while(1);
*/			
	
	//a more complicated option.  evaluate the average depth of all branches in the tree
	int lb=0, mb=0, rb=0;
	int ls=0, ms=0, rs=0;
	int llb=0, lrb=0, mlb=0, mrb=0, rlb=0, rrb=0;
	int lls=0, lrs=0, mls=0, mrs=0, rls=0, rrs=0;
	int lastRot=0;
	if(root->left->left!=NULL){
		root->left->left->CountSubtreeBranchesAndDepth(llb, lls, 3, true);
		root->left->right->CountSubtreeBranchesAndDepth(lrb, lrs, 3, true);
		lb=llb+lrb+2;
		ls=lls+lrs+4;
		}
	if(root->left->next->left!=NULL){
		root->left->next->left->CountSubtreeBranchesAndDepth(mlb, mls, 3, true);
		root->left->next->right->CountSubtreeBranchesAndDepth(mrb, mrs, 3, true);
		mb=mlb+mrb+2;
		ms=mls+mrs+4;
		}
	if(root->right->left!=NULL){
		root->right->left->CountSubtreeBranchesAndDepth(rlb, rls, 3, true);
		root->right->right->CountSubtreeBranchesAndDepth(rrb, rrs, 3, true);
		rb=rlb+rrb+2;
		rs=rls+rrs+4;
		}
/*	
	int dl=0, dm=0, dr=0;
	root->left->CalcDepth(dl);
	root->left->next->CalcDepth(dm);
	root->right->CalcDepth(dr);
*/
	do{
		int cur=ls+ms+rs+3;
		int rotLeft=(lls-llb+lrs-lrb+2+ms+mb+rs+rb+5);
		int rotMid=(mls-mlb+mrs-mrb+2+ls+lb+rs+rb+5);
		int rotRight=(rls-rlb+rrs-rrb+2+ms+mb+ls+lb+5);
		
		if(cur<=rotLeft&&cur<=rotMid&&cur<=rotRight) return;
		else if(rotLeft<rotMid&&rotLeft<rotRight){
			RotateNodesAtRoot(root->left);
			lastRot=1;
			}
		else if(rotMid<rotLeft&&rotMid<rotRight){
			RotateNodesAtRoot(root->left->next);
			lastRot=2;
			}
		else if(rotRight<cur){
			RotateNodesAtRoot(root->right);
			lastRot=3;
			}



		lb=mb=rb=ls=ms=rs=llb=lrb=mlb=mrb=rlb=rrb=lls=lrs=mls=mrs=rls=rrs=0;
	
		if(root->left->left!=NULL){
			root->left->left->CountSubtreeBranchesAndDepth(llb, lls, 3, true);
			root->left->right->CountSubtreeBranchesAndDepth(lrb, lrs, 3, true);
			lb=llb+lrb+2;
			ls=lls+lrs+4;
			}
		if(root->left->next->left!=NULL){
			root->left->next->left->CountSubtreeBranchesAndDepth(mlb, mls, 3, true);
			root->left->next->right->CountSubtreeBranchesAndDepth(mrb, mrs, 3, true);
			mb=mlb+mrb+2;
			ms=mls+mrs+4;
			}
		if(root->right->left!=NULL){
			root->right->left->CountSubtreeBranchesAndDepth(rlb, rls, 3, true);
			root->right->right->CountSubtreeBranchesAndDepth(rrb, rrs, 3, true);
			rb=rlb+rrb+2;
			rs=rls+rrs+4;
			}
/*		
		root->left->CalcDepth(dl);
		root->left->next->CalcDepth(dm);
		root->right->CalcDepth(dr);		



*/		}while(1);
	}

void Tree::SwapAndFreeNodes(TreeNode *cop){
	assert(cop->left);//only swap internal nodes
	int tofree=cop->nodeNum;
	//we need to actually swap the memory addresses of the nodes in the allnodes array so that all other node pointers in the
	//tree stay correct
	if(allNodes[tofree]->attached){
		//find a node to swap with
		int unused=FindUnusedNode(numTipsTotal+1);
		TreeNode *tempnode=allNodes[unused];
		//swap the adresses of the nodes
		allNodes[unused]=allNodes[tofree];
		allNodes[tofree]=tempnode;
		//now adjust the nodeNums and claIndeces
		int temp=allNodes[unused]->nodeNum;
		allNodes[unused]->nodeNum=allNodes[tofree]->nodeNum;
		allNodes[tofree]->nodeNum=temp;
		
		MakeNodeDirty(allNodes[unused]);
		MakeNodeDirty(allNodes[tofree]);
		/*
		temp=allNodes[unused]->claIndexDown;
		allNodes[unused]->claIndexDown=allNodes[tofree]->claIndexDown;
		allNodes[tofree]->claIndexDown=temp;
		temp=allNodes[unused]->claIndexUL;
		allNodes[unused]->claIndexUL=allNodes[tofree]->claIndexUL;
		allNodes[tofree]->claIndexUL=temp;
		temp=allNodes[unused]->claIndexUR;
		allNodes[unused]->claIndexUR=allNodes[tofree]->claIndexUR;
		allNodes[tofree]->claIndexUR=temp;
		*/
		//set the nodes to dirty
//		assert(0);
//		allNodes[tofree]->claIndex=claMan->SetDirty(allNodes[tofree]->nodeNum, allNodes[tofree]->claIndex, true);
//		allNodes[unused]->claIndex=claMan->SetDirty(allNodes[unused]->nodeNum, allNodes[unused]->claIndex, true);
		allNodes[unused]->attached=true;
		allNodes[tofree]->attached=true;//actual its not attached, but we need to mark it as such so it isn't used as a connector
		}
	else//this is odd, but if a node will need to be used to 
		//mimic nodenums in the subtree, but was already unattached,
		//we need to mark it as attached so that it isn't used for
		//some other purpose.
		allNodes[tofree]->attached=true;
			
	if(cop->left->left) SwapAndFreeNodes(cop->left);
	if(cop->right->left) SwapAndFreeNodes(cop->right);	
	}

double Tree::SumSiteLikes(const double *cla, const int *underflow_mult){
	
	int nSites=data->NChar();
	const int *countit=data->GetCounts();
	double siteL, Lk, totalL=0.0;
	
#ifdef UNIX
	madvise((void*)cla, sizeof(double)*16*nSites, MADV_SEQUENTIAL);
	madvise((void*)underflow_mult, sizeof(int)*nSites, MADV_SEQUENTIAL);
#endif
	if(mod->ProportionInvariant()==0.0){
		if(mod->NRateCats()==1){//no invariants or gamma
		for( int k = 0; k < nSites; k++ ){
				Lk =  mod->Pi(0) * cla[0] + mod->Pi(1) * cla[1] + mod->Pi(2) * cla[2] + mod->Pi(3) * cla[3];
				cla+=4;
				siteL = (log( Lk ) - underflow_mult[k]);
				totalL += (  *countit++ *  siteL);
				}
			}
		else{//gamma, no invariants
			for( int k = 0; k < nSites; k++ ){
				Lk =  mod->Pi(0) * cla[0] + mod->Pi(1) * cla[1] + mod->Pi(2) * cla[2] + mod->Pi(3) * cla[3];
				Lk +=  mod->Pi(0) * cla[4] + mod->Pi(1) * cla[5] + mod->Pi(2) * cla[6] + mod->Pi(3) * cla[7];
				Lk +=  mod->Pi(0) * cla[8] + mod->Pi(1) * cla[9] + mod->Pi(2) * cla[10] + mod->Pi(3) * cla[11];
				Lk +=  mod->Pi(0) * cla[12] + mod->Pi(1) * cla[13] + mod->Pi(2) * cla[14] + mod->Pi(3) * cla[15];
				cla+=16;
				//this is hard coded for 4 equal sized rate cats
				siteL = (log( Lk*.25 ) - underflow_mult[k]);
				totalL += (  *countit * siteL);
				countit++;
				}
			}		
		}
	else {
		double prI=mod->ProportionInvariant();
		int lastConst=data->LastConstant();
		const int *conBases=data->GetConstBases();
		
		if(mod->NRateCats()==1){//invariants without gamma
			assert(0);
			//this isn't valid :mod->Pi(conBases[k]), because the con bases are coded as 1 2 4 8 for amiguity
			for( int k = 0; k < nSites; k++ ){
				Lk =  mod->Pi(0) * cla[0] + mod->Pi(1) * cla[1] + mod->Pi(2) * cla[2] + mod->Pi(3) * cla[3];
				cla+=4;
				if(k > lastConst){
					siteL = log( Lk * (1.0-prI)) - underflow_mult[k];
					totalL += ( *countit++ * siteL);
					}
				else{
					siteL = log( Lk * (1.0-prI) + (prI * mod->Pi(conBases[k])) * exp((double)underflow_mult[k]));
					totalL += ( *countit++ * (siteL + underflow_mult[k]));
					}
				}
			}
		else{//gamma and invariants
			double scaledGammaProp=0.25 * (1.0-prI);
	//		ofstream sites("sites.log");
			for( int k = 0; k <= lastConst; k++ ){
				Lk = mod->Pi(0)*(cla[0] + cla[4] + cla[8] + cla[12]) + mod->Pi(1)*(cla[1] + cla[5] + cla[9] + cla[13]) + mod->Pi(2)*(cla[2] + cla[6] + cla[10] + cla[14]) + mod->Pi(3)*(cla[3] + cla[7] + cla[11] + cla[15]);	
				cla+=16;
				double btot=0.0;
				if(conBases[k]&1) btot+=mod->Pi(0);
				if(conBases[k]&2) btot+=mod->Pi(1);
				if(conBases[k]&4) btot+=mod->Pi(2);
				if(conBases[k]&8) btot+=mod->Pi(3);
				if(underflow_mult[k]==0)
					siteL = log( Lk * scaledGammaProp + (prI * btot));
				else 
					siteL = log( Lk * scaledGammaProp + (prI * btot) * exp((double)underflow_mult[k]));
				totalL += ( *countit++ * (siteL - underflow_mult[k]));
				assert(totalL<0);
				
	//			sites << siteL << "\t" << underflow_mult[k] << endl;
				}			
			
			for( int k = lastConst+1; k < nSites; k++ ){
				Lk = mod->Pi(0)*(cla[0] + cla[4] + cla[8] + cla[12]) + mod->Pi(1)*(cla[1] + cla[5] + cla[9] + cla[13]) + mod->Pi(2)*(cla[2] + cla[6] + cla[10] + cla[14]) + mod->Pi(3)*(cla[3] + cla[7] + cla[11] + cla[15]);	
				cla+=16;
				siteL = log( Lk * scaledGammaProp) - underflow_mult[k];
				totalL += ( *countit++ * siteL);
				assert(totalL<0);
	//			sites << siteL << "\t" << underflow_mult[k] << endl;
				}
			}
		}
	return totalL;
	}

void Tree::CalcBipartitions(){
	root->CalcBipartition();
	}
	
void Tree::OutputBipartitions(){
	ofstream out("biparts.log", ios::app);
	root->OutputBipartition(out);
	}

void Tree::SetDistanceBasedBranchLengthsAroundNode(TreeNode *nd){
	double D1, D2, D3, k1, k2, k3, k4, a, b, c;
	TreeNode *T1, *T2, *T3, *T4;

	FindNearestTerminalUp(nd->left, T1, k1);
	FindNearestTerminalUp(nd->right, T2, k2);
	FindNearestTerminalsDown(nd->anc, nd, T3, T4, k3, k4);
//	FindNearestTerminalUp(nd->, T2, k2);

	if(k4<k3){
		T3=T4;
		k3=k4;
		}

	D1=CalculatePDistance(T1->tipData, T2->tipData, data->NChar())/(1.0-mod->ProportionInvariant()) - k1 -k2;
	D2=CalculatePDistance(T1->tipData, T3->tipData, data->NChar())/(1.0-mod->ProportionInvariant()) - k1 -k3;
	D3=CalculatePDistance(T2->tipData, T3->tipData, data->NChar())/(1.0-mod->ProportionInvariant()) - k2 -k3;

	b=(D3-D2+D1)*0.5;
	if(b < min_brlen) b=min_brlen;
	a=D1-b;
	if(a < min_brlen) a=min_brlen;
	c=D2-a;
	if(c < min_brlen) c=min_brlen;
	
	nd->left->dlen=a;
	nd->right->dlen=b;
	nd->dlen=c;

	SweepDirtynessOverTree(nd->left);
	SweepDirtynessOverTree(nd);
	SweepDirtynessOverTree(nd->right);
	}

void Tree::FindNearestTerminalUp(TreeNode *start, TreeNode *&term, double &dist){
	dist=999999.9;
	int nodeDist=9999;
	sprRange.clear();
	sprRange.setseed(start->nodeNum);
	int range=10;
    for(int i = 0;i<range;i++){
      int j =  sprRange.total;
		for(int k=0; k < j; k++){
			if(sprRange.front[k]==i){
				TreeNode *cur=allNodes[sprRange.element[k]];
				if(cur->left!=NULL){
				    sprRange.addelement(cur->left->nodeNum, i+1, sprRange.pathlength[k]+cur->left->dlen);
				    sprRange.addelement(cur->right->nodeNum, i+1, sprRange.pathlength[k]+cur->right->dlen);
				    }
				else{
					//if(sprRange.pathlength[k]<dist){
					if(sprRange.front[k]<nodeDist){
						nodeDist=sprRange.front[k];
						term=cur;
						dist=sprRange.pathlength[k];
						}
					}
				}
		    }
		}
	}

void Tree::FindNearestTerminalsDown(TreeNode *start, TreeNode *from, TreeNode *&term1, TreeNode *&term2, double &dist1, double &dist2){
	dist1=dist2=999999.9;
	int nodeDist1=9999, nodeDist2=9999;
	sprRange.clear();
	if(from==start->left) sprRange.setseed(start->right->nodeNum, start->right->dlen);
	else sprRange.setseed(start->left->nodeNum, start->left->dlen);
	int range=10;
    for(int i = 0;i<range;i++){
      int j =  sprRange.total;
		for(int k=0; k < j; k++){
			if(sprRange.front[k]==i){
				TreeNode *cur=allNodes[sprRange.element[k]];
				if(cur->left!=NULL){
				    sprRange.addelement(cur->left->nodeNum, i+1, sprRange.pathlength[k]+cur->left->dlen);
				    sprRange.addelement(cur->right->nodeNum, i+1, sprRange.pathlength[k]+cur->right->dlen);
				    }
				else{
					//if(sprRange.pathlength[k]<dist1){
					if(sprRange.front[k]<nodeDist1){
						nodeDist1=sprRange.front[k];
						term1=cur;
						dist1=sprRange.pathlength[k];
						}
					}
				}
		    }
		}

	sprRange.clear();
	if(start->anc != NULL){		
		sprRange.setseed(start->anc->nodeNum, start->dlen);
		for(int i = 0;i<range;i++){
	      int j =  sprRange.total;
			for(int k=0; k < j; k++){
				if(sprRange.front[k]==i){
					TreeNode *cur=allNodes[sprRange.element[k]];
					if(cur->left!=NULL){
					    if(cur->left!=from->anc) sprRange.addelement(cur->left->nodeNum, i+1, sprRange.pathlength[k]+cur->left->dlen);
					    if(cur->right!=from->anc) sprRange.addelement(cur->right->nodeNum, i+1, sprRange.pathlength[k]+cur->right->dlen);
					    }
					else{
						//if(sprRange.pathlength[k]<dist2){
						if(sprRange.front[k]<nodeDist2){
							nodeDist2=sprRange.front[k];
							term2=cur;
							dist2=sprRange.pathlength[k];
							}
						}
					if(cur->anc) sprRange.addelement(cur->anc->nodeNum, i+1, sprRange.pathlength[k]+cur->dlen);
					else sprRange.addelement(cur->left->next->nodeNum, i+1, sprRange.pathlength[k]+cur->left->next->dlen);
					}
			    }
			}
		}
	else{
		if(from!=start->left->next) sprRange.setseed(start->left->next->nodeNum, start->left->next->dlen);
		else sprRange.setseed(start->right->nodeNum, start->right->dlen);
		int range=10;
	    for(int i = 0;i<range;i++){
	      int j =  sprRange.total;
			for(int k=0; k < j; k++){
				if(sprRange.front[k]==i){
					TreeNode *cur=allNodes[sprRange.element[k]];
					if(cur->left!=NULL){
					    sprRange.addelement(cur->left->nodeNum, i+1, sprRange.pathlength[k]+cur->left->dlen);
					    sprRange.addelement(cur->right->nodeNum, i+1, sprRange.pathlength[k]+cur->right->dlen);
					    }
					else{
						//if(sprRange.pathlength[k]<dist2){
						if(sprRange.front[k]<nodeDist2){
							nodeDist2=sprRange.front[k];
							term2=cur;
							dist2=sprRange.pathlength[k];
							}
						}
					}
			    }
			}
		}
	assert(term1 != term2);
	}

void Tree::OptimizeBranchesAroundNode(TreeNode *nd, double optPrecision, int subtreeNode){
	//depricated
	assert(0);
	//this function will optimize the three branches (2 descendents and one anc) connecting
	//to it.  It assumes that everything that is dirty has been marked so.
	//by default there is only a single optimization pass over the three nodes
/*	double precision1, precision2;

	if(subtreeNode==0) SetAllTempClasDirty();
	
	precision1=optPrecision;// * 0.5;
	if(optPrecision > .2) precision2=0.0;
	else precision2=precision1 * 0.5;
	
	if(nd != root){
		BrentOptimizeBranchLength(precision1, nd, false);
		BrentOptimizeBranchLength(precision1, nd->left, false);
		BrentOptimizeBranchLength(precision1, nd->right, false);
		}
	else{
		BrentOptimizeBranchLength(precision1, nd->left, false);
		BrentOptimizeBranchLength(precision1, nd->left->next, false);
		BrentOptimizeBranchLength(precision1, nd->right, false);	
		}
*/
/*	
	if(precision2 > 0){
		//if were're doing multiple optimization passes, only this stuff needs to be set dirty
		claMan->SetDirty(nd->nodeNum, nd->claIndex, true);
		claMan->SetTempDirty(nd->nodeNum, true);
		if(nd != root) claMan->SetTempDirty(nd->anc->nodeNum, true);

		if(nd != root){
			BrentOptimizeBranchLength(precision2, nd, false);
			BrentOptimizeBranchLength(precision2, nd->left, false);
			BrentOptimizeBranchLength(precision2, nd->right, false);
			}
		else {
			BrentOptimizeBranchLength(precision2, nd->left, false);
			BrentOptimizeBranchLength(precision2, nd->left->next, false);
			BrentOptimizeBranchLength(precision2, nd->right, false);			
			}
		}
*/		
/*	//these must be called after all optimization passes are done around this node
	TraceDirtynessToRoot(nd);
	if(subtreeNode==0)
		SetAllTempClasDirty();
	else SetTempClasDirtyWithinSubtree(subtreeNode);
*/	}

double Tree::ScoreAtSubtree(int subtreeNode){
	//depricated
	assert(0);
/*	TreeNode *effectiveRoot=allNodes[subtreeNode]->left;
	
	if(effectiveRoot->anc->anc){
		if(mod->NRateCats()==1);
//			TraceLikeUpFromRoot(effectiveRoot->anc->anc, effectiveRoot->anc, true);
//		else TraceLikeUpFromRootRateHet(effectiveRoot->anc->anc, effectiveRoot->anc, true);
		}	
	double score;
	if(mod->NRateCats()==1)
		score=BranchLike(effectiveRoot);
	else score=BranchLikeRateHet(effectiveRoot);
//	claMan->SetTempDirty(-1, true);
	return score;	
*/	
	return 0;
	}

void Tree::RerootHere(int newroot){
	//DJZ 1-5-05 adding functionality to adjust the direction of existing clas
	//so that they are still valid in the new context, rather than just dirtying everything
	//REMEMBER that the mutation_type of the individual this is called for needs to be 
	// "|= rerooted" so that the topo numbers are updated properly

	TreeNode *nroot=allNodes[newroot];

	TreeNode *prevnode=nroot;
	TreeNode *curnode=nroot->anc;
	TreeNode *nextnode=nroot->anc->anc;

	//first trace down to the old root and fix all the blens
	//Each branch with take the length of its descendent on that path
	//this will be easiest recursively
	nroot->FlipBlensToRoot(0);
	
	//now take the new root's current ancestor and make it the middle des
	//note that the existing cla directions at this node are still valid
	nroot->left->next=curnode;
	curnode->next=nroot->right;
	nroot->right->prev=curnode;
	curnode->prev=nroot->left;
	
	//this needs to work slightly differently if the old root is the anc of the new one
	if(curnode!=root){
		if(prevnode==curnode->left){
			curnode->left=curnode->anc;
			curnode->AdjustClasForReroot(UPLEFT);
			}
		else{
			curnode->right=curnode->anc;
			curnode->AdjustClasForReroot(UPRIGHT);
			}
		
		curnode->left->next=curnode->right;
		curnode->left->prev=NULL;
		curnode->right->prev=curnode->left;
		curnode->right->next=NULL;

		prevnode=curnode;
		curnode=nextnode;
		nextnode=nextnode->anc;
		}
		
	curnode->anc=prevnode;
	nroot->anc=NULL;
			
	while(curnode!=root){
		if(prevnode==curnode->left){
			curnode->left=nextnode;
			curnode->AdjustClasForReroot(UPLEFT);
			}
		else{
			curnode->right=nextnode;
			curnode->AdjustClasForReroot(UPRIGHT);
			}
			
		curnode->left->next=curnode->right;
		curnode->left->prev=NULL;
		curnode->right->prev=curnode->left;
		curnode->right->next=NULL;
				
		curnode->anc=prevnode;
		
		prevnode=curnode;
		curnode=nextnode;
		nextnode=nextnode->anc;
		}
	
	//now deal with the old root, which is now curnode
	if(prevnode==curnode->left){
		curnode->left=curnode->right->prev;
		curnode->left->prev=NULL;
		curnode->AdjustClasForReroot(UPLEFT);
		}
	else if(prevnode==curnode->left->next){
		curnode->left->next=curnode->right;
		curnode->right->prev=curnode->left;
		//clas don't need to be adjusted in this case
		}
	else{
		curnode->right=curnode->left->next;
		curnode->right->next=NULL;		
		curnode->AdjustClasForReroot(UPRIGHT);
		}
		
	curnode->anc=prevnode;
	
	//now we just need to make the newroot node0 and swap it with the old root, which means moving the
	//_data_ to node 0, not just swapping the memory addresses
	
	SwapNodeDataForReroot(nroot);

	root->CheckTreeFormation();
	}

void Tree::SwapNodeDataForReroot(TreeNode *nroot){
	TreeNode tempold;
	tempold.left=root->left;
	tempold.right=root->right;
	tempold.next=root->next;
	tempold.prev=root->prev;
	//note that we need to watch out here if the new root is currently the anc of the old root
	if(root->anc==nroot) tempold.anc=root;
	else tempold.anc=root->anc;
	tempold.dlen=root->dlen;
	tempold.claIndexDown=root->claIndexDown;
	tempold.claIndexUL=root->claIndexUL;
	tempold.claIndexUR=root->claIndexUR;
	
	TreeNode tempnew;
	tempnew.left=nroot->left;
	tempnew.right=nroot->right;
	tempnew.next=nroot->next;
	tempnew.prev=nroot->prev;
	tempnew.anc=nroot->anc;
	tempnew.dlen=nroot->dlen;	
	tempnew.claIndexDown=nroot->claIndexDown;
	tempnew.claIndexUL=nroot->claIndexUL;
	tempnew.claIndexUR=nroot->claIndexUR;

	root->left=tempnew.left;
	root->left->anc=root;
	root->right=tempnew.right;
	root->right->anc=root;
	root->left->next->anc=root;
	root->prev=root->next=NULL;
	root->anc=NULL;
	root->dlen=-1;
	root->claIndexDown=tempnew.claIndexDown;
	root->claIndexUL=tempnew.claIndexUL;
	root->claIndexUR=tempnew.claIndexUR;
	
	nroot->left=tempold.left;
	nroot->left->anc=nroot;
	nroot->right=tempold.right;
	nroot->next=tempold.next;
	if(nroot->next) nroot->next->prev=nroot;
	nroot->prev=tempold.prev;
	if(nroot->prev) nroot->prev->next=nroot;
	nroot->right->anc=nroot;
	nroot->anc=tempold.anc;
	nroot->claIndexDown=tempold.claIndexDown;
	nroot->claIndexUL=tempold.claIndexUL;
	nroot->claIndexUR=tempold.claIndexUR;
	
	if(nroot->anc->left==root){
		nroot->anc->left=nroot;
		nroot->prev=NULL;
		nroot->next=nroot->anc->right;
		nroot->next->prev=nroot;
		}
	else if(nroot->anc->right==root){
		nroot->anc->right=nroot;
		nroot->next=NULL;
		nroot->prev=nroot->anc->left;
		nroot->prev->next=nroot;
		}
	else{
		nroot->anc->left->next=nroot;
//		nroot->next=NULL;
		nroot->prev=nroot->anc->left;
//		nroot->prev->next=nroot;		
		}
	nroot->dlen=tempold.dlen;
	}

	
void Tree::MakeNodeDirty(TreeNode *nd){
	nd->claIndexDown=claMan->SetDirty(nd->claIndexDown);
	nd->claIndexUL=claMan->SetDirty(nd->claIndexUL);
	nd->claIndexUR=claMan->SetDirty(nd->claIndexUR);
	}
	
void Tree::RemoveTempClaReservations(){
	if(memLevel > 1){
		for(int i=numTipsTotal+1;i<numNodesTotal;i++){
			claMan->ClearTempReservation(allNodes[i]->claIndexDown);
			}
		}
	
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		claMan->ClearTempReservation(allNodes[i]->claIndexUR);
		}
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		claMan->ClearTempReservation(allNodes[i]->claIndexUL);
		}
	}

void Tree::ReclaimUniqueClas(){
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		if(claMan->GetNumAssigned(allNodes[i]->claIndexDown) == 1){
			claMan->ReclaimSingleCla(allNodes[i]->claIndexDown);
			}
		if(claMan->GetNumAssigned(allNodes[i]->claIndexUL) == 1){
			claMan->ReclaimSingleCla(allNodes[i]->claIndexUL);
			}
		if(claMan->GetNumAssigned(allNodes[i]->claIndexUR) == 1){
			claMan->ReclaimSingleCla(allNodes[i]->claIndexUR);
			}
		}
	}

void Tree::MarkUpwardClasToReclaim(int subtreeNode){
	//if we are somewhat low on clas, mark some reclaimable that were 
	//used tracing the likelihood upward for blen optimization
	assert(0);
	if(subtreeNode==0){
/*		if(memLevel==2){
			if(allNodes[0]->claIndexUL > 0)
				claMan->MarkReclaimable(allNodes[0]->claIndexUL, 2);			
			if(allNodes[0]->claIndexUR > 0)
				claMan->MarkReclaimable(allNodes[0]->claIndexUR, 2);		
			}
*/		for(int i=numTipsTotal+1;i<numNodesTotal;i++){
//			claMan->MarkReclaimable(allNodes[i]->claIndexUL, 2, false);
//			claMan->MarkReclaimable(allNodes[i]->claIndexUR, 2, false);
			}
		}
	else{
		for(int i=numTipsTotal+1;i<numNodesTotal;i++){
			if((allNodes[i]->nodeNum != subtreeNode) && (allNodes[i]->nodeNum != allNodes[subtreeNode]->anc->nodeNum)){
				if(allNodes[i]->claIndexUL > 0){
//					claMan->MarkReclaimable(allNodes[i]->claIndexUL, 2, false);
					}
				if(allNodes[i]->claIndexUR > 0){
//					claMan->MarkReclaimable(allNodes[i]->claIndexUR, 2, false);
					}
				}
			}
		}
	}

void Tree::MarkDownwardClasToReclaim(int subtreeNode){
	//if we're calling this, we must really be desperate for clas
	//this should only be called after the tree has been scored
	assert(0);

	if(subtreeNode==0){
		if(memLevel<3){
			for(int i=numTipsTotal+1;i<numNodesTotal;i++){
//				claMan->MarkReclaimable(allNodes[i]->claIndexDown, 1);
				}
			}
		else{
			for(int i=numTipsTotal+1;i<numNodesTotal;i++){
//				claMan->MarkReclaimable(allNodes[i]->claIndexDown, 1, false);
				}
			}
		}
	else{
		return;  //I think that this is safe, since in general many fewer node will be necessary in subtree mode
		for(int i=numTipsTotal+1;i<numNodesTotal;i++){
			if((allNodes[i]->nodeNum != subtreeNode) && (allNodes[i]->nodeNum != allNodes[subtreeNode]->anc->nodeNum)){
				if(allNodes[i]->claIndexUL > 0){
//					claMan->MarkReclaimable(allNodes[i]->claIndexUL, 1);
					}
				}
			}
		}
	}

void Tree::MarkClasNearTipsToReclaim(int subtreeNode){
	assert(0);
	if(subtreeNode==0){
		for(int i=1;i<numTipsTotal;i++){
//			claMan->MarkReclaimable(allNodes[i]->anc->claIndexDown, 1, false);
			}
		}
	else{
		return;  //I think that this is safe, since in general many fewer node will be necessary in subtree mode
		for(int i=numTipsTotal+1;i<numNodesTotal;i++){
			if((allNodes[i]->nodeNum != subtreeNode) && (allNodes[i]->nodeNum != allNodes[subtreeNode]->anc->nodeNum)){
				if(allNodes[i]->claIndexUL > 0){
//					claMan->MarkReclaimable(allNodes[i]->claIndexUL, 1);
					}
				}
			}
		}
	}

void Tree::OutputFirstClaAcrossTree(ofstream &deb, TreeNode *nd){
	int site=0;
	int index=16*site;

	
	if(nd->left!=NULL && claMan->IsDirty(nd->claIndexDown) == false)
		deb << nd->nodeNum << "\t0\t" << nd->claIndexDown << "\t" << claMan->GetCla(nd->claIndexDown)->arr[index] << "\t" << claMan->GetCla(nd->claIndexDown)->underflow_mult[site] <<"\n";
	if(nd->left!=NULL && claMan->IsDirty(nd->claIndexUL) == false)
		deb << nd->nodeNum << "\t1\t" << nd->claIndexUL << "\t" << claMan->GetCla(nd->claIndexUL)->arr[index] << "\t" << claMan->GetCla(nd->claIndexUL)->underflow_mult[site] <<"\n";
	if(nd->left!=NULL && claMan->IsDirty(nd->claIndexUR) == false)
		deb << nd->nodeNum << "\t2\t" << nd->claIndexUR << "\t" << claMan->GetCla(nd->claIndexUR)->arr[index] << "\t" << claMan->GetCla(nd->claIndexUR)->underflow_mult[site] <<"\n";
				
	if(nd->left!=NULL)
		OutputFirstClaAcrossTree(deb, nd->left);
	if(nd->next!=NULL)
		OutputFirstClaAcrossTree(deb, nd->next);
	}

void Tree::CountNumReservedClas(int &clean, int &tempRes, int&res){
	clean=0;
	tempRes=0;
	res=0;
	
	if(claMan->IsDirty(allNodes[0]->claIndexDown)==false){
		clean++;
		res += (claMan->IsClaReserved(allNodes[0]->claIndexDown));
		tempRes += (claMan->IsClaTempReserved(allNodes[0]->claIndexDown));
		}
	if(claMan->IsDirty(allNodes[0]->claIndexUL)==false){
		clean++;
		res += (claMan->IsClaReserved(allNodes[0]->claIndexUL));
		tempRes += (claMan->IsClaTempReserved(allNodes[0]->claIndexUL));
		}
	if(claMan->IsDirty(allNodes[0]->claIndexUR)==false){
		clean++;
		res += (claMan->IsClaReserved(allNodes[0]->claIndexUR));
		tempRes += (claMan->IsClaTempReserved(allNodes[0]->claIndexUR));
		}
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		if(claMan->IsDirty(allNodes[i]->claIndexDown)==false){
			clean++;
			res += (claMan->IsClaReserved(allNodes[i]->claIndexDown));
			tempRes += (claMan->IsClaTempReserved(allNodes[i]->claIndexDown));
			}
		if(claMan->IsDirty(allNodes[i]->claIndexUL)==false){
			clean++;
			res += (claMan->IsClaReserved(allNodes[i]->claIndexUL));
			tempRes += (claMan->IsClaTempReserved(allNodes[i]->claIndexUL));
			}
		if(claMan->IsDirty(allNodes[i]->claIndexUR)==false){
			clean++;
			res += (claMan->IsClaReserved(allNodes[i]->claIndexUR));
			tempRes += (claMan->IsClaTempReserved(allNodes[i]->claIndexUR));
			}
		}
	}

void Tree::SetupClasForSubtreeMode(int subtreeNode){
	TreeNode *subnode=allNodes[subtreeNode];
	
	claMan->ReserveCla(subnode->claIndexDown, false);
	claMan->ReserveCla(subnode->claIndexUL, false);
	claMan->ReserveCla(subnode->claIndexUR, false);
	
	if(subnode->anc != root){
		if(subnode->anc->left==subnode) claMan->ReserveCla(subnode->anc->claIndexUL, false);
		else if(subnode->anc->right==subnode) claMan->ReserveCla(subnode->anc->claIndexUR, false);
		}
	
	DirtyNodesOutsideOfSubtree(root, subtreeNode);
	}
	
void Tree::DirtyNodesOutsideOfSubtree(TreeNode *nd, int subtreeNode){

	if(nd != root){
		claMan->ReclaimSingleCla(nd->claIndexDown);
		claMan->ReclaimSingleCla(nd->claIndexUL);
		claMan->ReclaimSingleCla(nd->claIndexUR);
		}
	
	if(nd->left->left != NULL && nd->left->nodeNum != subtreeNode && nd->left->nodeNum != allNodes[subtreeNode]->anc->nodeNum){
		DirtyNodesOutsideOfSubtree(nd->left, subtreeNode);
		}
	if(nd->right->left != NULL && nd->right->nodeNum != subtreeNode && nd->right->nodeNum != allNodes[subtreeNode]->anc->nodeNum){
		DirtyNodesOutsideOfSubtree(nd->right, subtreeNode);
		}
	if(nd->anc==NULL && nd->left->next->left != NULL && nd->left->next->nodeNum != subtreeNode && nd->left->next->nodeNum != allNodes[subtreeNode]->anc->nodeNum){
		DirtyNodesOutsideOfSubtree(nd->left->next, subtreeNode);
		}
	}

void Tree::OutputValidClaIndeces(){
	ofstream cla("claind.log");
	if(claMan->IsDirty(allNodes[0]->claIndexDown)==false){
		cla << "0\t" << allNodes[0]->claIndexDown << "\t" << claMan->GetNumAssigned(allNodes[0]->claIndexDown) << "\t" << claMan->GetReclaimLevel(allNodes[0]->claIndexDown) << "\t" << claMan->IsClaReserved(allNodes[0]->claIndexDown) <<"\n";
		}
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		cla << i << "\t" << allNodes[i]->claIndexDown << "\t" << claMan->GetNumAssigned(allNodes[i]->claIndexDown) << "\t" << claMan->GetReclaimLevel(allNodes[i]->claIndexDown) << "\t" << claMan->IsClaReserved(allNodes[i]->claIndexDown) << "\n";
		}
	cla.close();
	}

void Tree::GetInternalStateString(char *string, int nodeNum){
	assert(0);
//	Score(nodeNum);
//	InferStatesFromCla(string, claMan->GetTempCla()->arr, data->NChar());	
	}

void Tree::InferAllInternalStateProbs(char *ofprefix){
	char filename[80];
	sprintf(filename, "%s.internalstates.log", ofprefix);
	ofstream out(filename);
	out.precision(5);
	RecursivelyCalculateInternalStateProbs(root, out);
	out.close();
	}

void Tree::RecursivelyCalculateInternalStateProbs(TreeNode *nd, ofstream &out){
	if(nd->left != NULL) RecursivelyCalculateInternalStateProbs(nd->left, out);
	if(nd->next != NULL) RecursivelyCalculateInternalStateProbs(nd->next, out);
	
	if(nd->left != NULL){
		int wholeTreeIndex=ConditionalLikelihoodRateHet(ROOT, nd, true);
		vector<InternalState *> *stateProbs=InferStatesFromCla(claMan->GetCla(wholeTreeIndex)->arr, data->NChar());

		char subtreeString[5000];
		nd->MakeNewickForSubtree(subtreeString);		
		out << "node " << nd->nodeNum << "\t" << subtreeString << "\t";
		char *loc=subtreeString;
		NxsString temp;
		
		while(*loc){
			if(isdigit(*loc) == false) out << *loc++;
			else{
				while(isdigit(*loc))
					temp += *loc++;
				out << data->TaxonLabel(atoi(temp.c_str())-1);
				temp="";
				}
			}
		out << "\n";
		
		for(int s=0;s<data->GapsIncludedNChar();s++){
			out << s+1 << "\t";
			if(data->Number(s) > -1)
				(*stateProbs)[data->Number(s)]->Output(out);
			else out << "Entirely uninformative character (gaps,N's or ?'s)\n";
			}
		
		claMan->ClearTempReservation(wholeTreeIndex);
		claMan->DecrementCla(wholeTreeIndex);
		
		for(vector<InternalState*>::iterator delit=stateProbs->begin();delit!=stateProbs->end();delit++){
			delete *(delit);
			}
		delete stateProbs;
		}
	}
	
void Tree::ClaReport(ofstream &cla){
	int totDown=0;
	int totUL=0;
	int totUR=0;
	
	cla << "root\t" << claMan->GetReclaimLevel(root->claIndexDown) << "\t" << claMan->GetNumAssigned(root->claIndexDown)<< "\t" << claMan->GetClaNumber(root->claIndexDown);
	cla << "\n\t" << claMan->GetReclaimLevel(root->claIndexUL) << "\t" << claMan->GetNumAssigned(root->claIndexUL) << "\t" << claMan->GetClaNumber(root->claIndexUL);
	cla << "\n\t" << claMan->GetReclaimLevel(root->claIndexUR)  << "\t" << claMan->GetNumAssigned(root->claIndexUR)  << "\t" << claMan->GetClaNumber(root->claIndexUR) << "\n";
//	cla << "\t" << claMan->GetNumAssigned(root->claIndexDown) << "\t" << claMan->GetNumAssigned(root->claIndexUL) << "\t" << claMan->GetNumAssigned(root->claIndexUR)  << "\n";
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		TreeNode *n=allNodes[i];
	cla << i << "\t" << claMan->GetReclaimLevel(n->claIndexDown) << "\t" << claMan->GetNumAssigned(n->claIndexDown) << "\t" << claMan->GetClaNumber(n->claIndexDown);
	cla << "\n\t" << claMan->GetReclaimLevel(n->claIndexUL) << "\t" << claMan->GetNumAssigned(n->claIndexUL) << "\t" << claMan->GetClaNumber(n->claIndexUL);
	cla << "\n\t" << claMan->GetReclaimLevel(n->claIndexUR)  << "\t" << claMan->GetNumAssigned(n->claIndexUR)  << "\t" << claMan->GetClaNumber(n->claIndexUR) << "\n";
		totDown += claMan->GetReclaimLevel(n->claIndexDown);
		totUL += claMan->GetReclaimLevel(n->claIndexUL);
		totUR += claMan->GetReclaimLevel(n->claIndexUR);
		}
	cla << "tots\t" << totDown << "\t" << totUL << "\t" << totUR << endl;
//	cla.close();
	}
	
double Tree::CountClasInUse(){
	double inUse=0.0;
	
	if(claMan->IsDirty(root->claIndexDown) == false) inUse += 1.0/claMan->GetNumAssigned(root->claIndexDown);
	if(claMan->IsDirty(root->claIndexUL) == false) inUse += 1.0/claMan->GetNumAssigned(root->claIndexUL);
	if(claMan->IsDirty(root->claIndexUR) == false) inUse += 1.0/claMan->GetNumAssigned(root->claIndexUR);
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		TreeNode *n=allNodes[i];	
		if(claMan->IsDirty(n->claIndexDown) == false) inUse += 1.0/claMan->GetNumAssigned(n->claIndexDown);
		if(claMan->IsDirty(n->claIndexUL) == false) inUse += 1.0/claMan->GetNumAssigned(n->claIndexUL);
		if(claMan->IsDirty(n->claIndexUR) == false) inUse += 1.0/claMan->GetNumAssigned(n->claIndexUR);		
		}
	return inUse;
	}
	
double Tree::GetScorePartialTerminalRateHet(const CondLikeArray *partialCLA, const double *prmat, const char *Ldata){

	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *partial=partialCLA->arr;
	int *underflow_mult=partialCLA->underflow_mult;

	int nchar=data->NChar();
#ifdef UNIX
	madvise(partial, nchar*16*sizeof(double), MADV_SEQUENTIAL);
#endif

	double siteL, totallnL=0.0;
	double La, Lc, Lg, Lt;
	
	//gamma and invariants
	const int *countit=data->GetCounts();
	int lastConst=data->LastConstant();
	const int *conBases=data->GetConstBases();
	double prI=mod->ProportionInvariant();
	double scaledGammaProp=(1.0-prI) / mod->NRateCats();
	
	for(int i=0;i<nchar;i++){
		La=Lc=Lg=Lt=0.0;
		if(*Ldata > -1){ //no ambiguity
			for(int i=0;i<4;i++){
				La  += prmat[(*Ldata)+16*i] * partial[0];
				Lc  += prmat[(*Ldata+4)+16*i] * partial[1];
				Lg  += prmat[(*Ldata+8)+16*i] * partial[2];
				Lt  += prmat[(*Ldata+12)+16*i] * partial[3];
				partial += 4;
				}
			Ldata++;
			}
			
		else if(*Ldata == -4){ //total ambiguity
			for(int i=0;i<4;i++){
				La += partial[0];
				Lc += partial[1];
				Lg += partial[2];
				Lt += partial[3];
				partial += 4;
				}
			Ldata++;
			}
		else{ //partial ambiguity
			char nstates=-1 * *(Ldata++);
			for(int i=0;i<nstates;i++){
				for(int i=0;i<4;i++){
					La += prmat[(*Ldata)+16*i]  * partial[4*i];
					Lc += prmat[(*Ldata+4)+16*i] * partial[1+4*i];
					Lg += prmat[(*Ldata+8)+16*i]* partial[2+4*i];
					Lt += prmat[(*Ldata+12)+16*i]* partial[3+4*i];
					}
				Ldata++;
				}
			partial+=16;
			}
		if(mod->NoPinvInModel() == false && i<=lastConst){
			double btot=0.0;
			if(conBases[i]&1) btot+=mod->Pi(0);
			if(conBases[i]&2) btot+=mod->Pi(1);
			if(conBases[i]&4) btot+=mod->Pi(2);
			if(conBases[i]&8) btot+=mod->Pi(3);
			if(underflow_mult[i]==0)
				siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)) * scaledGammaProp + prI*btot);
			else 
				siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)) * scaledGammaProp + (prI*btot*exp((double)underflow_mult[i])));
			}
		else
			siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)) * scaledGammaProp);
		
		totallnL += (*countit++ * (log(siteL) - underflow_mult[i]));
		}
	return totallnL;
	}
	
double Tree::GetScorePartialInternalRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const double *prmat){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *CL1=childCLA->arr;
	double *partial=partialCLA->arr;
	int *underflow_mult1=partialCLA->underflow_mult;
	int *underflow_mult2=childCLA->underflow_mult;

	int nchar=data->NChar();

#ifdef UNIX
	madvise(partial, nchar*16*sizeof(double), MADV_SEQUENTIAL);
	madvise(CL1, nchar*16*sizeof(double), MADV_SEQUENTIAL);
#endif

	double siteL, totallnL=0.0;
	double La, Lc, Lg, Lt;

	//gamma and invariants
	const int *countit=data->GetCounts();
	int lastConst=data->LastConstant();
	const int *conBases=data->GetConstBases();
	double prI=mod->ProportionInvariant();
	double scaledGammaProp=(1.0-prI) / mod->NRateCats();

//	ofstream sit("sitelikes.log");

	for(int i=0;i<nchar;i++){
		La=Lc=Lg=Lt=0.0;
		for(int r=0;r<4;r++){
			int rOff=r*16;
			La += ( prmat[rOff ]*CL1[0]+prmat[rOff + 1]*CL1[1]+prmat[rOff + 2]*CL1[2]+prmat[rOff + 3]*CL1[3]) * partial[0];
			Lc += ( prmat[rOff + 4]*CL1[0]+prmat[rOff + 5]*CL1[1]+prmat[rOff + 6]*CL1[2]+prmat[rOff + 7]*CL1[3]) * partial[1];
			Lg += ( prmat[rOff + 8]*CL1[0]+prmat[rOff + 9]*CL1[1]+prmat[rOff + 10]*CL1[2]+prmat[rOff + 11]*CL1[3]) * partial[2];
			Lt += ( prmat[rOff + 12]*CL1[0]+prmat[rOff + 13]*CL1[1]+prmat[rOff + 14]*CL1[2]+prmat[rOff + 15]*CL1[3]) * partial[3];
			partial+=4;
			CL1+=4;
			}
		if(mod->NoPinvInModel() == false && i<=lastConst){
			double btot=0.0;
			if(conBases[i]&1) btot+=mod->Pi(0);
			if(conBases[i]&2) btot+=mod->Pi(1);
			if(conBases[i]&4) btot+=mod->Pi(2);
			if(conBases[i]&8) btot+=mod->Pi(3);
			if(underflow_mult1[i] + underflow_mult2[i] == 0)
				siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)) * scaledGammaProp + prI*btot);
			else 
				siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)) * scaledGammaProp + (prI*btot*exp((double)underflow_mult1[i]+underflow_mult2[i])));
			}
		else
			siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)) * scaledGammaProp);
		
		totallnL += (*countit++ * (log(siteL) - underflow_mult1[i] - underflow_mult2[i]));
//		sit << siteL << "\t" << underflow_mult1[i] << "\t" << underflow_mult2[i] << "\t" << totallnL << "\n";
		}
//	sit.close();
	return totallnL;
	}
	
	
