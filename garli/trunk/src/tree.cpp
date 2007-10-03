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

#include <algorithm>
#include <vector>
#include <list>
#include <cassert>
#ifdef UNIX
	#include <sys/mman.h>
#endif

using namespace std;

#include "defs.h"
#include "mlhky.h"
#include "clamanager.h"
#include "funcs.h"
#include "stopwatch.h"
#include "model.h"
#include "tree.h"
#include "datamatr.h"
#include "reconnode.h"

#include "utility.h"
Profiler ProfIntInt   ("ClaIntInt     ");
Profiler ProfIntTerm  ("ClaIntTerm    ");
Profiler ProfTermTerm ("ClaTermTerm   ");
Profiler ProfRescale  ("Rescale       ");
Profiler ProfScoreInt ("ScoreInt      ");
Profiler ProfScoreTerm("ScoreTerm     ");
Profiler ProfEQVectors("EQVectors     ");

FLOAT_TYPE precalcThresh[30];
FLOAT_TYPE precalcMult[30];
int precalcIncr[30] = {1, 3, 5, 7, 10, 12, 14, 17, 19, 21, 24, 26, 28, 30, 33, 35, 37, 40, 42, 44, 47, 49, 51, 53, 56, 58, 60, 63, 65, 67};

extern rng rnd;
extern bool output_tree;
extern bool uniqueSwapTried;

#ifdef VARIABLE_OPTIMIZATION
ofstream var("variable.log");
ofstream uni("unique.log");
#endif

#ifdef OUTPUT_UNIQUE_TREES
ofstream uni("unique.log");
#endif

//external global variables
extern int calcCount;
extern int optCalcs;
extern ofstream opt;
extern ofstream optsum;
extern int memLevel;
extern ModelSpecification modSpec;

//Tree static definitions
FLOAT_TYPE Tree::meanBrlenMuts;
FLOAT_TYPE Tree::alpha;
FLOAT_TYPE Tree::min_brlen;	  // branch lengths never below this value
FLOAT_TYPE Tree::max_brlen;	  
FLOAT_TYPE Tree::exp_starting_brlen;    // expected starting branch length
ClaManager *Tree::claMan;
list<TreeNode *> Tree::nodeOptVector;
const HKYData *Tree::data;
unsigned Tree::rescaleEvery;
FLOAT_TYPE Tree::treeRejectionThreshold;
vector<Constraint> Tree::constraints;
AttemptedSwapList Tree::attemptedSwaps;
FLOAT_TYPE Tree::uniqueSwapBias;
FLOAT_TYPE Tree::distanceSwapBias;

float Tree::uniqueSwapPrecalc[500];
float Tree::distanceSwapPrecalc[1000];

Bipartition *Tree::outgroup = NULL;

void InferStatesFromCla(char *states, FLOAT_TYPE *cla, int nchar);
FLOAT_TYPE CalculateHammingDistance(const char *str1, const char *str2, int nchar);
void SampleBranchLengthCurve(FLOAT_TYPE (*func)(TreeNode*, Tree*, FLOAT_TYPE, bool), TreeNode *thisnode, Tree *thistree);
FLOAT_TYPE CalculatePDistance(const char *str1, const char *str2, int nchar);
inline FLOAT_TYPE CallBranchLike(TreeNode *thisnode, Tree *thistree, FLOAT_TYPE blen, bool brak);
bool Rescale(CondLikeArray *destCLA, int nsites);
bool RescaleRateHet(CondLikeArray *destCLA, int nsites);


void Tree::SetTreeStatics(ClaManager *claMan, const HKYData *data, const GeneralGamlConfig *conf){
	Tree::claMan=claMan;
	Tree::data=data;
#ifdef SINGLE_PRECISION_FLOATS
	Tree::rescaleEvery=6;
#else
	Tree::rescaleEvery=16;
#endif
	Tree::uniqueSwapBias = conf->uniqueSwapBias;
	Tree::distanceSwapBias = conf->distanceSwapBias;
	for(int i=0;i<500;i++){
		Tree::uniqueSwapPrecalc[i] = (float) pow(Tree::uniqueSwapBias, i);
		if(Tree::uniqueSwapPrecalc[i] != Tree::uniqueSwapPrecalc[i]) Tree::uniqueSwapPrecalc[i]=0.0f;
		}
	for(int i=0;i<1000;i++){
		Tree::distanceSwapPrecalc[i] = (float) pow(Tree::distanceSwapBias, i);
		if(Tree::distanceSwapPrecalc[i] != Tree::distanceSwapPrecalc[i]) Tree::distanceSwapPrecalc[i]=0.0f;	
		}

	Tree::meanBrlenMuts	= conf->meanBrlenMuts;
	Tree::alpha		= conf->gammaShapeBrlen;
	Tree::treeRejectionThreshold = conf->treeRejectionThreshold;
	Tree::min_brlen = conf->minBrlen;
	Tree::max_brlen = conf->maxBrlen;
	Tree::exp_starting_brlen = conf->startingBrlen;
#ifdef SINGLE_PRECISION_FLOATS
	for(int i=0;i<30;i++){
		precalcThresh[i] = exp((FLOAT_TYPE)(-precalcIncr[i] - 1));
		precalcMult[i] =  exp((FLOAT_TYPE)(precalcIncr[i]));
		}
#endif
	//deal with the outgroup specification, if there is one
	if(conf->outgroupString.length() > 0){
		//NxsString s(conf->outgroupString.c_str());
		const char *o = conf->outgroupString.c_str();
		vector<int> nums;
		int pos1=0, pos2;
		while(pos1 < conf->outgroupString.size()){
			pos2 = conf->outgroupString.find(" ", pos1+1);
			string tax = conf->outgroupString.substr(pos1, pos2);
			nums.push_back(atoi(tax.c_str()));
			pos1 = pos2;
			}

		if(outgroup) outgroup->ClearBipartition();
		else outgroup = new Bipartition();
		outgroup->BipartFromNodenums(nums);
		
		//if the outgroup consists of multiple taxa, make a constraint that the ingroup be monophyletic
/*		if(nums.size() > 1){
			Constraint out(outgroup, true);
			out.Standardize();
			out.SetAsOutgroup();

			int num=1;
			for(vector<Constraint>::iterator con=constraints.begin();con!=constraints.end();con++){
				if((*con).IsCompatibleWithConstraint(&out) == false) throw ErrorException("Specified outgroup is not compatible with constraint number %d!", num);
				num++;
				}
			constraints.push_back(out);
			}
*/		}
	}

//DJZ 4-28-04
//adding the ability to read in treestrings in which the internal node numbers are specified.  I'd like to make the
//internal numbers be specified the way that internal node labels are according to the newick format, ie directly after
//the closing paren that represents the internal node.  But, that makes going from string -> tree annoying
//because by the time the internal node number would be read the treeNode structure would have already been created.
//So, the internal node numbers will go just BEFORE the opening paren that represents that node
//Example:  50(1:.05, 2:.02):.1 signifies a node numbered 50 that is ancestral to 1 and 2.
Tree::Tree(const char* s, bool numericalTaxa , bool allowPolytomies /*=false*/){
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
				s++;
				temp->dlen = (FLOAT_TYPE) atof( info.c_str() );
				if(temp->dlen < min_brlen){
					outman.UserMessage("->Branch of length %s is less than min of %.1e.  Setting to min.", info.c_str(), min_brlen);
					temp->dlen = min_brlen;
					}
				else if (temp->dlen > max_brlen){
					outman.UserMessage("->Branch of length %f is greater than max of %.0f.  Setting to max.", temp->dlen, max_brlen);
					temp->dlen = max_brlen;
					}
				}
			else if(*s==','||*s==')'){
				temp->dlen=Tree::exp_starting_brlen;
#ifdef STOCHASTIC_STARTING_BLENS
				temp->dlen *= rnd.gamma(1.0);
#endif
				}
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
					numTipsAdded++;
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
						temp->dlen = (FLOAT_TYPE)atof( info.c_str() );
						if(temp->dlen < min_brlen){
							outman.UserMessage("->Branch of length %s is less than min of %.1e.  Setting to min.", info.c_str(), min_brlen);
							temp->dlen = min_brlen;
							}
						else if (temp->dlen > max_brlen){
							outman.UserMessage("->Branch of length %f is greater than max of %.0f.  Setting to max.", temp->dlen, max_brlen);
							temp->dlen = max_brlen;
							}
						}
					else{
						temp->dlen = Tree::exp_starting_brlen;
#ifdef STOCHASTIC_STARTING_BLENS
						temp->dlen *= rnd.gamma(1.0);
#endif
						}
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
	if(allowPolytomies == false) root->CheckforPolytomies();
	root->CheckTreeFormation();
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
		if(modSpec.IsCodon() || modSpec.IsCodonAminoAcid())
			allNodes[i]->tipData=(char *) static_cast<const CodonData*>(data)->GetCodonRow(i-1);
		else if(modSpec.IsAminoAcid())
			allNodes[i]->tipData=(char *)(data)->GetRow(i-1);
		else{
			allNodes[i]->tipData=data->GetAmbigString(i-1);
#ifdef OPEN_MP
			allNodes[i]->ambigMap=data->GetAmbigToCharMap(i-1);
#endif
			}
		}
	
	numTipsAdded=0;
	numNodesAdded=1;//root
	numTipsTotal=data->NTax();
	numNodesTotal=2*data->NTax()-2;
	lnL=0.0;

	calcs=0;
	numBranchesAdded=0;
	taxtags=new int[numTipsTotal+1];

#ifdef EQUIV_CALCS
	//need to do the root too, since that node is sometimes stolen
	allNodes[0]->tipData = new char[data->NChar()];
	for(int i=data->NTax()+1;i<numNodesTotal;i++){
		allNodes[i]->tipData = new char[data->NChar()];
		}
	dirtyEQ=true;
#endif
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
		if(modSpec.IsCodon() || modSpec.IsCodonAminoAcid())
			allNodes[i]->tipData=(char *) static_cast<const CodonData*>(data)->GetCodonRow(i-1);
		else if(modSpec.IsAminoAcid())
			allNodes[i]->tipData=(char *)(data)->GetRow(i-1);
		else{
			allNodes[i]->tipData=data->GetAmbigString(i-1);
#ifdef OPEN_MP
			allNodes[i]->ambigMap=data->GetAmbigToCharMap(i-1);
#endif
			}
		}
	
	numTipsAdded=0;
	numNodesAdded=1;//root
	numTipsTotal=data->NTax();
	numNodesTotal=2*data->NTax()-1;
	lnL=0.0;
	numBranchesAdded=0;
	taxtags=new int[numTipsTotal+1];
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

int Tree::BrlenMutate(){
	//random_binomial is now called with the mean number of blen muts, which is easiser to specify across datasets
	//than is a per branch probability
	int numBrlenMuts;
	if(rnd.uniform() < 0.05){//do a whole tree rescale occasionally
		ScaleWholeTree();
		numBrlenMuts = numNodesTotal - 1;
		}
	else{
		do{
			numBrlenMuts=rnd.random_binomial(numNodesTotal-1, meanBrlenMuts);
			}while(numBrlenMuts==0);
		for(int i=0;i<numBrlenMuts;i++){
			int branch=GetRandomNonRootNode();
			allNodes[branch]->dlen*=rnd.gamma( Tree::alpha );
			SweepDirtynessOverTree(allNodes[branch]);
			allNodes[branch]->dlen = (allNodes[branch]->dlen > min_brlen ? (allNodes[branch]->dlen < max_brlen ? allNodes[branch]->dlen : max_brlen) : min_brlen);
			}
		}
	return numBrlenMuts;
	}

void Tree::PerturbAllBranches(){
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		allNodes[i]->dlen*=rnd.gamma(100);
		}
	MakeAllNodesDirty();
	}

void Tree::ScaleWholeTree(FLOAT_TYPE factor/*=-1.0*/){
	if(factor==-1.0) factor = rnd.gamma( Tree::alpha );
	//9-12-06 Stupid!  Why the hell was this only scaling the internals?
	//for(int i=numTipsTotal;i<numNodesTotal;i++){
	for(int i=1;i<numNodesTotal;i++){
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
		numBrlenMuts=rnd.random_binomial((int)subtreeMemberNodes.size(), meanBrlenMuts);
		}while(numBrlenMuts==0);
	for(int i=0;i<numBrlenMuts;i++){
		int branch=subtreeMemberNodes[(int)(rnd.uniform()*subtreeMemberNodes.size())];//can't mutate the root
		allNodes[branch]->dlen*=rnd.gamma( Tree::alpha );
		SweepDirtynessOverTree(allNodes[branch]);
		allNodes[branch]->dlen = (allNodes[branch]->dlen > min_brlen ? (allNodes[branch]->dlen < max_brlen ? allNodes[branch]->dlen : max_brlen) : min_brlen);
		}	
	return numBrlenMuts;
	}

void Tree::MakeTrifurcatingRoot(bool reducenodes, bool clasAssigned ){
	//reducenodes should only =1 if this function is called after generating a random tree
	//or after reading in a tree with a bifurcating root.  DO NOT call with reducenodes=1 if
	//this is being used after one of the initial root branches was pruned off

	//clasAssigned should be true if the clas have been assigned to the nodes by the claManager.
	//(ie, not right after tree creation)
	TreeNode *t1, *t2, *removedNode;
	FLOAT_TYPE l;
	assert(root->left->next==root->right);
	if(root->left->IsInternal()){
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
	nd->dlen = Tree::exp_starting_brlen;
#ifdef STOCHASTIC_STARTING_BLENS
	nd->dlen *= rnd.gamma(1.0);
#endif
	if(nd->dlen < min_brlen) nd->dlen = min_brlen;
	else if(nd->dlen > max_brlen) nd->dlen = max_brlen;
	
	nd->next=nd->prev=NULL;//in case this node was connected in some other tree
	
	//Make sure that the root has 3 decendents
	if(numBranchesAdded<3)
		{root->AddDes(nd);
		}
	else
		{// If we're not adding directly to the root node, then we will need
		// a connector node and make the new terminal its left des
		TreeNode* connector=allNodes[placeInAllNodes++];
		numNodesAdded++;
		connector->dlen = Tree::exp_starting_brlen;
#ifdef STOCHASTIC_STARTING_BLENS
		connector->dlen *= rnd.gamma(1.0);
#endif
		nd->dlen = (nd->dlen > min_brlen ? nd->dlen : min_brlen);
		connector->left=connector->right=NULL;
		connector->AddDes(nd);
		
		//select a branch to break with the connector
		int k = rnd.random_int( numBranchesAdded ) + 1;
		TreeNode* otherDes = root->FindNode( k );

		// replace puts connection in the tree where otherDes had been
		otherDes->SubstituteNodeWithRespectToAnc(connector);
		
		//add otherDes back to the tree as the sister to the new tip
		connector->AddDes(otherDes);
		numBranchesAdded++;//numBranchesAdded needs to be incremented twice because a total of two branches have been added
		}
	numBranchesAdded++;
	numNodesAdded++;
	numTipsAdded++;
	}

void Tree::AddRandomNodeWithConstraints(int nodenum, int &placeInAllNodes, Bipartition &mask){
	//the trick here with the constraints is that only a subset of the taxa will be in the
	//growing tree.  To properly determine bipartition comptability a mask consisting of only
	//the present taxa with need to be used

	assert(nodenum>0 && nodenum<=numTipsTotal);  //should be adding a terminal
	TreeNode* nd=allNodes[nodenum];
	Bipartition temp;
	mask += temp.TerminalBipart(nodenum);
	nd->dlen = Tree::exp_starting_brlen;
#ifdef STOCHASTIC_STARTING_BLENS
	nd->dlen *= rnd.gamma(1.0);
#endif
	if(nd->dlen < min_brlen) nd->dlen = min_brlen;
	else if(nd->dlen > max_brlen) nd->dlen = max_brlen;
	
	nd->next=nd->prev=NULL;//in case this node was connected in some other tree
	
	//Make sure that the root has 3 decendents
	if(numBranchesAdded<3)
		{root->AddDes(nd);
		}
	else
		{// If we're not adding directly to the root node, then we will need
		// a connector node and make the new terminal its left des
		TreeNode* connector=allNodes[placeInAllNodes++];
		numNodesAdded++;
		connector->dlen = Tree::exp_starting_brlen;
#ifdef STOCHASTIC_STARTING_BLENS
		connector->dlen *= rnd.gamma(1.0);
#endif
		connector->dlen = (connector->dlen > min_brlen ? connector->dlen : min_brlen);
		connector->left=connector->right=NULL;
		connector->AddDes(nd);
		
		//select a branch to break with the connector
		int k;
		TreeNode *otherDes;
		bool compat;

		CalcBipartitions();
		nd->CalcBipartition();
		nd->bipart->Standardize();

		do{
			k = rnd.random_int( numBranchesAdded ) + 1;
			otherDes = root->FindNode( k );
			compat=true;
			for(vector<Constraint>::iterator conit=constraints.begin();conit!=constraints.end();conit++){
				//BACKBONE
				if((*conit).IsBackbone()){
					assert((*conit).IsPositive());
					//if the taxon being added isn't in the backbone, it can go anywhere
					if((*conit).GetBackboneMask()->ContainsTaxon(nd->nodeNum)){
						compat=AllowedByPositiveBackboneConstraintWithMask(&(*conit), &mask, nd, otherDes);
						//need to do this because we screw with the biparts while checking for their compatability
						//with a backbone constraint
						CalcBipartitions();
						}
					if(compat == false) break;
					}
				else{
					if((*conit).IsPositive())
						compat=AllowedByPositiveConstraintWithMask(&(*conit), &mask, nd, otherDes);
					else
						compat=AllowedByNegativeConstraintWithMask(&(*conit), &mask, nd, otherDes);
					if(compat==false) break;
					}
				}
			}while(compat == false);

		// replace puts connection in the tree where otherDes had been
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
	
void Tree::RecombineWith( Tree *t, bool sameModel, FLOAT_TYPE optPrecision ){
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
	connector->dlen=broken->dlen*ZERO_POINT_FIVE;
	broken->dlen-=connector->dlen;

	TraceDirtynessToRoot(connector);
	OptimizeBranchesAroundNode(connector, optPrecision, 0);
	}

TreeNode *Tree::ContainsBipartition(const Bipartition *bip){
	//note that this doesn't work for terminals (but there's no reason to call for them anyway)
	//find a taxon that appears "on" in the bipartition

	//turning this back on
	int tax=bip->FirstPresentTaxon();
	
	//now start moving down the tree from taxon 1 until a bipart that
	//conflicts or a match is found
	//TreeNode *nd=allNodes[1]->anc;
	TreeNode *nd=allNodes[tax]->anc;
	while(nd->anc){
		if(nd->bipart->IsASubsetOf(bip) == false) return NULL;
		else if(nd->bipart->EqualsEquals(bip)) return nd;
		else nd=nd->anc;
		}
	return NULL;
	}

TreeNode *Tree::ContainsBipartitionOrComplement(const Bipartition *bip){
	//this version will detect if the same bipartition exists in the trees, even
	//if it is in different orientation, which could happen due to rooting
	//differences

	//find a taxon that appears "on" in the bipartition
	int tax=bip->FirstPresentTaxon();
	
	//now start moving down the tree from that taxon until a bipart that
	//conflicts or a match is found
	//7/17/07 changing this to start from the trivial terminal branch, rather
	//then its anc
	TreeNode *nd=allNodes[tax];
	while(nd->anc){
		if(nd->bipart->IsASubsetOf(bip) == false) break;
		else if(nd->bipart->EqualsEquals(bip)) return nd;
		else nd=nd->anc;
		}
		
	//find a taxon that is NOT "on" in the bipartition
	tax=bip->FirstNonPresentTaxon();
	
	//now start moving down the tree from that taxon until a bipart that
	//conflicts or a match is found
	//7/17/07 changing this to start from the trivial terminal branch, rather
	//then its anc
	nd=allNodes[tax];
	while(nd->anc){
		//if(nd->bipart->ComplementIsASubsetOf(bip) == false){
		if(bip->IsASubsetOf(nd->bipart) == false){
			return NULL;
			}
		else if(nd->bipart->EqualsEquals(bip)) return nd;
		else nd=nd->anc;
		}
		
	return NULL;
	}

TreeNode *Tree::ContainsMaskedBipartitionOrComplement(const Bipartition *bip, const Bipartition *mask){
	//as in ContainsMaskedBipartitionOrComplement, but bits not on in the
	//mask are ignored

	//find a taxon that appears "on" in the bipartition and is on in the mask
	Bipartition temp = *bip;
	temp.AndEquals(mask);
	int tax=temp.FirstPresentTaxon();
	
	//now start moving down the tree from that taxon until we find a
	//match or reach the root
	TreeNode *nd=allNodes[tax]->anc;
	temp = *bip;
	temp.Complement();
	while(nd->anc){
		if(nd->bipart->MaskedEqualsEquals(bip, mask)) return nd;
		if(nd->bipart->MaskedEqualsEquals(&temp, mask)) return nd;
		else nd=nd->anc;
		}

	//find a taxon that is NOT "on" in the bipartition
	temp = *bip;
	temp.Complement();
	temp.AndEquals(mask);
	tax=temp.FirstPresentTaxon();
	
	//now start moving down the tree from that taxon until we find a
	//match or reach the root
	nd=allNodes[tax]->anc;
	temp = *bip;
	temp.Complement();
	while(nd->anc){
		if(nd->bipart->MaskedEqualsEquals(bip, mask)) return nd;
		if(nd->bipart->MaskedEqualsEquals(&temp, mask)) return nd;
		else nd=nd->anc;
		}
		
	return NULL;
	}

int Tree::SubtreeBasedRecombination( Tree *t, int recomNodeNum, bool sameModel, FLOAT_TYPE optPrecision){
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
	assert(other->IsNotRoot());
	bool identical;
	
	if(other->IsRoot() == false){
		if(other->IsTerminal()) return true;
		identical=(ContainsBipartition(other->bipart) != NULL);
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
	//overall topology, but assumes the same rooting
	bool identical;
	
	if(other->IsRoot() == false){
		if(other->IsTerminal()) return true;
		identical= (ContainsBipartition(other->bipart) != NULL);
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

bool Tree::IdenticalTopologyAllowingRerooting(const TreeNode *other){
	//this is intitially called with the root, it will detect any difference in the 
	//overall topology
	bool identical;
	
	if(other->IsRoot() == false){
		if(other->IsTerminal()) return true;
		identical= (ContainsBipartitionOrComplement(other->bipart) != NULL);
		if(identical==true){
			identical=IdenticalTopologyAllowingRerooting(other->left);
			if(identical==true)
				identical=IdenticalTopologyAllowingRerooting(other->right);
			}
		}
	else{
		TreeNode *nd=other->left;
		while(nd != NULL){
			identical=IdenticalTopologyAllowingRerooting(nd);
			if(identical == false){
				return identical;
				}
			nd=nd->next;
			}
		}
	return identical;
	}

int Tree::BipartitionBasedRecombination( Tree *t, bool sameModel, FLOAT_TYPE optPrecision){
	//find a bipartition that is shared between the trees
	TreeNode *tonode, *fromnode;
	bool found=false;
	int tries=0;
	while(!found && (++tries<50)){
		int i;
		do{
			i=GetRandomInternalNode();
			//WTF!!!  How did this work?
			}while((allNodes[i]->left->IsTerminal() && allNodes[i]->right->IsTerminal()));
			//}while((t->allNodes[i]->left->IsTerminal() && t->allNodes[i]->right->IsTerminal()));
		//fromnode=t->ContainsBipartition(allNodes[i]->bipart);
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
			FLOAT_TYPE toscore, fromscore;
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
		//OptimizeBranchLength(optPrecision, tonode, true);
		OptimizeBranchesWithinRadius(tonode, optPrecision, 0, NULL);

		Score(tonode->nodeNum);
		}
	else return -1;
	return 1;
	}
	
//this is essentially a version of TopologyMutator that goes through cut nodes in order
//and for each cut node goes through the broken nodes in order.  It the swaps are performed
//on a temporary tree
void Tree::DeterministicSwapperByCut(Individual *source, double optPrecision, int range, bool furthestFirst){

	TreeNode *cut;
	int swapNum=0;
	double bestScore = -DBL_MAX;
	
	Individual tempIndiv;
	tempIndiv.treeStruct=new Tree();
	
	tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, source);

	//ensure that the starting tree is optimal up to the required precision
	FLOAT_TYPE imp = 999.9;
	do{
		imp = tempIndiv.treeStruct->OptimizeAllBranches(optPrecision);
		}while(imp > 0.0);

	outman.UserMessage("starting score:%f", tempIndiv.treeStruct->lnL);

	char str[50];
	if(furthestFirst) sprintf(str, "determImpsCutR.%d.%f.tre", range, optPrecision);
	else sprintf(str, "determImpsCut.%d.%f.tre", range, optPrecision);
	ofstream better(str);
	better.precision(9);
	data->BeginNexusTreesBlock(better);

	if(furthestFirst) sprintf(str, "determCutR%d.%f.log", range, optPrecision);
	else sprintf(str, "determCut%d.%f.log", range, optPrecision);
	FILE *log = fopen(str, "w");
	
	//allocate a treeString
	double taxsize=log10((double) ((double)data->NTax())*data->NTax()*2);
	int stringSize=(int)((data->NTax()*2)*(10+DEF_PRECISION));
	char *treeString=new char[stringSize];
	stringSize--;
	treeString[stringSize]='\0';
	bool newBest=false;
	attemptedSwaps.ClearAttemptedSwaps();
	int startC, c=1;

	/*
	AttemptedSwapList availableSwaps;

	Bipartition proposed;
	for(int i=1;i<numNodesTotal;i++){
		cut=allNodes[i];
		GatherValidReconnectionNodes(range, cut, NULL);
		for(list<ReconNode>::iterator b = sprRang.begin();b != sprRang.end();b++){
			ReconNode *broken = &(*b);
			proposed.FillWithXORComplement(cut->bipart, allNodes[broken->nodeNum]->bipart);
			availableSwaps.AddSwap(proposed, cut->nodeNum, broken->nodeNum, broken->reconDist);
			}
		}
*/	

	int acceptedSwaps = 0;
	startC = c;

	while(1){
		cut=allNodes[c];
		GatherValidReconnectionNodes(range, cut, NULL);
		if(furthestFirst) sprRang.Reverse();

		for(list<ReconNode>::iterator b = sprRang.begin();b != sprRang.end();b++){
			ReconNode *broken = &(*b);
			
			//log the swap about to be performed.  Although this func goes through the swaps in order,
			//there will be duplication because of the way that NNIs are performed.  Two different cut
			//nodes can be reconnected with an NNI such that the same topology results
			bool unique=false;
			Bipartition proposed;
			proposed.FillWithXORComplement(cut->bipart, allNodes[broken->nodeNum]->bipart);
			unique = attemptedSwaps.AddSwap(proposed, cut->nodeNum, broken->nodeNum, broken->reconDist);

			if(unique){
				swapNum++;
				if(swapNum %100 == 0) fprintf(log, "%d\t%d\t%f\n", swapNum, acceptedSwaps, lnL);
				if(broken->withinCutSubtree == true){
					tempIndiv.treeStruct->ReorientSubtreeSPRMutate(cut->nodeNum, broken, optPrecision);
					}
				else{
					tempIndiv.treeStruct->SPRMutate(cut->nodeNum, broken, optPrecision, 0);
					}

				if(tempIndiv.treeStruct->lnL > (lnL+optPrecision)){

					outman.UserMessage("%f\t%f\t%d\t%d", tempIndiv.treeStruct->lnL, lnL - tempIndiv.treeStruct->lnL, c, b->reconDist);
					source->CopySecByRearrangingNodesOfFirst(source->treeStruct, &tempIndiv, true);
					lnL = tempIndiv.treeStruct->lnL;

					tempIndiv.treeStruct->root->MakeNewick(treeString, false, true);
					better << "tree " << c << "." << b->reconDist << "= [&U][" << lnL << "]" << treeString << ";" << endl;
					newBest = true;
					acceptedSwaps++;
					attemptedSwaps.ClearAttemptedSwaps();
					break;
					}
				else{
					tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, source, true);
					}
				}
			}
		c++;
		if(c == numNodesTotal) c = 1;
		if(newBest == true){
			startC = c;
			newBest = false;
			}
		else if(c == startC){
			outman.UserMessage("done. %d swaps, %d accepted", swapNum, acceptedSwaps);
			break;
			}
		}

	better << "end;";
	better.close();
	delete []treeString;
	fclose(log);

	tempIndiv.treeStruct->RemoveTreeFromAllClas();
	delete tempIndiv.treeStruct;
	tempIndiv.treeStruct=NULL;
	}

//this is essentially a version of TopologyMutator that goes through cut nodes in order
//and for each cut node goes through the broken nodes in order.  It the swaps are performed
//on a temporary tree
void Tree::DeterministicSwapperByDist(Individual *source, double optPrecision, int range, bool furthestFirst){

	TreeNode *cut;
	int swapNum=0;
	
	Individual tempIndiv;
	tempIndiv.treeStruct=new Tree();
	
	tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, source);

	//ensure that the starting tree is optimal up to the required precision
	FLOAT_TYPE imp = 999.9;
	do{
		imp = tempIndiv.treeStruct->OptimizeAllBranches(optPrecision);
		}while(imp > 0.0);

	outman.UserMessage("starting score:%f", tempIndiv.treeStruct->lnL);

	char str[50];
	if(furthestFirst) sprintf(str, "determImpsDistR.%d.%f.tre", range, optPrecision);
	else sprintf(str, "determImpsDist.%d.%f.tre", range, optPrecision);

	ofstream better(str);
	better.precision(9);
	data->BeginNexusTreesBlock(better);

	if(furthestFirst) sprintf(str, "determDistR%d.%f.log", range, optPrecision);
	else sprintf(str, "determDist%d.%f.log", range, optPrecision);
	FILE *log = fopen(str, "w");
	
	//allocate a treeString
	double taxsize=log10((double) ((double)data->NTax())*data->NTax()*2);
	int stringSize=(int)((data->NTax()*2)*(10+DEF_PRECISION));
	char *treeString=new char[stringSize];
	stringSize--;
	treeString[stringSize]='\0';
	bool newBest=false;
	attemptedSwaps.ClearAttemptedSwaps();
	int startC, c=1;

	/*
	AttemptedSwapList availableSwaps;

	Bipartition proposed;
	for(int i=1;i<numNodesTotal;i++){
		cut=allNodes[i];
		GatherValidReconnectionNodes(range, cut, NULL);
		for(list<ReconNode>::iterator b = sprRang.begin();b != sprRang.end();b++){
			ReconNode *broken = &(*b);
			proposed.FillWithXORComplement(cut->bipart, allNodes[broken->nodeNum]->bipart);
			availableSwaps.AddSwap(proposed, cut->nodeNum, broken->nodeNum, broken->reconDist);
			}
		}
*/	

	int currentDist;
	if(furthestFirst) currentDist = range;
	else currentDist = 1;
	int acceptedSwaps = 0;
	startC = c;
	do{
		cut=allNodes[c];
		//outman.UserMessageNoCR("cut=%d ", c);
		GatherValidReconnectionNodes(range, cut, NULL);

		for(list<ReconNode>::iterator b = sprRang.GetFirstNodeAtDist(currentDist);b != sprRang.end() && b->reconDist == currentDist;b++){
			ReconNode *broken = &(*b);
			
			//log the swap about to be performed.  Although this func goes through the swaps in order,
			//there will be duplication because of the way that NNIs are performed.  Two different cut
			//nodes can be reconnected with an NNI such that the same topology results
			bool unique=false;
			Bipartition proposed;
			proposed.FillWithXORComplement(cut->bipart, allNodes[broken->nodeNum]->bipart);
			unique = attemptedSwaps.AddSwap(proposed, cut->nodeNum, broken->nodeNum, broken->reconDist);

			if(unique){
				swapNum++;
				if(broken->withinCutSubtree == true){
					tempIndiv.treeStruct->ReorientSubtreeSPRMutate(cut->nodeNum, broken, optPrecision);
					}
				else{
					tempIndiv.treeStruct->SPRMutate(cut->nodeNum, broken, optPrecision, 0);
					}

				if(tempIndiv.treeStruct->lnL > (lnL+optPrecision)){
					outman.UserMessage("%f\t%f\t%d\t%d", tempIndiv.treeStruct->lnL, lnL - tempIndiv.treeStruct->lnL, c, b->reconDist);
				
					source->CopySecByRearrangingNodesOfFirst(source->treeStruct, &tempIndiv, true);
					lnL = tempIndiv.treeStruct->lnL;

					tempIndiv.treeStruct->root->MakeNewick(treeString, false, true);
					better << "tree " << c << "." << b->reconDist << "= [&U][" << lnL << "]" << treeString << ";" << endl;
					newBest = true;
					acceptedSwaps++;
					attemptedSwaps.ClearAttemptedSwaps();
					break;
					}
				else{
					tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, source, true);
					}
				if(swapNum %100 == 0) fprintf(log, "%d\t%d\t%f\n", swapNum, acceptedSwaps, lnL);
				}
			}
		c++;
		if(c == numNodesTotal) c = 1;
		if(newBest == true){
			startC = c;
			if(furthestFirst)
				currentDist = range;
			else
				currentDist = 1;
			newBest = false;
			}	
		else if(c == startC){
			if(furthestFirst) 
				currentDist--;
			else currentDist++;
			outman.UserMessage("dist = %d", currentDist);
			}
		}while(currentDist <= range && currentDist > 0);

	outman.UserMessage("done. %d swaps, %d accepted", swapNum, acceptedSwaps);

	better << "end;";
	better.close();
	delete []treeString;
	fclose(log);

	tempIndiv.treeStruct->RemoveTreeFromAllClas();
	delete tempIndiv.treeStruct;
	tempIndiv.treeStruct=NULL;
	}


//this is essentially a version of TopologyMutator that goes through cut nodes in order
//and for each cut node goes through the broken nodes in order.  It the swaps are performed
//on a temporary tree
void Tree::DeterministicSwapperRandom(Individual *source, double optPrecision, int range){

	TreeNode *cut;
	int swapNum=0;
	double bestScore = -DBL_MAX;
	
	Individual tempIndiv;
	tempIndiv.treeStruct=new Tree();
	
	tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, source);

	//ensure that the starting tree is optimal up to the required precision
	FLOAT_TYPE imp = 999.9;
	do{
		imp = tempIndiv.treeStruct->OptimizeAllBranches(optPrecision);
		}while(imp > 0.0);

	outman.UserMessage("starting score:%f", tempIndiv.treeStruct->lnL);

	char str[50];
	sprintf(str, "determImpsRand.%d.%f.tre", range, optPrecision);
	ofstream better(str);
	better.precision(9);
	data->BeginNexusTreesBlock(better);
	
	sprintf(str, "determRand%d.%f.log", range, optPrecision);
	FILE *log = fopen(str, "w");

	//allocate a treeString
	double taxsize=log10((double) ((double)data->NTax())*data->NTax()*2);
	int stringSize=(int)((data->NTax()*2)*(10+DEF_PRECISION));
	char *treeString=new char[stringSize];
	stringSize--;
	treeString[stringSize]='\0';
	bool newBest=false;
	attemptedSwaps.ClearAttemptedSwaps();
	int c=1;

	/*
	AttemptedSwapList availableSwaps;

	Bipartition proposed;
	for(int i=1;i<numNodesTotal;i++){
		cut=allNodes[i];
		GatherValidReconnectionNodes(range, cut, NULL);
		for(list<ReconNode>::iterator b = sprRang.begin();b != sprRang.end();b++){
			ReconNode *broken = &(*b);
			proposed.FillWithXORComplement(cut->bipart, allNodes[broken->nodeNum]->bipart);
			availableSwaps.AddSwap(proposed, cut->nodeNum, broken->nodeNum, broken->reconDist);
			}
		}
*/	

	int *completed = new int[numNodesTotal];
	for(int i=1;i<numNodesTotal;i++) completed[i] = 0;

	int acceptedSwaps = 0;

	while(1){
		int attempts = 0;
		do{
			c = GetRandomNonRootNode();
			if(attempts++ > numNodesTotal){
				int n=1;
				while(completed[n] && n < numNodesTotal) n++;
				if(n == numNodesTotal){
					outman.UserMessage("done. %d swaps, %d accepted", swapNum, acceptedSwaps);

					better << "end;";
					better.close();
					delete []treeString;
					fclose(log);

					tempIndiv.treeStruct->RemoveTreeFromAllClas();
					delete tempIndiv.treeStruct;
					tempIndiv.treeStruct=NULL;
					return;
					}
				}
			}while(completed[c]);
		cut=allNodes[c];
		//outman.UserMessageNoCR("cut=%d ", c);
		GatherValidReconnectionNodes(range, cut, NULL);

		//for(list<ReconNode>::iterator b = sprRang.GetFirstNodeAtDist(currentDist);b != sprRang.end() && b->reconDist == currentDist;b++){
		bool noSwapFound = true;
		listIt b;
		while(sprRang.size() > 0){
			b = sprRang.NthElement(rnd.random_int(sprRang.size()));
			ReconNode *broken = &(*b);
			
			//log the swap about to be performed.  Although this func goes through the swaps in order,
			//there will be duplication because of the way that NNIs are performed.  Two different cut
			//nodes can be reconnected with an NNI such that the same topology results
			bool unique=false;
			Bipartition proposed;
			proposed.FillWithXORComplement(cut->bipart, allNodes[broken->nodeNum]->bipart);
			unique = attemptedSwaps.AddSwap(proposed, cut->nodeNum, broken->nodeNum, broken->reconDist);

			if(unique){
				swapNum++;
				if(broken->withinCutSubtree == true){
					tempIndiv.treeStruct->ReorientSubtreeSPRMutate(cut->nodeNum, broken, optPrecision);
					}
				else{
					tempIndiv.treeStruct->SPRMutate(cut->nodeNum, broken, optPrecision, 0);
					}

				if(tempIndiv.treeStruct->lnL > (lnL+optPrecision)){
					outman.UserMessage("%f\t%f\t%d\t%d", tempIndiv.treeStruct->lnL, lnL - tempIndiv.treeStruct->lnL, c, b->reconDist);
					source->CopySecByRearrangingNodesOfFirst(source->treeStruct, &tempIndiv, true);
					lnL = tempIndiv.treeStruct->lnL;

					tempIndiv.treeStruct->root->MakeNewick(treeString, false, true);
					better << "tree " << c << "." << b->reconDist << "= [&U][" << lnL << "]" << treeString << ";" << endl;
					newBest = true;
					acceptedSwaps++;
					attemptedSwaps.ClearAttemptedSwaps();
					for(int i=0;i<numNodesTotal;i++) completed[i] = 0;
					break;
					}
				else{
					tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, source, true);
					break;
					}
				}
			else sprRang.RemoveElement(b);
			}
		if(swapNum %100 == 0) fprintf(log, "%d\t%d\t%f\n", swapNum, acceptedSwaps, lnL);
		if(sprRang.size() == 0){
			outman.UserMessage("completed %d", c);
			completed[c] = 1;
			}
		}
	}

//this function now returns the reconnection distance, with it being negative if its a
//subtree reorientation swap
int Tree::TopologyMutator(FLOAT_TYPE optPrecision, int range, int subtreeNode){
	//All topology mutations go through here now.  Range will be 1 in the case of NNI's
	//Range will be some small number in the case of limSPR's and will be 999999 in the case
	//of random SPR's
	TreeNode *cut;
	ReconNode *broken;
	bool unique;

#ifdef EQUIV_CALCS
	dirtyEQ = true;
#endif
	
	int err=0;
	int ret=0;
	int tryNum = 0;
	do{
		do{
			cut=allNodes[GetRandomNonRootNode()];
			GatherValidReconnectionNodes(range, cut, NULL);
			}while(sprRang.size()==0);

		if((uniqueSwapBias == 1.0 && distanceSwapBias == 1.0) || range < 0)
			broken = sprRang.RandomReconNode();
		else{//only doing this on limSPR and NNI
			err = AssignWeightsToSwaps(cut);
			err = err && (tryNum++ < 5);
#ifdef SWAP_BASED_TERMINATION
			if(!err){
#else
			if(1){
#endif
				sprRang.CalcProbsFromWeights();
				broken = sprRang.ChooseNodeByWeight();
				}
			}

#ifdef SWAP_BASED_TERMINATION
		if(!err){
#else
		if(1){
#endif
			//log the swap about to be performed
			if( ! ((uniqueSwapBias == 1.0 && distanceSwapBias == 1.0) || range < 0)){
				Bipartition proposed;
				proposed.FillWithXORComplement(cut->bipart, allNodes[broken->nodeNum]->bipart);
				unique = attemptedSwaps.AddSwap(proposed, cut->nodeNum, broken->nodeNum, broken->reconDist);
				uniqueSwapTried = uniqueSwapTried || unique;
				//uniqueSwapTried = uniqueSwapTried || attemptedSwaps.AddSwap(proposed, cut->nodeNum, broken->nodeNum, broken->reconDist);
				}

			if(broken->withinCutSubtree == true){
				#ifdef OPT_DEBUG
					optsum << "reorientSPR\t" << broken->reconDist << "\t" << range << "\n";
				#endif
				#ifdef VARIABLE_OPTIMIZATION
					if(unique == true) ReorientSubtreeSPRMutateDummy(cut->nodeNum, broken, optPrecision);
					else return broken->reconDist * -1;
				#endif
				ReorientSubtreeSPRMutate(cut->nodeNum, broken, optPrecision);	
				ret=broken->reconDist * -1;
				}
			else{
				#ifdef OPT_DEBUG
					optsum << "SPR\t" << broken->reconDist << "\t" << range << "\n";
				#endif
				#ifdef VARIABLE_OPTIMIZATION
					if(unique == true) err=SPRMutateDummy(cut->nodeNum, broken, optPrecision, subtreeNode);
					else return broken->reconDist;
				#endif
				err=SPRMutate(cut->nodeNum, broken, optPrecision, subtreeNode);
				ret=broken->reconDist;
				}
			#ifdef OUTPUT_UNIQUE_TREES
			if(unique == true){
				output_tree = true;
				//uni.precision(9);
				if(broken->withinCutSubtree == false) uni << "SPR" << "\t" << broken->reconDist << "\t" << lnL << "\t" << cut->nodeNum << "\t" << broken->nodeNum << "\n";
				else uni << "reSPR" << "\t" << broken->reconDist << "\t" << lnL << "\t" << cut->nodeNum << "\t" << broken->nodeNum << "\n";
				}
			#endif
			}
		}while(err);

#ifndef NDEBUG
	CalcBipartitions();

	for(vector<Constraint>::iterator conit=constraints.begin();conit!=constraints.end();conit++){
		if((*conit).IsBackbone()){
			assert((*conit).IsPositive());
			assert(ContainsMaskedBipartitionOrComplement((*conit).GetBipartition(), (*conit).GetBackboneMask()) != NULL);
			}
		else{
			if((*conit).IsPositive())
				assert(ContainsBipartitionOrComplement((*conit).GetBipartition()) != NULL);
			else assert(ContainsBipartitionOrComplement((*conit).GetBipartition()) == NULL);
			}
		}
#endif
	return ret;
	}

void Tree::GatherValidReconnectionNodes(int maxDist, TreeNode *cut, const TreeNode *subtreeNode){
	/* 7/11/06 making this function more multipurpose
	It now assumes that the cut branch has NOT YET BEEN DETACHED. This is important so that
	when branches are chosen without a viable reconnection due to a constraint another cut
	can be chosen without having the put the tree back together again
	1.	Gather all nodes within maxRange.  This can include nodes that are des of the 
		cut node.  In this case the portion of the tree containing the root is considered
		the subtree to be reattached, and the swap would be done by ReorientSubtreeSPRMutate
	2.	Keep information on the potential reconnection nodes, including reconnection distance and
		branchlength distance.  This allows for various schemes of differentially weighting the 
		swaps.
	3.	filter out reconnection nodes incompatible with constraints
	*/
	sprRang.clear();
	const TreeNode *center=cut->anc;
	
	//add the descendent branches
	if(center->left != cut) 
		sprRang.AddNode(center->left->nodeNum, 0, (float) center->left->dlen);
	if(center->left->next != cut) 
		sprRang.AddNode(center->left->next-> nodeNum, 0, (float) center->left->next->dlen);
	
	//add either the center node itself or the third descendent in the case of the root
	if(center->IsNotRoot()){
		if(center->anc != subtreeNode)
			sprRang.AddNode(center->nodeNum, 0, (float) center->dlen);
		}
	else{
		if(center->left->next->next != cut)
			sprRang.AddNode(center->left->next->next->nodeNum, 0, (float) center->left->next->next->dlen);
		}
	
	assert(sprRang.size() == 2);
	
	for(int curDist = 0; curDist < maxDist || maxDist < 0; curDist++){
		list<ReconNode>::iterator it=sprRang.GetFirstNodeAtDist(curDist);
		if(it == sprRang.end()){
			break; //need this to break out of loop when curDist exceeds any branches in the tree
			}
		for(; it != sprRang.end() && it->reconDist == curDist; it++){
			TreeNode *cur=allNodes[it->nodeNum];
			assert(cur->IsNotRoot());
			
			if(cur->left!=NULL && cur->left!=cut) 
			    sprRang.AddNode(cur->left->nodeNum, curDist+1, (float) (it->pathlength + cur->left->dlen));
			if(cur->right!=NULL && cur->right!=cut) 
		    	sprRang.AddNode(cur->right->nodeNum, curDist+1, (float) (it->pathlength + cur->right->dlen));
			if(cur->next!=NULL && cur->next!=cut){
			    sprRang.AddNode(cur->next->nodeNum, curDist+1, (float) (it->pathlength + cur->next->dlen));
			    if(cur->next->next!=NULL && cur->next->next!=cut){//if cur is the left descendent of the root
			    	sprRang.AddNode(cur->next->next->nodeNum, curDist+1, (float) (it->pathlength + cur->next->next->dlen));
			    	}
			    }
			if(cur->prev!=NULL && cur->prev!=cut){
			    sprRang.AddNode(cur->prev->nodeNum, curDist+1, (float) (it->pathlength + cur->prev->dlen));
			    if(cur->prev->prev!=NULL && cur->prev->prev!=cut){//if cur is the right descendent of the root
			    	sprRang.AddNode(cur->prev->prev->nodeNum, curDist+1, (float) (it->pathlength + cur->prev->prev->dlen));
			    	}
			    }
		    if(cur->anc->nodeNum != 0){//if the anc is not the root, add it.
		    	if(cur->anc!=subtreeNode){
			    	sprRang.AddNode(cur->anc->nodeNum, curDist+1, (float) (it->pathlength + cur->anc->dlen));
			 		}
			 	}
		    }
		}
	
	if(maxDist != 1 && cut->IsInternal()){
		//Gather nodes within the cut subtree to allow SPRs in which the portion of the tree containing
		//the root is considered the subtree to be reattached
		//start by adding cut's left and right
		sprRang.AddNode(cut->left->nodeNum, 0, (float) cut->left->dlen, true);
		sprRang.AddNode(cut->right->nodeNum, 0, (float) cut->right->dlen, true);

		for(int curDist = 0; curDist < maxDist || maxDist < 0; curDist++){
			list<ReconNode>::iterator it=sprRang.GetFirstNodeAtDistWithinCutSubtree(curDist);	
			if(it == sprRang.end()){
				break; //need this to break out of loop when curDist exceeds any branches in the tree
				}
			for(; it != sprRang.end() && it->reconDist == curDist; it++){
				TreeNode *cur=allNodes[it->nodeNum];
				
				if(cur->left!=NULL) 
					sprRang.AddNode(cur->left->nodeNum, curDist+1, (float) (it->pathlength + cur->left->dlen), true);
				if(cur->right!=NULL) 
		    		sprRang.AddNode(cur->right->nodeNum, curDist+1, (float) (it->pathlength + cur->right->dlen), true);
				if(cur->next!=NULL){
					sprRang.AddNode(cur->next->nodeNum, curDist+1, (float) (it->pathlength + cur->next->dlen), true);
					}
				}
			}
		}

    //remove general unwanted nodes from the subset
	sprRang.RemoveNodesOfDist(0); //remove branches adjacent to cut
//	if(maxDist != 1)
//		sprRang.RemoveNodesOfDist(1); //remove branches equivalent to NNIs
	
	CalcBipartitions();

#ifdef CONSTRAINTS
	//now deal with constraints, if any
	if(constraints.size() > 0){
		Bipartition scratch;

		for(vector<Constraint>::iterator conit=constraints.begin();conit!=constraints.end();conit++){
			if(sprRang.size() != 0){
				listIt it=sprRang.begin();
				do{
					TreeNode* broken=allNodes[it->nodeNum];
					//if(AllowedByConstraint(&(*conit), cut, broken, scratch) == false) it=sprRang.RemoveElement(it);
					if(AllowedByConstraint(&(*conit), cut, &*it, scratch) == false) it=sprRang.RemoveElement(it);
					else it++;
					}while(it != sprRang.end());
				}
			else return;
			}
		}
#endif
	}


bool Tree::AssignWeightsToSwaps(TreeNode *cut){
	//Assign weights to each swap (reconnection node) based on 
	//some criterion

	Bipartition proposed;
	list<Swap>::iterator thisSwap;
	bool someUnique = false;

	for(listIt it = sprRang.begin();it != sprRang.end();it++){
		bool found;
		proposed.FillWithXORComplement(cut->bipart, allNodes[(*it).nodeNum]->bipart);
		Swap tmp=Swap(proposed, cut->nodeNum, (*it).nodeNum, (*it).reconDist);
		thisSwap = attemptedSwaps.FindSwap(tmp, found);

		if(found == false){
			someUnique = true;
			if((*it).reconDist - 1 < 1000)
				(*it).weight = distanceSwapPrecalc[(*it).reconDist - 1];
			else (*it).weight = 0.0;
			}
		else{
			if((*it).reconDist - 1 < 1000 && (*thisSwap).Count() < 500)
				(*it).weight = uniqueSwapPrecalc[(*thisSwap).Count()] * distanceSwapPrecalc[(*it).reconDist - 1];
			else (*it).weight = 0.0;
			}
		}
	return someUnique==false;
	}

int Tree::SPRMutateDummy(int cutnum, ReconNode *broke, FLOAT_TYPE optPrecision, int subtreeNode){
	//this is just a spoof version of SPRMutate that will perform the same mutation
	//several times with different optimiation settings, but will otherwise 
	//maintain exactly the same program flow because it resets the seed

#ifndef VARIABLE_OPTIMIZATION
	assert(0);
#else
	Individual tempIndiv;
	tempIndiv.treeStruct=new Tree();
	
	Individual sourceIndiv;
	sourceIndiv.treeStruct=this;
	sourceIndiv.mod->CopyModel(this->mod);
		
	int savedSeed;
	
	var.precision(10);
	var << "SPR" << "\t" << broke->reconDist << "\t" << lnL << "\t";
	
	tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, &sourceIndiv);
//	FLOAT_TYPE prec[5]={(FLOAT_TYPE).01, (FLOAT_TYPE).5, (FLOAT_TYPE).01, (FLOAT_TYPE).01, (FLOAT_TYPE).01};
	FLOAT_TYPE origThresh = treeRejectionThreshold;
/*	
	treeRejectionThreshold = 10000;
	for(int i=0;i<1;i++){
		savedSeed = rnd.seed();
		optCalcs = 0;
		tempIndiv.treeStruct->SPRMutate(cutnum, broke, optPrecision, 0);
		var << tempIndiv.treeStruct->lnL << "\t" << optCalcs << "\t";
		optCalcs = 0;
		rnd.set_seed(savedSeed);
		tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, &sourceIndiv, true);
		}
	treeRejectionThreshold = -10000;
	for(int i=0;i<1;i++){
		savedSeed = rnd.seed();
		optCalcs = 0;
		tempIndiv.treeStruct->SPRMutate(cutnum, broke, optPrecision, 0);
		var << tempIndiv.treeStruct->lnL << "\t" << optCalcs << "\t";
		optCalcs = 0;
		rnd.set_seed(savedSeed);
		tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, &sourceIndiv, true);
		}
*/
/*	for(int i=0;i<1;i++){
		treeRejectionThreshold = origThresh;
		savedSeed = rnd.seed();
		optCalcs = 0;
		tempIndiv.treeStruct->SPRMutate(cutnum, broke, optPrecision, 0);
		var << tempIndiv.treeStruct->lnL << "\t" << optCalcs << "\t";
		optCalcs = 0;
		rnd.set_seed(savedSeed);
		tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, &sourceIndiv, true);
		}
*/
	treeRejectionThreshold = origThresh;
	tempIndiv.treeStruct->RemoveTreeFromAllClas();
	delete tempIndiv.treeStruct;
	tempIndiv.treeStruct=NULL;
	sourceIndiv.treeStruct=NULL;
	optCalcs = 0;
	SPRMutate(cutnum, broke, optPrecision, 0);
	var << lnL << "\t" << optCalcs << "\n";
	optCalcs = 0;
#endif
	return 1;
	}


// 7/21/06 This function is now called by TopologyMutator to actually do the rearrangement
//It has the cut and broken nodenums passed in.  It also does NNI's
int Tree::SPRMutate(int cutnum, ReconNode *broke, FLOAT_TYPE optPrecision, int subtreeNode){
	//if the optPrecision passed in is < 0 it means that we're just trying to 
	//make the tree structure for some reason, but don't have CLAs allocated
	//and don't intend to do blen opt
	bool createTopologyOnly=false;
	if(optPrecision < 0.0) createTopologyOnly=true;

	TreeNode* cut = allNodes[cutnum];
	TreeNode *broken = allNodes[broke->nodeNum];
	TreeNode *connector=NULL;
	TreeNode* prunePoint=NULL;
	TreeNode *sib;
	//note that this assignment of the sib can be overridden below if cut is attached to the root or the subtreeNode
	if(cut->next!=NULL) sib=cut->next;
	else sib=cut->prev;

	//determine who the connector node will be.  It will be cut->anc unless that is the root
	//if cut->anc is the root, connector will be one of cut's siblings, which is freed when
	//the basal trichotomy is reestablished after removing cut.
	if(cut->anc->IsNotRoot()){
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
		if(root->left!=cut && root->left->IsInternal()) connector = root->left;
		else if(root->left->next!=cut && root->left->next->IsInternal()) connector = root->left->next;
		else if(root->right!=cut && root->right->IsInternal()) connector = root->right;
		else{//this should be quite rare, and means that the three descendents of the root
			//are cut and two terminals, so no viable swap exists, just try again
			return -1;
			}
		}
	
	//all clas below cut will need to be recalced
	if(createTopologyOnly == false) SweepDirtynessOverTree(cut);
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
			SetBranchLength(replaceForConn, min(max_brlen, replaceForConn->dlen+connector->dlen));
			connector->SubstituteNodeWithRespectToAnc(replaceForConn);
		   	}
		else{//cut is attached to the subtreeNode, so we will have to use it's sib as the connector
			//connector's two children become the subtreeNodes new children, and connector's dlen gets added to subtreeNodes
			TreeNode *subnode=allNodes[subtreeNode];
			SetBranchLength(subnode, min(max_brlen, subnode->dlen+connector->dlen));
			SweepDirtynessOverTree(connector);
			subnode->left=connector->left;
			subnode->right=connector->right;
			connector->left->anc=subnode;
			connector->right->anc=subnode;
			}
		}
	else{//cut is connected to the root so we need to steal a non terminal sib node as the connector
		if(createTopologyOnly == false) MakeNodeDirty(root);
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
			SetBranchLength(root->left, min(max_brlen, root->left->dlen+connector->dlen));
			sib=root->left;
			}
		else{
			SetBranchLength(root->right, min(max_brlen, root->right->dlen+connector->dlen));
			sib=root->right;
			}
			
		//add the connectors two desccendants as descendants of the root
		assert(connector->right==connector->left->next);
		connector->SubstituteNodeWithRespectToAnc(connector->left);
		root->AddDes(connector->right);
		}
	
	//establish correct topology for connector and cut nodes
	if(createTopologyOnly == false) MakeNodeDirty(connector);
	cut->anc=connector;
	connector->left=connector->right=cut;
	connector->next=connector->prev=connector->anc=cut->next=cut->prev=NULL;

	broken->SubstituteNodeWithRespectToAnc(connector);
	connector->AddDes(broken);
	assert(connector->right == broken);

	SetBranchLength(connector, max(min_brlen, broken->dlen*ZERO_POINT_FIVE));
	SetBranchLength(broken, connector->dlen);

	if(createTopologyOnly == false){
		SweepDirtynessOverTree(connector, cut);
		if(broke->reconDist > 1)
			OptimizeBranchesWithinRadius(connector, optPrecision, subtreeNode, sib);
		else 
			OptimizeBranchesWithinRadius(connector, optPrecision, subtreeNode, NULL);
		}
	return 0;
}

void Tree::ReorientSubtreeSPRMutateDummy(int oroot, ReconNode *nroot, FLOAT_TYPE optPrecision){
	//this is just a spoof version of SPRMutate that will perform the same mutation
	//several times with different optimiation settings, but will otherwise 
	//maintain exactly the same program flow because it resets the seed
#ifndef VARIABLE_OPTIMIZATION
	assert(0);
#else
	Individual tempIndiv;
	tempIndiv.treeStruct=new Tree();
	
	Individual sourceIndiv;
	sourceIndiv.treeStruct=this;
	sourceIndiv.mod->CopyModel(this->mod);
		
	int savedSeed;
	var.precision(10);
	var << "reSPR" << "\t" << nroot->reconDist << "\t" << lnL << "\t";
	//FLOAT_TYPE prec[5]={(FLOAT_TYPE).01, (FLOAT_TYPE).5, (FLOAT_TYPE).01, (FLOAT_TYPE).01, (FLOAT_TYPE).01};
	
	tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, &sourceIndiv);
	FLOAT_TYPE origThresh = treeRejectionThreshold;
/*		
	treeRejectionThreshold = 10000;

	for(int i=0;i<1;i++){
		savedSeed = rnd.seed();
		optCalcs = 0;
		tempIndiv.treeStruct->ReorientSubtreeSPRMutate(oroot, nroot, optPrecision);
		var << tempIndiv.treeStruct->lnL << "\t" << optCalcs << "\t";
		optCalcs = 0;
		rnd.set_seed(savedSeed);
		tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, &sourceIndiv, true);
		}
	treeRejectionThreshold = -10000;
	for(int i=0;i<1;i++){
		savedSeed = rnd.seed();
		optCalcs = 0;
		tempIndiv.treeStruct->ReorientSubtreeSPRMutate(oroot, nroot, optPrecision);
		var << tempIndiv.treeStruct->lnL << "\t" << optCalcs << "\t";
		optCalcs = 0;
		rnd.set_seed(savedSeed);
		tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, &sourceIndiv, true);
		}
*/
		/*	for(int i=0;i<1;i++){
		treeRejectionThreshold = origThresh;
		savedSeed = rnd.seed();
		optCalcs = 0;
		tempIndiv.treeStruct->ReorientSubtreeSPRMutate(oroot, nroot, optPrecision);
		var << tempIndiv.treeStruct->lnL << "\t" << optCalcs << "\t";
		optCalcs = 0;
		rnd.set_seed(savedSeed);
		tempIndiv.CopySecByRearrangingNodesOfFirst(tempIndiv.treeStruct, &sourceIndiv, true);
		}
*/
	treeRejectionThreshold = origThresh;
	tempIndiv.treeStruct->RemoveTreeFromAllClas();
	delete tempIndiv.treeStruct;
	tempIndiv.treeStruct=NULL;
	sourceIndiv.treeStruct=NULL;
	optCalcs = 0;
	ReorientSubtreeSPRMutate(oroot, nroot, optPrecision);
	var << lnL << "\t" << optCalcs << "\n";
	optCalcs = 0;
#endif
	}

void Tree::ReorientSubtreeSPRMutate(int oroot, ReconNode *nroot, FLOAT_TYPE optPrecision){
	//this is used to allow the other half of SPR rearrangements in which
	//the part of the tree containing the root is considered the subtree
	//to be attached.  Terminology is VERY confusing here. newRoot is the 
	//branch to be bisected (rooted at).  oldRoot is the node that is at the 
	//base of the subtree currently.  After the rearrangement it will still
	//be at the base of the subtree, but in the middle of a different branch

	//if the optPrecision passed in is < 0 it means that we're just trying to 
	//make the tree structure for some reason, but don't have CLAs allocated
	//and don't intend to do blen opt
	bool createTopologyOnly=false;
	if(optPrecision < 0.0) createTopologyOnly=true;

	TreeNode *newroot=allNodes[nroot->nodeNum];
	TreeNode *oldroot=allNodes[oroot];

	//these are the only blens that need to be dealt with specially
	FLOAT_TYPE fusedBlen = min(max_brlen, oldroot->left->dlen + oldroot->right->dlen);
	FLOAT_TYPE dividedBlen = max(ZERO_POINT_FIVE * newroot->dlen, min_brlen);

	//first detatch the subtree and make it free floating.  This will
	//leave oroot in its place and fuse two branches in the subtree
	//into a branch connecting one of oroots des to its other des
	//This makes that des a tricotomy with a NULL anc.  Then the rotating 
	//begins.
	if(createTopologyOnly == false){
		SweepDirtynessOverTree(oldroot->left);
		SweepDirtynessOverTree(oldroot->right);
		}
	
	TreeNode *prunePoint;
	TreeNode *tempRoot;
	if(oldroot->left->IsInternal()){
		tempRoot=oldroot->left;
		prunePoint=oldroot->right;
		}
	else{
		tempRoot=oldroot->right;
		prunePoint=oldroot->left;
		}

	tempRoot->AddDes(prunePoint);
	//prunePoint->dlen=fusedBlen;
	SetBranchLength(prunePoint, fusedBlen);
	tempRoot->anc=NULL;

	if(createTopologyOnly == false) MakeNodeDirty(tempRoot);
	
	//collect each of the nodes that will need to be flipped
	vector<TreeNode *> path;
	path.reserve(10);
	TreeNode *tmp=newroot->anc;
	while(tmp){
		path.push_back(tmp);
		tmp=tmp->anc;
		}
	reverse(path.begin(),path.end());

	for(vector<TreeNode*>::iterator it=path.begin();(it+1)!=path.end();it++){
		(*it)->MoveDesToAnc(*(it+1));
		}

	//now disconnect the oldroot
	oldroot->left = NULL;
	oldroot->right = NULL;
	
	//and add the new des
	TreeNode *oldanc=newroot->anc;
	oldanc->RemoveDes(newroot);
	oldroot->AddDes(oldanc);
	oldroot->AddDes(newroot);

	SetBranchLength(oldroot->left, dividedBlen);
	SetBranchLength(oldroot->right, dividedBlen);
 
	if(createTopologyOnly == false){
		SweepDirtynessOverTree(newroot);
		SweepDirtynessOverTree(oldroot);
		SweepDirtynessOverTree(tempRoot);
		SweepDirtynessOverTree(prunePoint);
		if(nroot->reconDist > 1) OptimizeBranchesWithinRadius(oldroot, optPrecision, 0, prunePoint);
		else OptimizeBranchesWithinRadius(oldroot, optPrecision, 0, NULL);
		}
	}

void Tree::LoadConstraints(ifstream &con, int nTaxa){
	string temp;//=new char[numTipsTotal + 100];
	Constraint constr;
	int conNum=0;
	do{
		getline(con, temp);
		//getline works strangely on some compilers.  temp should end with ; or \0 , but 
		//might end with \r or \n
		size_t len=temp.length();
		char last=temp.c_str()[len-1];
		if(last == '\r' || last == '\n'){
			temp.erase(len-1, 1);
			len--;
			}
		if(temp[0] != '\0'){
			if(temp[0] != '+' && temp[0] != '-') throw ErrorException("constraint string must start with \'+\' (positive constraint) or \'-\' (negative constraint)");
			if(temp[1] == '.' || temp[1] == '*'){//if individual biparts are specified in *. format
				if(len != nTaxa+1) throw ErrorException("constraint # %d does not have the correct number of characters!\n(has %d) constraint strings must start with \n\'+\' (positive constraint) or \'-\' (negative constraint)\nfollowed by either a ...*** type specification\nor a constraint in newick format", conNum, len);
				constr.ReadDotStarConstraint(temp.c_str());
				constraints.push_back(constr);
				conNum++;
				}
			else if(temp[1] == '('){//if a constraint tree in parenthetical notation is used
				bool numericalTaxa=true;
				for(unsigned i=0;i<len;i++){//see if we are dealing with a treestring with taxa as # or names
					if(isalpha(temp[i])){
						numericalTaxa=false;
						break;
						}
					}
				bool pos;
				if(temp[0] == '+') pos=true;
				else pos=false;
				Tree contree(temp.c_str()+1, numericalTaxa, true);
				contree.CalcBipartitions();
				vector<Bipartition> bip;
				contree.root->GatherConstrainedBiparitions(bip);
				if(pos==false && (bip.size() > 1)) throw ErrorException("Sorry, GARLI can currently only handle a single negatively (conversely) constrainted branch (bipartition):-(");
				//BACKBONE - see if all taxa appear in this constraint or if its a backbone
				if(contree.numTipsAdded < contree.numTipsTotal){
					Bipartition mask = *(contree.root->bipart);
					//complement the mask if necessary
					TreeNode *n=contree.root;
					while(n->IsInternal()) n = n->left;
					if(mask.ContainsTaxon(n->nodeNum) == false) mask.Complement();

					for(vector<Bipartition>::iterator bit=bip.begin();bit!=bip.end();bit++){
						constraints.push_back(Constraint(&(*bit), &mask, pos));
						conNum++;
						}					
					}
				else{
					for(vector<Bipartition>::iterator bit=bip.begin();bit!=bip.end();bit++){
						constraints.push_back(Constraint(&(*bit), pos));
						conNum++;
						}
					}
				}
			else{
				throw ErrorException("problem with constraint # %d\nconstraint strings must start with \n\'+\' (positive constraint) or \'-\' (negative constraint)\nfollowed by either a ...*** type specification\nor a constraint in newick format", conNum, len);
				}
			}
		}while(con.eof() == false);

	//make sure the constraints are compatible with each other!
	if(conNum > 1){
		for(vector<Constraint>::iterator first=constraints.begin();first!=constraints.end();first++){
			for(vector<Constraint>::iterator sec=first+1;sec!=constraints.end();sec++){
				if((*first).IsPositive() != (*sec).IsPositive()) throw ErrorException("cannot mix positive and negative constraints!");
				if(((*first).IsPositive()==false) && ((*sec).IsPositive()==false)) throw ErrorException("Sorry, GARLI can currently only handle a single negatively (conversely) constrainted branch :-(");
				if((*first).IsCompatibleWithConstraint(&(*sec)) == false) throw ErrorException("constraints are not compatible with one another!");
			}
			}
		}
	}

bool Tree::AllowedByConstraint(Constraint *constr, TreeNode *cut, ReconNode *broken, Bipartition &proposed) {
	//Bipartition proposed;
	proposed.FillWithXORComplement(cut->bipart, allNodes[broken->nodeNum]->bipart);

	if(constr->IsBackbone()){
		bool compat = constr->IsCompatibleWithConstraint(&proposed);
		if(compat == false) return compat;

		//this screws up the biparitions, so they need to be recalculated before returning
		if(cut->anc->IsNotRoot()) cut->anc->RecursivelyAddOrRemoveSubtreeFromBipartitions(cut->bipart);
		if(allNodes[broken->nodeNum]->anc->IsNotRoot()) allNodes[broken->nodeNum]->anc->RecursivelyAddOrRemoveSubtreeFromBipartitions(cut->bipart);
		
		if(root->left->IsInternal()) compat=RecursiveAllowedByPositiveConstraintWithMask(constr, constr->GetBackboneMask(), root->left);
		if(compat == true) if(root->left->next->IsInternal()) compat=RecursiveAllowedByPositiveConstraintWithMask(constr, constr->GetBackboneMask(), root->left->next);
		if(compat == true) if(root->right->IsInternal()) compat=RecursiveAllowedByPositiveConstraintWithMask(constr, constr->GetBackboneMask(), root->right);

		CalcBipartitions();//fix the bipartitions
		return compat;
		}

	else{
		if(constr->IsPositive())
			return constr->IsCompatibleWithConstraint(&proposed);
		else{//this is trickier with a negative constraint
			bool compat=constr->IsCompatibleWithConstraint(&proposed);
			if(compat==false) return compat;
			else{//here we need to check if the removal of the cut subtree would create the unallowed bipart
				//I think this could happen about anywhere in the tree, so the cleanest way i see to check
				//is to actually make the tree and verify that it doesn't have the constrained bipart
				
				Tree propTree;
				propTree.MimicTopo(this);

				if(broken->withinCutSubtree == false)  propTree.SPRMutate(cut->nodeNum, broken, -1.0, 0);
				else propTree.ReorientSubtreeSPRMutate(cut->nodeNum, broken, -1.0);

				propTree.CalcBipartitions();

				bool containsBip = (propTree.ContainsBipartitionOrComplement(constr->GetBipartition()) != NULL);
				return (containsBip == false);
				}
			}
		}
	}

bool Tree::AllowedByPositiveConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *cut, TreeNode *broken){
	Bipartition proposed;
	proposed.FillWithXORComplement(cut->bipart, broken->bipart);
	bool compat = constr->IsCompatibleWithConstraintWithMask(&proposed, mask);
	if(compat==false) return compat;
	//this is a little sneaky here.  Cut has not been added to the tree, but since we are going up from broken
	//and it is present in the mask it will effectively appear in biparts in that direction
	else if(broken->IsInternal()){
		compat=RecursiveAllowedByPositiveConstraintWithMask(constr, mask, broken);
		}
	return compat;
	}

bool Tree::AllowedByPositiveBackboneConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *cut, TreeNode *broken){
	Bipartition jointMask=*(constr->GetBackboneMask());
	jointMask.AndEquals(mask);

	Bipartition proposed;
	proposed.FillWithXORComplement(cut->bipart, broken->bipart);
	bool compat = constr->IsCompatibleWithConstraintWithMask(&proposed, &jointMask);
	if(compat==false) return compat;
	else{
		if(broken->anc->IsNotRoot()) broken->anc->RecursivelyAddOrRemoveSubtreeFromBipartitions(cut->bipart);

		if(root->left->IsInternal()) compat=RecursiveAllowedByPositiveConstraintWithMask(constr, &jointMask, root->left);
		if(compat==false) return compat;
		if(root->left->next->IsInternal()) compat=RecursiveAllowedByPositiveConstraintWithMask(constr, &jointMask, root->left->next);
		if(compat==false) return compat;
		if(root->right->IsInternal()) compat=RecursiveAllowedByPositiveConstraintWithMask(constr, &jointMask, root->right);
		if(compat==false) return compat;
		}
	return compat;
	}

bool Tree::AllowedByNegativeConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *cut, TreeNode *broken){
	Bipartition proposed;
	proposed.FillWithXORComplement(cut->bipart, broken->bipart);
	bool compat = constr->IsCompatibleWithConstraintWithMask(&proposed, mask);
	if(compat==true) return compat;
	else if(broken->IsInternal()) compat=RecursiveAllowedByNegativeConstraintWithMask(constr, mask, broken);
	return compat;
	}

bool Tree::RecursiveAllowedByPositiveConstraintWithMask(Constraint *constr, const Bipartition *mask, TreeNode *nd){
	bool compat = constr->IsCompatibleWithConstraintWithMask(nd->bipart, mask);
	if(compat==false) return compat;
	else if(nd->left->IsInternal()) compat=RecursiveAllowedByPositiveConstraintWithMask(constr, mask, nd->left);

	if(compat==false) return compat;
	else if(nd->left->next->IsInternal()) compat=RecursiveAllowedByPositiveConstraintWithMask(constr, mask, nd->left->next);
	
	if(compat==false) return compat;
	if(nd->left->next->next != NULL)
		if(nd->left->next->next->IsInternal())
			compat=RecursiveAllowedByPositiveConstraintWithMask(constr, mask, nd->left->next);

	return compat;
	}

bool Tree::RecursiveAllowedByNegativeConstraintWithMask(Constraint *constr, Bipartition *mask, TreeNode *nd){
	bool compat = constr->IsCompatibleWithConstraintWithMask(nd->bipart, mask);
	if(compat==true) return compat;
	else if(nd->left->IsInternal()) compat=RecursiveAllowedByNegativeConstraintWithMask(constr, mask, nd->left);

	if(compat==true) return compat;
	else if(nd->left->next->IsInternal()) compat=RecursiveAllowedByNegativeConstraintWithMask(constr, mask, nd->left->next);
	
	if(compat==true) return compat;
	if(nd->left->next->next != NULL)
		if(nd->left->next->next->IsInternal())
			compat=RecursiveAllowedByNegativeConstraintWithMask(constr, mask, nd->left->next);

	return compat;
	}

//DJZ 8-11-04  This version is only for the master doing SPRs on nodes that aren't in a subtree when subtree
//mode is on.  Basically the only difference is that if the ancestor of the cut node is the root, we need to 
//choose one of the other nonSubtree nodes to make a connector to avoid screwing up the subtree partitioning
void Tree::SPRMutate(int cutnum, int broknum, FLOAT_TYPE optPrecision, const vector<int> &nonSubNodes)
{	assert( numBranchesAdded > 3 );
	assert(0);//needst to be verified
	TreeNode* cut = allNodes[cutnum];
	assert(cut!=NULL);
	SweepDirtynessOverTree(cut->anc);
	TreeNode *connector;

	if(cut->anc->IsNotRoot()){
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

	if(broken->dlen*ZERO_POINT_FIVE > min_brlen){
		connector->dlen=broken->dlen*ZERO_POINT_FIVE;
		broken->dlen-=connector->dlen;
		}
	else connector->dlen=broken->dlen=min_brlen;

	SweepDirtynessOverTree(connector, cut);
	MakeNodeDirty(connector);
	
#ifdef OPT_DEBUG
	opt << "SPR\n";
#endif
	OptimizeBranchesWithinRadius(connector, optPrecision, 0, NULL);
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
	
	if(from->left->IsInternal()) CopyClaIndecesInSubtree(from->left, remove);
	if(from->right->IsInternal()) CopyClaIndecesInSubtree(from->right, remove);
	}

void Tree::DirtyNodesInSubtree(TreeNode *nd){
	
	MakeNodeDirty(nd);
	if(nd->left->IsInternal()) DirtyNodesInSubtree(nd->left);
	if(nd->right->IsInternal()) DirtyNodesInSubtree(nd->right);
	
	}

bool RescaleRateHet(CondLikeArray *destCLA, int nsites, int nRateCats){


		FLOAT_TYPE *destination=destCLA->arr;
		int *underflow_mult=destCLA->underflow_mult;

		//check if any clas are getting close to underflow
#ifdef UNIX
		madvise(destination, sizeof(FLOAT_TYPE)*4*nRateCats*nsites, MADV_SEQUENTIAL);
		madvise(underflow_mult, sizeof(int)*nsites, MADV_SEQUENTIAL);
#endif
		bool reduceRescale=false;
		FLOAT_TYPE small1, large1, small2, large2;
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
			small1= (destination[0] < destination[2] ? destination[0] : destination[2]);
			large1= (destination[0] > destination[2] ? destination[0] : destination[2]);
			small2= (destination[1] < destination[3] ? destination[1] : destination[3]);
			large2= (destination[1] > destination[3] ? destination[1] : destination[3]);
			small1 = (small1 < small2 ? small1 : small2);
			large1= (large1 > large2 ? large1 : large2);

			for(int r=1;r<nRateCats;r++){
				small2= (destination[0 + r*4] < destination[2 + r*4] ? destination[0 + r*4] : destination[2 + r*4]);
				large2= (destination[0 + r*4] > destination[2 + r*4] ? destination[0 + r*4] : destination[2 + r*4]);
				small1 = (small1 < small2 ? small1 : small2);
				large1= (large1 > large2 ? large1 : large2);				
				small2= (destination[1 + r*4] < destination[3 + r*4] ? destination[1 + r*4] : destination[3 + r*4]);
				large2= (destination[1 + r*4] > destination[3 + r*4] ? destination[1 + r*4] : destination[3 + r*4]);
				small1 = (small1 < small2 ? small1 : small2);
				large1= (large1 > large2 ? large1 : large2);	
				}			
#endif
#ifdef SINGLE_PRECISION_FLOATS
			if(large1< 1.0f){
				if(large1 < 1e-30f){
					throw(1);
					}
				int index=0;
				while(precalcThresh[index] > large1){
					index++;
					}
				int incr = precalcIncr[index];
				int incr2=((int) -log(large1));
				underflow_mult[i]+=incr;
				FLOAT_TYPE mult=precalcMult[index];
				//FLOAT_TYPE mult=exp((FLOAT_TYPE)incr);
				for(int r=0;r<nRateCats;r++){
					for(int q=0;q<4;q++){
						destination[r*4 + q]*=mult;
						assert(destination[r*4 +q] == destination[r*4 +q]);
						assert(destination[r*4 +q] < 1e10);
						}
					}
				}
#else
			if(large1< 1e-5){
				if(large1 < 1e-150){
					throw(1);
					}
				int incr=((int) -log(large1))+5;
				underflow_mult[i]+=incr;
				FLOAT_TYPE mult=exp((FLOAT_TYPE)incr);
				for(int r=0;r<nRateCats;r++){
					for(int q=0;q<4;q++){
						destination[r*4 + q]*=mult;
						assert(destination[r*4 +q] == destination[r*4 +q]);
						assert(destination[r*4 +q] < 1e50);
						}
					}
				}
#endif
			destination+= 4*nRateCats;
			}

		destCLA->rescaleRank=0;
		return reduceRescale;
		}

bool RescaleRateHetNState(CondLikeArray *destCLA, int nsites, int nstates, int nRateCats){

		FLOAT_TYPE *destination=destCLA->arr;
		int *underflow_mult=destCLA->underflow_mult;

		//check if any clas are getting close to underflow
#ifdef UNIX
		madvise(destination, sizeof(FLOAT_TYPE)*nstates*nRateCats*nsites, MADV_SEQUENTIAL);
		madvise(underflow_mult, sizeof(int)*nsites, MADV_SEQUENTIAL);
#endif
		bool reduceRescale=false;
		FLOAT_TYPE small1, large1;
		for(int i=0;i<nsites;i++){
			small1 = (destination[0] < destination[1]) ? destination[0] :  destination[1];
			large1 = (destination[0] > destination[1]) ? destination[0] :  destination[1];
			for(int s=2;s<nstates*nRateCats;s++){
				small1 = (destination[s] < small1) ? destination[s] : small1;
				large1 = (destination[s] > large1) ? destination[s] : large1;
				}

			assert(large1 < 1e10);
			if(large1< 1e-5){
				if(large1 < 1e-150){
					throw(1);
					}
				int incr=((int) -log(large1))+5;
				underflow_mult[i]+=incr;
				FLOAT_TYPE mult=exp((FLOAT_TYPE)incr);
				for(int q=0;q<nstates*nRateCats;q++){
					destination[q]*=mult;
					assert(destination[q] == destination[q]);
					assert(destination[q] < 1e50);
					}
				}
			destination+= nstates*nRateCats;
			}

		destCLA->rescaleRank=0;
		return reduceRescale;
		}

int Tree::ConditionalLikelihoodRateHet(int direction, TreeNode* nd, bool fillFinalCLA /*=false*/){
	//note that fillFinalCLA just refers to whether we actually want to calc a CLA
	//representing the contribution of the entire tree vs just calcing the score
	//The only reason I can think of for doing that is to calc internal state probs
	//the fuction will then return a pointer to the CLA

	assert(this != NULL);

	calcCount++;
	
	int nsites = data->NChar();

	CondLikeArray *destCLA=NULL;

	TreeNode* Lchild, *Rchild;
	CondLikeArray *LCLA=NULL, *RCLA=NULL, *partialCLA=NULL;

	vector<FLOAT_TYPE> Rprmat(mod->NStates() * mod->NStates() * mod->NRateCats());
	vector<FLOAT_TYPE> Lprmat(mod->NStates() * mod->NStates() * mod->NRateCats());

	if(direction != ROOT){
		//the only complicated thing here will be to set up the two children depending on the direction
		//get all of the clas, underflow mults and pmat set up here, then the actual calc loops below 
		//won't depend on direction
		if(direction==DOWN){
			Lchild=nd->left;
			Rchild=nd->right;
			
			if(Lchild->IsInternal())
				LCLA=GetClaDown(Lchild);
			if(Rchild->IsInternal())
				RCLA=GetClaDown(Rchild);
			
			mod->CalcPmat(Lchild->dlen, &Lprmat[0], false);
			mod->CalcPmat(Rchild->dlen, &Rprmat[0], false);
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

				mod->CalcPmat(nd->dlen, &Lprmat[0], false);
			
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
				if(Lchild->IsInternal())
					LCLA=GetClaDown(Lchild);

				mod->CalcPmat(Lchild->dlen, &Lprmat[0], false);
				}
			
			if(Rchild->IsInternal())
				RCLA=GetClaDown(Rchild);
			
			mod->CalcPmat(Rchild->dlen, &Rprmat[0], false);
			}

		if(direction==DOWN) destCLA=GetClaDown(nd, false);
		else if(direction==UPRIGHT) destCLA=GetClaUpRight(nd, false);
		else if(direction==UPLEFT) destCLA=GetClaUpLeft(nd, false);

		if(LCLA!=NULL && RCLA!=NULL){
			//two internal children
			ProfIntInt.Start();
#ifdef EQUIV_CALCS
			if(direction==DOWN){
				CalcFullCLAInternalInternalEQUIV(destCLA, LCLA, RCLA, &Lprmat[0], &Rprmat[0], nsites,  mod->NRateCats(), nd->left->tipData, nd->right->tipData);
				}
			else 
#endif
			if(modSpec.IsNucleotide() == false)
				CalcFullCLAInternalInternalNState(destCLA, LCLA, RCLA, &Lprmat[0], &Rprmat[0], nsites,  mod->NRateCats(), mod->NStates());
			else
				CalcFullCLAInternalInternal(destCLA, LCLA, RCLA, &Lprmat[0], &Rprmat[0], nsites,  mod->NRateCats());
			ProfIntInt.Stop();
			}

		else if(LCLA==NULL && RCLA==NULL){
			//two terminal children
			ProfTermTerm.Start();
			if(modSpec.IsNucleotide() == false)
				CalcFullCLATerminalTerminalNState(destCLA, &Lprmat[0], &Rprmat[0], Lchild->tipData, Rchild->tipData, nsites,  mod->NRateCats(), mod->NStates());
			else
				CalcFullCLATerminalTerminal(destCLA, &Lprmat[0], &Rprmat[0], Lchild->tipData, Rchild->tipData, nsites,  mod->NRateCats());
			ProfTermTerm.Stop();
			}

		else{
			//one terminal, one internal
			ProfIntTerm.Start();

			if(modSpec.IsNucleotide() == false){
				if(LCLA==NULL)
					CalcFullCLAInternalTerminalNState(destCLA, RCLA, &Rprmat[0], &Lprmat[0], Lchild->tipData, nsites,  mod->NRateCats(), mod->NStates());
				else 
					CalcFullCLAInternalTerminalNState(destCLA, LCLA, &Lprmat[0], &Rprmat[0], Rchild->tipData, nsites,  mod->NRateCats(), mod->NStates());
				}
			else{
#ifdef OPEN_MP
				if(LCLA==NULL)
					CalcFullCLAInternalTerminal(destCLA, RCLA, &Rprmat[0], &Lprmat[0], Lchild->tipData, nsites,  mod->NRateCats(), Lchild->ambigMap);
				else
					CalcFullCLAInternalTerminal(destCLA, LCLA, &Lprmat[0], &Rprmat[0], Rchild->tipData, nsites,  mod->NRateCats(), Rchild->ambigMap);
				}
#else
				if(LCLA==NULL)
					CalcFullCLAInternalTerminal(destCLA, RCLA, &Rprmat[0], &Lprmat[0], Lchild->tipData, nsites,  mod->NRateCats());
				else 
					CalcFullCLAInternalTerminal(destCLA, LCLA, &Lprmat[0], &Rprmat[0], Rchild->tipData, nsites,  mod->NRateCats());
				}
#endif
			ProfIntTerm.Stop();
			}
		}
	
	if(direction==ROOT){
		//at the root we need to include the contributions of 3 branches.  Check if we have a
		//valid CLA that already represents two of these three. If so we can save a bit of
		//computation.  This will mainly be the case during blen optimization, when when we 
		//only change one of the branches again and again.
		vector<FLOAT_TYPE> prmat(mod->NStates() * mod->NStates() * mod->NRateCats());
		TreeNode *child;
		CondLikeArray *childCLA=NULL;
		int *childUnderMult=NULL;
		
		if(claMan->IsDirty(nd->claIndexUL) == false){
			partialCLA=GetClaUpLeft(nd, false);
			child=nd->left;
			if(child->IsInternal()){
				childCLA=GetClaDown(child, true);
				}
			mod->CalcPmat(child->dlen, &prmat[0], false);
			}
		else if(claMan->IsDirty(nd->claIndexUR) == false){
			partialCLA=GetClaUpRight(nd, false);
			child=nd->right;
			if(child->IsInternal()){
				childCLA=GetClaDown(child, true);
				}
			mod->CalcPmat(child->dlen, &prmat[0], false);
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
				mod->CalcPmat(nd->dlen, &prmat[0], false);
				}
			else{
				child=nd->left->next;
				if(child->IsInternal()){
					childCLA=GetClaDown(child, true);
					}
				mod->CalcPmat(child->dlen, &prmat[0], false);
				}
			}	
		
		if(fillFinalCLA==false){
			if(childCLA!=NULL){//if child is internal
				ProfScoreInt.Start();
				if(modSpec.IsNucleotide())
					lnL = GetScorePartialInternalRateHet(partialCLA, childCLA, &prmat[0]);
				else
					lnL = GetScorePartialInternalNState(partialCLA, childCLA, &prmat[0]);
					
				ProfScoreInt.Stop();
				}	
			else{
				ProfScoreTerm.Start();
				if(modSpec.IsNucleotide())
					lnL = GetScorePartialTerminalRateHet(partialCLA, &prmat[0], child->tipData);
				else
					lnL = GetScorePartialTerminalNState(partialCLA, &prmat[0], child->tipData);

				ProfScoreTerm.Stop();
				}
			}
		
		else{
			//this is only for inferring internal states
			//careful!  This will have to be returned manually!!
			int wholeTreeIndex=claMan->AssignClaHolder();
			claMan->FillHolder(wholeTreeIndex, ROOT);
			claMan->ReserveCla(wholeTreeIndex);
			if(childCLA!=NULL)//if child is internal
				CalcFullCLAPartialInternalRateHet(claMan->GetCla(wholeTreeIndex), childCLA, &prmat[0], partialCLA, nsites,  mod->NRateCats());
			else
				CalcFullCLAPartialTerminalRateHet(claMan->GetCla(wholeTreeIndex), partialCLA, &prmat[0], child->tipData, nsites,  mod->NRateCats());
				
			return wholeTreeIndex;
			}
		}

	if(direction != ROOT)
		if(destCLA->rescaleRank >= rescaleEvery){
			ProfRescale.Start();
			if(modSpec.IsNucleotide())
				RescaleRateHet(destCLA, nsites, mod->NRateCats());
			else
				RescaleRateHetNState(destCLA, nsites, mod->NStates(), mod->NRateCats());

			ProfRescale.Stop();
			}
	return -1;
	}

int Tree::Score(int rootNodeNum /*=0*/){

	TreeNode *rootNode=allNodes[rootNodeNum];

#ifdef EQUIV_CALCS
	if(dirtyEQ){
		ProfEQVectors.Start();
		root->SetEquivalentConditionalVectors(data);
		ProfEQVectors.Stop();
		dirtyEQ=false;
		}
#endif

	bool scoreOK=true;
	do{
		try{
			scoreOK=true;
/*			if(mod->NRateCats()==1)
				ConditionalLikelihood( ROOT, rootNode);
			else
*/				ConditionalLikelihoodRateHet( ROOT, rootNode);
			}
#if defined(NDEBUG)
			catch(int){
#else
			catch(int err){
#endif
				assert(err==1);
				scoreOK=false;
				MakeAllNodesDirty();
				rescaleEvery -= 2;
				ofstream resc("rescale.log", ios::app);
				resc << "rescale reduced to " << rescaleEvery << endl;
				resc.close();
				if(rescaleEvery<2) if(rescaleEvery<2) throw(ErrorException("Problem with rescaling in branchlength optimization.  Please report this error to zwickl@nescent.org."));
				}
		}while(scoreOK==false);

	return 1;
	}
/*
FLOAT_TYPE Tree::SubTreeScore( TreeNode *nd){
	//calculates the likelihood of the tree above the node passed in
	FLOAT_TYPE lnL = 0.0;
	int nSites = data->NChar();
	int ck;

	if(claMan->IsDirty(nd->claIndexDown)){
		if(mod->NRateCats()==1)
		ConditionalLikelihood( DOWN, nd);
		else
			ConditionalLikelihoodRateHet(DOWN, nd);
		}

	FLOAT_TYPE *cla=claMan->GetCla(nd->claIndexDown)->arr;
	int *underflow_mult=claMan->GetCla(nd->claIndexDown)->underflow_mult;

	// loop over all patterns
	long FLOAT_TYPE Lk;
	FLOAT_TYPE siteL;
	int ufcount=0;
	const int *countit=data->GetCounts();
	if(mod->PropInvar()==0.0){
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
		FLOAT_TYPE prI=mod->PropInvar();
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
					siteL = log( Lk * (1.0-prI) + (prI * mod->Pi(conBases[k])) * exp((FLOAT_TYPE)underflow_mult[k]));
					lnL += ( *countit++ * (siteL + underflow_mult[k]));
					}
				}
			}
		else{//gamma and invariants
			FLOAT_TYPE scaledGammaProp=0.25 * (1.0-prI);
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
					siteL = log( Lk * scaledGammaProp + (prI * mod->Pi(conBases[k])) * exp((FLOAT_TYPE)underflow_mult[k]));
					lnL += ( *countit++ * (siteL + underflow_mult[k]));
					}
				}
			}
		}
	return lnL;
	}
*/
/*
FLOAT_TYPE Tree::SubTreeScoreRateHet( TreeNode *nd){
	//calculates the likelihood of the tree above the node passed in
	FLOAT_TYPE sublnL = 0.0;
	int nSites = data->NChar();
	int ck;

	if(claMan->IsDirty(nd->claIndexDown))
		ConditionalLikelihoodRateHet(DOWN, nd);


	FLOAT_TYPE *cla=claMan->GetCla(nd->claIndexDown)->arr;
	int *underflow_mult=claMan->GetCla(nd->claIndexDown)->underflow_mult;

	// loop over all patterns
	long FLOAT_TYPE Lk;
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

	//this will be the case if we are simply making the tree structure but
	//never intend to score it
	if(nd->IsInternal() && nd->claIndexDown == -1){
		return;
		}

	if(from==NULL){
		//if this is the branch where the dirtyness starts
		if(nd->IsInternal()){
			nd->claIndexUL=claMan->SetDirty(nd->claIndexUL);
			nd->claIndexUR=claMan->SetDirty(nd->claIndexUR);
			if(nd->left->IsInternal()) SweepDirtynessOverTree(nd->left, nd);
			if(nd->right->IsInternal()) SweepDirtynessOverTree(nd->right, nd);
			}
		if(nd->anc!=NULL) SweepDirtynessOverTree(nd->anc, nd);	
		}
	else{
	//if the change was below, invalidating clas above, also if the change
	//was on the path connecting to the central des of the root
		if(from==nd->anc || (nd->IsRoot() && from==nd->left->next)){
			nd->claIndexUL=claMan->SetDirty(nd->claIndexUL);
			nd->claIndexUR=claMan->SetDirty(nd->claIndexUR);
			if(nd->left->IsInternal()) SweepDirtynessOverTree(nd->left, nd);
			if(nd->right->IsInternal()) SweepDirtynessOverTree(nd->right, nd);
			}
		else if(from==nd->left){
			nd->claIndexUR=claMan->SetDirty(nd->claIndexUR);
			nd->claIndexDown=claMan->SetDirty(nd->claIndexDown);
			if(nd->right->IsInternal()) SweepDirtynessOverTree(nd->right, nd);
			if(nd->anc!=NULL) SweepDirtynessOverTree(nd->anc, nd);		
			else if(nd->left->next->IsInternal()) SweepDirtynessOverTree(nd->left->next, nd);
			}
		else if(from==nd->right){
			nd->claIndexUL=claMan->SetDirty(nd->claIndexUL);
			nd->claIndexDown=claMan->SetDirty(nd->claIndexDown);
			if(nd->left->IsInternal()) SweepDirtynessOverTree(nd->left, nd);
			if(nd->anc!=NULL) SweepDirtynessOverTree(nd->anc, nd);		
			else if(nd->left->next->IsInternal()) SweepDirtynessOverTree(nd->left->next, nd);
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

void Tree::RotateNodesAtRoot(TreeNode *newroot){
	//DZ 11-3-02 This can be used to rebalance the tree
	//I'm assuming that this will be called with one of the des of the root;
	assert(newroot->anc==root);
	assert(newroot->IsInternal());
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
	if(root->left->IsInternal()){
		root->left->left->CountSubtreeBranchesAndDepth(llb, lls, 3, true);
		root->left->right->CountSubtreeBranchesAndDepth(lrb, lrs, 3, true);
		lb=llb+lrb+2;
		ls=lls+lrs+4;
		}
	if(root->left->next->IsInternal()){
		root->left->next->left->CountSubtreeBranchesAndDepth(mlb, mls, 3, true);
		root->left->next->right->CountSubtreeBranchesAndDepth(mrb, mrs, 3, true);
		mb=mlb+mrb+2;
		ms=mls+mrs+4;
		}
	if(root->right->IsInternal()){
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
	
		if(root->left->IsInternal()){
			root->left->left->CountSubtreeBranchesAndDepth(llb, lls, 3, true);
			root->left->right->CountSubtreeBranchesAndDepth(lrb, lrs, 3, true);
			lb=llb+lrb+2;
			ls=lls+lrs+4;
			}
		if(root->left->next->IsInternal()){
			root->left->next->left->CountSubtreeBranchesAndDepth(mlb, mls, 3, true);
			root->left->next->right->CountSubtreeBranchesAndDepth(mrb, mrs, 3, true);
			mb=mlb+mrb+2;
			ms=mls+mrs+4;
			}
		if(root->right->IsInternal()){
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

void Tree::CalcBipartitions(bool standardize /*=true*/){
	root->CalcBipartition();
	if(standardize)
		root->StandardizeBipartition();
	}
	
void Tree::OutputBipartitions(){
	ofstream out("biparts.log", ios::app);
	root->OutputBipartition(out);
	}
/*
void Tree::SetDistanceBasedBranchLengthsAroundNode(TreeNode *nd){
	FLOAT_TYPE D1, D2, D3, k1, k2, k3, k4, a, b, c;
	TreeNode *T1, *T2, *T3, *T4;

	FindNearestTerminalUp(nd->left, T1, k1);
	FindNearestTerminalUp(nd->right, T2, k2);
	FindNearestTerminalsDown(nd->anc, nd, T3, T4, k3, k4);
//	FindNearestTerminalUp(nd->, T2, k2);

	if(k4<k3){
		T3=T4;
		k3=k4;
		}
#ifdef FLEX_RATES
	assert(0);
#else
	D1=CalculatePDistance(T1->tipData, T2->tipData, data->NChar())/(1.0-mod->PropInvar()) - k1 -k2;
	D2=CalculatePDistance(T1->tipData, T3->tipData, data->NChar())/(1.0-mod->PropInvar()) - k1 -k3;
	D3=CalculatePDistance(T2->tipData, T3->tipData, data->NChar())/(1.0-mod->PropInvar()) - k2 -k3;
#endif
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

void Tree::FindNearestTerminalUp(TreeNode *start, TreeNode *&term, FLOAT_TYPE &dist){
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

void Tree::FindNearestTerminalsDown(TreeNode *start, TreeNode *from, TreeNode *&term1, TreeNode *&term2, FLOAT_TYPE &dist1, FLOAT_TYPE &dist2){
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
*/
void Tree::OptimizeBranchesAroundNode(TreeNode *nd, FLOAT_TYPE optPrecision, int subtreeNode){
	//depricated
	assert(0);
	//this function will optimize the three branches (2 descendents and one anc) connecting
	//to it.  It assumes that everything that is dirty has been marked so.
	//by default there is only a single optimization pass over the three nodes
/*	FLOAT_TYPE precision1, precision2;

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

void Tree::OutputNthClaAcrossTree(ofstream &deb, TreeNode *nd, int site){
	//int site=0;
	int nstates = mod->NStates();
	int index=nstates * mod->NRateCats() * site;
	
	if(nd->IsInternal() && claMan->IsDirty(nd->claIndexDown) == false){
		deb << nd->nodeNum << "\t0\t" << nd->claIndexDown << "\t";
		for(int i=0;i<nstates*mod->NRateCats();i++) deb << claMan->GetCla(nd->claIndexDown)->arr[index+i] << "\t";
		deb << claMan->GetCla(nd->claIndexDown)->underflow_mult[site] <<"\n";
		}
	if(nd->IsInternal() && claMan->IsDirty(nd->claIndexUL) == false){
		deb << nd->nodeNum << "\t1\t" << nd->claIndexUL << "\t";
		for(int i=0;i<nstates*mod->NRateCats();i++) deb << claMan->GetCla(nd->claIndexUL)->arr[index+i] << "\t";
		deb << claMan->GetCla(nd->claIndexUL)->underflow_mult[site] <<"\n";
		}
	if(nd->IsInternal() && claMan->IsDirty(nd->claIndexUR) == false){
		deb << nd->nodeNum << "\t2\t" << nd->claIndexUR << "\t";
		for(int i=0;i<nstates*mod->NRateCats();i++) deb << claMan->GetCla(nd->claIndexUR)->arr[index+i] << "\t";
		deb << claMan->GetCla(nd->claIndexUR)->underflow_mult[site] <<"\n";
		}

	if(nd->IsInternal())
		OutputNthClaAcrossTree(deb, nd->left, site);
	if(nd->next!=NULL)
		OutputNthClaAcrossTree(deb, nd->next, site);
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
	
	if(nd->left->IsInternal() && nd->left->nodeNum != subtreeNode && nd->left->nodeNum != allNodes[subtreeNode]->anc->nodeNum){
		DirtyNodesOutsideOfSubtree(nd->left, subtreeNode);
		}
	if(nd->right->IsInternal() && nd->right->nodeNum != subtreeNode && nd->right->nodeNum != allNodes[subtreeNode]->anc->nodeNum){
		DirtyNodesOutsideOfSubtree(nd->right, subtreeNode);
		}
	if(nd->IsRoot() && nd->left->next->IsInternal() && nd->left->next->nodeNum != subtreeNode && nd->left->next->nodeNum != allNodes[subtreeNode]->anc->nodeNum){
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

void Tree::InferAllInternalStateProbs(const char *ofprefix){
	char filename[80];
	sprintf(filename, "%s.internalstates.log", ofprefix);
	ofstream out(filename);
	out.precision(5);
	RecursivelyCalculateInternalStateProbs(root, out);
	out.close();
	}

void Tree::RecursivelyCalculateInternalStateProbs(TreeNode *nd, ofstream &out){
	if(nd->IsInternal()) RecursivelyCalculateInternalStateProbs(nd->left, out);
	if(nd->IsInternal()) RecursivelyCalculateInternalStateProbs(nd->next, out);
	
	if(nd->IsInternal()){
		int wholeTreeIndex=ConditionalLikelihoodRateHet(ROOT, nd, true);
		vector<InternalState *> *stateProbs=InferStatesFromCla(claMan->GetCla(wholeTreeIndex)->arr, data->NChar(), mod->NRateCats());

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
	
FLOAT_TYPE Tree::CountClasInUse(){
	FLOAT_TYPE inUse=0.0;
	
	if(claMan->IsDirty(root->claIndexDown) == false) inUse += ONE_POINT_ZERO/claMan->GetNumAssigned(root->claIndexDown);
	if(claMan->IsDirty(root->claIndexUL) == false) inUse += ONE_POINT_ZERO/claMan->GetNumAssigned(root->claIndexUL);
	if(claMan->IsDirty(root->claIndexUR) == false) inUse += ONE_POINT_ZERO/claMan->GetNumAssigned(root->claIndexUR);
	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		TreeNode *n=allNodes[i];	
		if(claMan->IsDirty(n->claIndexDown) == false) inUse += ONE_POINT_ZERO/claMan->GetNumAssigned(n->claIndexDown);
		if(claMan->IsDirty(n->claIndexUL) == false) inUse += ONE_POINT_ZERO/claMan->GetNumAssigned(n->claIndexUL);
		if(claMan->IsDirty(n->claIndexUR) == false) inUse += ONE_POINT_ZERO/claMan->GetNumAssigned(n->claIndexUR);		
		}
	return inUse;
	}
	
FLOAT_TYPE Tree::GetScorePartialTerminalNState(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const char *Ldata){

	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	FLOAT_TYPE *partial=partialCLA->arr;
	int *underflow_mult=partialCLA->underflow_mult;

	int nstates = mod->NStates();
	int nRateCats = mod->NRateCats();
	int nchar = data->NChar();
	const int *countit=data->GetCounts();

	const FLOAT_TYPE *rateProb=mod->GetRateProbs();

#ifdef UNIX
	madvise(partial, nchar*nstates*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	FLOAT_TYPE siteL, totallnL=0.0;

	FLOAT_TYPE *freqs = new FLOAT_TYPE[nstates];
	for(int i=0;i<nstates;i++) freqs[i]=mod->StateFreq(i);

#undef OUTPUT_SITELIKES

#ifdef OUTPUT_SITELIKES
	ofstream sit("sitelikes.log");
	sit.precision(15);
#endif

	if(nRateCats == 1){
		for(int i=0;i<nchar;i++){
			siteL = 0.0;
			if(*Ldata < nstates){ //no ambiguity
				for(int from=0;from<nstates;from++){
					siteL += prmat[(*Ldata)+nstates*from] * partial[from] * freqs[from];
					}
				partial += nstates;
				Ldata++;
				}
				
			else if(*Ldata == nstates){ //total ambiguity
				for(int from=0;from<nstates;from++){
					siteL += partial[from] * freqs[from];
					}
				partial += nstates;
				Ldata++;
				}
			else{ //partial ambiguity
				assert(0);
				}

			totallnL += (countit[i] * (log(siteL) - underflow_mult[i]));
	#ifndef OUTPUT_SITELIKES
			}
	#else
//			for(int q=0;q<countit[i];q++)
				sit << siteL << "\t" << underflow_mult[i] << "\t" << totallnL << "\n";
			}
		sit.close();
	#endif
		}
	else{

		vector<FLOAT_TYPE> stateContrib(nstates);

		for(int i=0;i<nchar;i++){
			for(int f=0;f<nstates;f++) stateContrib[f] = 0.0;
			if(*Ldata < nstates){ //no ambiguity
				for(int rate=0;rate<nRateCats;rate++){
					for(int from=0;from<nstates;from++){
						stateContrib[from] += prmat[rate*nstates*nstates + (*Ldata)+nstates*from] * partial[from] * rateProb[rate];
						}
					partial += nstates;
					}
				Ldata++;
				}
				
			else if(*Ldata == nstates){ //total ambiguity
				for(int rate=0;rate<nRateCats;rate++){
					for(int from=0;from<nstates;from++){
						stateContrib[from] += partial[from] * rateProb[rate];
						}
					partial += nstates;
					}
				Ldata++;
				}

			else{ //partial ambiguity
				assert(0);
				}

			siteL = 0.0;
			for(int from=0;from<nstates;from++)
				siteL += stateContrib[from] * freqs[from];
			totallnL += *countit++ * (log(siteL) - underflow_mult[i]);
	#ifndef OUTPUT_SITELIKES
			}
	#else
			for(int q=0;q<*(countit-1);q++)
				sit << siteL << "\t" << underflow_mult[i] << "\t" << totallnL << "\n";
			}
		sit.close();
	#endif
		}


	delete []freqs;
	return totallnL;
	}

FLOAT_TYPE Tree::GetScorePartialTerminalRateHet(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const char *Ldata){

	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	FLOAT_TYPE *partial=partialCLA->arr;
	int *underflow_mult=partialCLA->underflow_mult;
	const int nRateCats=mod->NRateCats();

	int nchar=data->NChar();
#ifdef UNIX
	madvise(partial, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	FLOAT_TYPE siteL, totallnL=0.0;
	FLOAT_TYPE La, Lc, Lg, Lt;
	
	//gamma and invariants
	const int *countit=data->GetCounts();
	
	const FLOAT_TYPE *rateProb=mod->GetRateProbs();

	int lastConst=data->LastConstant();
	const int *conBases=data->GetConstBases();
	FLOAT_TYPE prI=mod->PropInvar();
//	FLOAT_TYPE scaledGammaProp=(1.0-prI) / mod->NRateCats();

	FLOAT_TYPE freqs[4];
	for(int i=0;i<4;i++) freqs[i]=mod->StateFreq(i);

#undef OUTPUT_SITELIKES

#ifdef OUTPUT_SITELIKES
	ofstream sit("sitelikes.log");
#endif

	for(int i=0;i<nchar;i++){
		La=Lc=Lg=Lt=0.0;
		if(*Ldata > -1){ //no ambiguity
			for(int i=0;i<nRateCats;i++){

				La  += prmat[(*Ldata)+16*i] * partial[0] * rateProb[i];
				Lc  += prmat[(*Ldata+4)+16*i] * partial[1] * rateProb[i];
				Lg  += prmat[(*Ldata+8)+16*i] * partial[2] * rateProb[i];
				Lt  += prmat[(*Ldata+12)+16*i] * partial[3] * rateProb[i];
				partial += 4;
				}
			Ldata++;
			}
			
		else if(*Ldata == -4){ //total ambiguity
			for(int i=0;i<nRateCats;i++){
				La += partial[0] * rateProb[i];
				Lc += partial[1] * rateProb[i];
				Lg += partial[2] * rateProb[i];
				Lt += partial[3] * rateProb[i];
				partial += 4;
				}
			Ldata++;
			}
		else{ //partial ambiguity
			char nstates=-1 * *(Ldata++);
			for(int i=0;i<nstates;i++){
				for(int i=0;i<nRateCats;i++){
					La += prmat[(*Ldata)+16*i]  * partial[4*i] * rateProb[i];
					Lc += prmat[(*Ldata+4)+16*i] * partial[1+4*i] * rateProb[i];
					Lg += prmat[(*Ldata+8)+16*i]* partial[2+4*i] * rateProb[i];
					Lt += prmat[(*Ldata+12)+16*i]* partial[3+4*i] * rateProb[i];
					}
				Ldata++;
				}
			partial+=4*nRateCats;
			}
		if((mod->NoPinvInModel() == false) && (i<=lastConst)){
			FLOAT_TYPE btot=0.0;
			if(conBases[i]&1) btot+=freqs[0];
			if(conBases[i]&2) btot+=freqs[1];
			if(conBases[i]&4) btot+=freqs[2];
			if(conBases[i]&8) btot+=freqs[3];
			if(underflow_mult[i]==0)
				siteL  = ((La*freqs[0]+Lc*freqs[1]+Lg*freqs[2]+Lt*freqs[3]) + prI*btot);
			else 
				siteL  = ((La*freqs[0]+Lc*freqs[1]+Lg*freqs[2]+Lt*freqs[3]) + (prI*btot*exp((FLOAT_TYPE)underflow_mult[i])));
			}
		else
			siteL  = ((La*freqs[0]+Lc*freqs[1]+Lg*freqs[2]+Lt*freqs[3]));

		totallnL += (*countit++ * (log(siteL) - underflow_mult[i]));
#ifndef OUTPUT_SITELIKES
		}
#else
		sit << siteL << "\t" << underflow_mult[i] << "\t" << totallnL << "\n";
		}
	sit.close();
#endif
	return totallnL;
	}
	
FLOAT_TYPE Tree::GetScorePartialInternalRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	FLOAT_TYPE *CL1=childCLA->arr;
	FLOAT_TYPE *partial=partialCLA->arr;
	int *underflow_mult1=partialCLA->underflow_mult;
	int *underflow_mult2=childCLA->underflow_mult;

	int nchar=data->NChar();
	const int nRateCats=mod->NRateCats();


#ifdef UNIX
	madvise(partial, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise(CL1, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	FLOAT_TYPE siteL, totallnL=0.0;
	FLOAT_TYPE La, Lc, Lg, Lt;

	//gamma and invariants
	const int *countit=data->GetCounts();

	const FLOAT_TYPE *rateProb=mod->GetRateProbs();
	FLOAT_TYPE prI=mod->PropInvar();
	int lastConst=data->LastConstant();
	const int *conBases=data->GetConstBases();

	FLOAT_TYPE freqs[4];
	for(int i=0;i<4;i++) freqs[i]=mod->StateFreq(i);

#undef OUTPUT_SITELIKES

#ifdef OUTPUT_SITELIKES
	ofstream sit("sitelikes.log");
#endif

	for(int i=0;i<nchar;i++){
		La=Lc=Lg=Lt=0.0;
		for(int r=0;r<nRateCats;r++){
			int rOff=r*16;
			La += ( prmat[rOff ]*CL1[0]+prmat[rOff + 1]*CL1[1]+prmat[rOff + 2]*CL1[2]+prmat[rOff + 3]*CL1[3]) * partial[0] * rateProb[r];
			Lc += ( prmat[rOff + 4]*CL1[0]+prmat[rOff + 5]*CL1[1]+prmat[rOff + 6]*CL1[2]+prmat[rOff + 7]*CL1[3]) * partial[1] * rateProb[r];
			Lg += ( prmat[rOff + 8]*CL1[0]+prmat[rOff + 9]*CL1[1]+prmat[rOff + 10]*CL1[2]+prmat[rOff + 11]*CL1[3]) * partial[2] * rateProb[r];
			Lt += ( prmat[rOff + 12]*CL1[0]+prmat[rOff + 13]*CL1[1]+prmat[rOff + 14]*CL1[2]+prmat[rOff + 15]*CL1[3]) * partial[3] * rateProb[r];
			partial+=4;
			CL1+=4;
			}
		if((mod->NoPinvInModel() == false) && (i<=lastConst)){
			FLOAT_TYPE btot=0.0;
			if(conBases[i]&1) btot+=freqs[0];
			if(conBases[i]&2) btot+=freqs[1];
			if(conBases[i]&4) btot+=freqs[2];
			if(conBases[i]&8) btot+=freqs[3];
			if(underflow_mult1[i] + underflow_mult2[i] == 0)
				siteL  = ((La*freqs[0]+Lc*freqs[1]+Lg*freqs[2]+Lt*freqs[3]) + prI*btot);
			else 
				siteL  = ((La*freqs[0]+Lc*freqs[1]+Lg*freqs[2]+Lt*freqs[3]) + (prI*btot*exp((FLOAT_TYPE)underflow_mult1[i]+underflow_mult2[i])));
			}
		else
			siteL  = ((La*freqs[0]+Lc*freqs[1]+Lg*freqs[2]+Lt*freqs[3]));	
		
		totallnL += (*countit++ * (log(siteL) - underflow_mult1[i] - underflow_mult2[i]));
#ifndef OUTPUT_SITELIKES
		}
#else
		sit << siteL << "\t" << underflow_mult1[i] << "\t" << underflow_mult2[i] << "\t" << totallnL << "\n";
		}
	sit.close();
#endif
	return totallnL;
	}


FLOAT_TYPE Tree::GetScorePartialInternalNState(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat){
	//this function assumes that the pmat is arranged with nstates^2 entries for the
	//first rate, followed by nstate^2 for the second, etc.
	FLOAT_TYPE *CL1=childCLA->arr;
	FLOAT_TYPE *partial=partialCLA->arr;
	int *underflow_mult1=partialCLA->underflow_mult;
	int *underflow_mult2=childCLA->underflow_mult;

	int nchar=data->NChar();
	const int *countit=data->GetCounts();
	const int nRateCats = mod->NRateCats();
	int nstates = mod->NStates();

	const FLOAT_TYPE *rateProb=mod->GetRateProbs();

#ifdef UNIX
	madvise(partial, nchar*nstates*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise(CL1, nchar*nstates*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	FLOAT_TYPE totallnL=0.0;

	FLOAT_TYPE *freqs = new FLOAT_TYPE[nstates];
	for(int i=0;i<nstates;i++) freqs[i]=mod->StateFreq(i);

#undef OUTPUT_SITELIKES

#ifdef OUTPUT_SITELIKES
	ofstream sit("sitelikes.log");
#endif

	if(nRateCats == 1){
		for(int i=0;i<nchar;i++){
			FLOAT_TYPE siteL = 0.0;
			for(int from=0;from<nstates;from++){
				FLOAT_TYPE temp = 0.0;
				for(int to=0;to<nstates;to++){
					temp += prmat[from*nstates + to]*CL1[to];
					}
				siteL += temp * partial[from] * freqs[from];
				}
			totallnL += (*countit++ * (log(siteL) - underflow_mult1[i] - underflow_mult2[i]));
			CL1 += nstates;
			partial += nstates;
#ifndef OUTPUT_SITELIKES
			}
#else
//			for(int q=0;q<*(countit-1);q++)
				sit << siteL << "\t" << underflow_mult1[i] << "\t" << underflow_mult2[i] << "\t" << totallnL << "\n";
			}
		sit.close();
#endif
		}
	else{//this is a lot nastier in the case of rate het, and I don't see a clean way to do it except by writing the
		//likelihood contributions of the various states to a temp array, instead of keeping a running total
		vector<FLOAT_TYPE> stateContrib(nstates);

		for(int i=0;i<nchar;i++){
			FLOAT_TYPE siteL = 0.0;
			for(int f=0;f<nstates;f++) stateContrib[f] = 0.0;
			for(int rate=0;rate<nRateCats;rate++){
				for(int from=0;from<nstates;from++){
					FLOAT_TYPE temp = 0.0;
					for(int to=0;to<nstates;to++){
						temp += prmat[rate*nstates*nstates + from*nstates + to]*CL1[to];
						}
					stateContrib[from] += temp * partial[from] * rateProb[rate];
					}
				CL1 += nstates;
				partial += nstates;
				}
			siteL = 0.0;
			for(int from=0;from<nstates;from++)
				siteL += stateContrib[from] * freqs[from];
			totallnL += (*countit++ * (log(siteL) - underflow_mult1[i] - underflow_mult2[i]));
#ifndef OUTPUT_SITELIKES
			}
#else
			for(int q=0;q<*(countit-1);q++)
				sit << siteL << "\t" << underflow_mult1[i] << "\t" << underflow_mult2[i] << "\t" << totallnL << "\n";
			}
		sit.close();
#endif
		}

	delete []freqs;
	return totallnL;
	}

void Tree::LocalMove(){
	assert(0);
	//This is not working
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
			if(v->IsRoot()){
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
	if(v->IsRoot()){
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
				if(c->IsRoot()){
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
	FLOAT_TYPE m;
	FLOAT_TYPE changing_blens[3];
//	FLOAT_TYPE new_blens[3];
	changing_blens[0]=a->dlen;
	changing_blens[1]=u->dlen;
	if(cPosition==3){
		changing_blens[2]=v->dlen;
		}
	else {
		changing_blens[2]=c->dlen;
		}
	m=changing_blens[0]+changing_blens[1]+changing_blens[2];
	FLOAT_TYPE r=rnd.uniform();



//	FLOAT_TYPE tuning=.25;
//	FLOAT_TYPE tuning=.1;
	FLOAT_TYPE mprime=m*exp((FLOAT_TYPE).5*(rnd.uniform()-(FLOAT_TYPE).5));
	FLOAT_TYPE x, y;
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

void Tree::NNIMutate(int node, int branch, FLOAT_TYPE optPrecision, int subtreeNode){

	assert(0);
	TreeNode* connector=NULL;
	TreeNode* cut=NULL;
	TreeNode* broken=NULL;
	TreeNode* sib=NULL;

	assert(node<numNodesTotal);
	connector = allNodes[node];
	assert(connector->IsInternal());
	
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

	if(broken->dlen*.5 > min_brlen){
		connector->dlen=broken->dlen*(FLOAT_TYPE).5;
		broken->dlen-=connector->dlen;
		}
	else connector->dlen=broken->dlen=min_brlen;

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

	OptimizeBranchesWithinRadius(connector, optPrecision, subtreeNode, NULL);
	}

/*
void Tree::OutputBinaryFormattedTree(ofstream &out) const{
	
	for(int i=0;i<numNodesTotal;i++){
		allNodes[i]->OutputBinaryNodeInfo(out);
		}
	out.write((char*) &lnL, sizeof(FLOAT_TYPE));
	out.write((char*) &numTipsTotal, sizeof(numTipsTotal));
	out.write((char*) &numTipsAdded, sizeof(numTipsAdded));
	out.write((char*) &numNodesAdded, sizeof(numNodesAdded));
	out.write((char*) &numBranchesAdded, sizeof(numBranchesAdded));
	out.write((char*) &numNodesTotal, sizeof(numNodesTotal));
	}
*/

void Tree::OutputBinaryFormattedTree(OUTPUT_CLASS &out) const{
	for(int i=0;i<numNodesTotal;i++){
		allNodes[i]->OutputBinaryNodeInfo(out);
		}
	out.WRITE_TO_FILE(&lnL, sizeof(FLOAT_TYPE), 1);
	out.WRITE_TO_FILE(&numTipsTotal, sizeof(numTipsTotal), 1);
	out.WRITE_TO_FILE(&numTipsAdded, sizeof(numTipsAdded), 1);
	out.WRITE_TO_FILE(&numNodesAdded, sizeof(numNodesAdded), 1);
	out.WRITE_TO_FILE(&numBranchesAdded, sizeof(numBranchesAdded), 1);
	out.WRITE_TO_FILE(&numNodesTotal, sizeof(numNodesTotal), 1);
	}

void Tree::ReadBinaryFormattedTree(FILE *in){
	int dum;

	fread((char*) &dum, sizeof(dum), 1, in);
	allNodes[0]->left = allNodes[dum];

	fread((char*) &dum, sizeof(dum), 1, in);
	allNodes[0]->right = allNodes[dum];

	fread((char*) &dum, sizeof(dum), 1, in);
	if(dum == 0) allNodes[0]->prev = NULL;
	else allNodes[0]->prev = allNodes[dum];
	
	fread((char*) &dum, sizeof(dum), 1, in);
	if(dum == 0) allNodes[0]->next = NULL;
	else allNodes[0]->next = allNodes[dum];
	
	fread((char*) &dum, sizeof(dum), 1, in);
	if(dum == 0) allNodes[0]->anc = NULL;
	else allNodes[0]->anc = allNodes[dum];
	
	fread((char*) &allNodes[0]->dlen, sizeof(FLOAT_TYPE), 1, in);

//	double d;
	for(int i=1;i<=numTipsTotal;i++){
		fread(&dum, sizeof(dum), 1, in);
		if(dum == 0) allNodes[i]->prev = NULL;
		else allNodes[i]->prev = allNodes[dum];

		fread(&dum, sizeof(dum), 1, in);
		if(dum == 0) allNodes[i]->next = NULL;
		else allNodes[i]->next = allNodes[dum];

		//all non-root nodes will have an anc, which might be nodenum 0 (the root)
		//so, don't test for zero here
		fread(&dum, sizeof(dum), 1, in);
		allNodes[i]->anc = allNodes[dum];

		fread(&(allNodes[i]->dlen), sizeof(FLOAT_TYPE), 1, in);
		}

	for(int i=numTipsTotal+1;i<numNodesTotal;i++){
		fread((char*) &dum, sizeof(dum), 1, in);
		allNodes[i]->left = allNodes[dum];
		
		fread((char*) &dum, sizeof(dum), 1, in);
		allNodes[i]->right = allNodes[dum];
		
		fread((char*) &dum, sizeof(dum), 1, in);
		if(dum == 0) allNodes[i]->prev = NULL;
		else allNodes[i]->prev = allNodes[dum];
		
		fread((char*) &dum, sizeof(dum), 1, in);
		if(dum == 0) allNodes[i]->next = NULL;
		else allNodes[i]->next = allNodes[dum];
		
		//all non-root nodes will have an anc, which might be nodenum 0 (the root)
		//so, don't test for zero here
		fread((char*) &dum, sizeof(dum), 1, in);
		allNodes[i]->anc = allNodes[dum];
		
		fread((char*) &allNodes[i]->dlen, sizeof(FLOAT_TYPE), 1, in);
		}

	fread((char*) &lnL, sizeof(FLOAT_TYPE), 1, in);
	fread((char*) &numTipsTotal, sizeof(numTipsTotal), 1, in);
	fread((char*) &numTipsAdded, sizeof(numTipsAdded), 1, in);
	fread((char*) &numNodesAdded, sizeof(numNodesAdded), 1, in);
	fread((char*) &numBranchesAdded, sizeof(numBranchesAdded), 1, in);
	fread((char*) &numNodesTotal, sizeof(numNodesTotal), 1, in);
	}
