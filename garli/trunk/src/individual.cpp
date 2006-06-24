// GARLI version 0.94 source code
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

#include <iosfwd>
#include <iomanip>

using namespace std;

#include "memchk.h"
#include "set.h"
#include "adaptation.h"
#include "model.h"
#include "tree.h"
#include "parameters.h"
#include "population.h"
#include "condlike.h"
#include "mlhky.h"
#include "treenode.h"
#include "individual.h"
#include "funcs.h"
#include "errorexception.h"

extern int memLevel;
extern int calcCount;

#define MUTUALLY_EXCLUSIVE_MUTS

#undef VARIABLE_OPTIMIZATION

//
//
// Methods for class Individual
//
//
Individual::Individual() : dirty(1), fitness(0.0), params(0)
	,reproduced(false), willreproduce(false), parent(-1),
	 willrecombine(false), recombinewith(-1), topo(-1), mutated_brlen(0), 
	 mutation_type(0), accurateSubtrees(0)/*, mutated_rates(0), mutated_pi(0), mutated_alpha(0), mutated_propinvar(0)*/{
	 
 	treeStruct=NULL;
	mod=new Model(6);
}

Individual::~Individual(){
	if(treeStruct!=NULL)
		delete treeStruct;
	if(mod!=NULL) delete mod;
	}

void Individual::CopySecByStealingFirstTree(Individual * sourceOfTreePtr, const Individual *sourceOfInformation){
	CopyNonTreeFields(sourceOfInformation);
	treeStruct=sourceOfTreePtr->treeStruct;
	treeStruct->CopyBranchLens(sourceOfInformation->treeStruct);
	treeStruct->CopyClaIndeces(sourceOfInformation->treeStruct,1);
	dirty=false;
}

void Individual::CopySecByRearrangingNodesOfFirst(Tree * sourceOfTreePtr, const Individual *sourceOfInformation, bool CLAassigned /*=false*/){
	CopyNonTreeFields(sourceOfInformation);
	treeStruct=sourceOfTreePtr;

	for(int i=treeStruct->getNumTipsTotal()+1;i<(2*treeStruct->getNumTipsTotal()-2);i++)
		treeStruct->allNodes[i]->attached=false;
		
	//DZ 10-28 changing this	
	treeStruct->MimicTopo(sourceOfInformation->treeStruct);
	treeStruct->CopyClaIndeces(sourceOfInformation->treeStruct,CLAassigned);
	dirty=false;
	treeStruct->lnL=sourceOfInformation->fitness;
	treeStruct->mod=mod;
}

void Individual::NonSubtreeMutate(const ParallelManager *pMan, double optPrecision, Adaptation *adap)
{//this version is used only by the master when subtree mode is active
//it will make a mutation on one of the nodes that are not contained within
//a subtree, which are in a vector that is passed in

	//because we don't do model mutations during subtree mode, factor the modelMutateProb out
  	double effectiveTopoProb=adap->topoMutateProb / (1.0/(1.0-adap->modelMutateProb));
	double r = rnd.uniform();

#ifndef MUTUALLY_EXCLUSIVE_MUTS
	if(adap->branchOptPrecision != adap->minOptPrecision || r >= effectiveTopoProb){
#else
	if(r >= effectiveTopoProb){
#endif
	 	mutated_brlen=treeStruct->BrlenMutateSubset(pMan->nonSubtreeNodesforSPR);
		if(mutated_brlen > 0){
			mutation_type |= brlen;
			dirty=true;
			}
		}

  if(r < effectiveTopoProb){
	  double r = rnd.uniform();
	  if(r<(adap->randNNIprob/(1.0-adap->randSPRprob)) && (pMan->nonSubtreeNodesforNNI.size() > 0)){
	    int randint1;
	    do{
	    	randint1 = pMan->nonSubtreeNodesforNNI[(int)(pMan->nonSubtreeNodesforNNI.size() *  rnd.uniform())];
	    	}while(randint1<=params->data->NTax());
	    int branch = rnd.uniform() < .5;
	    treeStruct->NNIMutate(randint1,branch,optPrecision, 0);
	    mutation_type |= randNNI;
	    if(treeStruct->lnL !=-1.0){
		    fitness=treeStruct->lnL;
		    dirty=false;
	    	}
   		else dirty=true;
	  }
	  
	  else if(pMan->nonSubtreeNodesforSPR.size() > 3){
	 	int randint1, randint2;
	   	bool done;
	    do{
	    	done=false;
	    	randint1 = pMan->nonSubtreeNodesforSPR[(int)(pMan->nonSubtreeNodesforSPR.size() *  rnd.uniform())];
	    	randint2 = pMan->nonSubtreeNodesforSPR[(int)(pMan->nonSubtreeNodesforSPR.size() *  rnd.uniform())];
	    	//check that the cut node (randint1) is not an ancestor of the attachment node (randint2)
	    	TreeNode *tmp=treeStruct->allNodes[randint2];
	    	while((tmp->nodeNum != 0) && (tmp->nodeNum != randint1)){
	    		tmp=tmp->anc;
	    		}
	    	if(tmp->nodeNum==0) done=true;
	    	
	    	//check if the nodes are siblings 
	    	tmp=treeStruct->allNodes[randint1]->anc;
	    	if(tmp->left->nodeNum==randint2) done=false;
	    	if(tmp->left->next->nodeNum==randint2) done=false;
	    	if(tmp->left->next->next != NULL)
	    		if(tmp->left->next->next->nodeNum == randint2) done=false;
		       	
	    	}while(done == false || treeStruct->allNodes[randint1]->anc->nodeNum==randint2); 
	 
	    treeStruct->SPRMutate(randint1, randint2, optPrecision, pMan->nonSubtreeNodesforNNI);
	    mutation_type |= limSPR;
	    if(treeStruct->lnL !=-1.0){
		    fitness=treeStruct->lnL;
		    dirty=false;
	    	}
   		else dirty=true;
	  	}
	}
/*  
  else{
    assert(TaxonSwapList.size>0);
    double s2, s1 = params->rnd.uniform();
    int randint2, randint1 = TaxonSwapList.size * s1 + 1;
    if(randint1>TaxonSwapList.size) randint1 = TaxonSwapList.size;
    do{
      s2 = params->rnd.uniform();
      randint2 = TaxonSwapList.size * s2 + 1;
    }while(randint2==randint1);

    if(randint2>TaxonSwapList.size) randint2 = TaxonSwapList.size;

    treeStruct->TaxonSwap(randint1, randint2, optPrecision);
    mutation_type |= taxonSwap;
  } 
*/
	CalcFitness(0);

	treeStruct->calcs=calcCount;
	calcCount=0;
}

void Individual::SubtreeMutate(int subdomain, double optPrecision, vector<int> const &subtreeMemberNodes, Adaptation *adap){
  //this version is used only by remotes when they have had a subtree defined for them
  //it will mutate only within that subtree, and because we know that the next mutation
  //will also be within that subtree we can get away without recalculating some likelihood
  //arrays

	//because we don't do model mutations during subtree mode, factor the modelMutateProb out
  	double effectiveTopoProb=adap->topoMutateProb / (1.0/(1.0-adap->modelMutateProb));
	double r = rnd.uniform();
#ifndef MUTUALLY_EXCLUSIVE_MUTS
	if(adap->branchOptPrecision != adap->minOptPrecision || r > effectiveTopoProb){
#else
	if(r >= effectiveTopoProb){
#endif
		mutated_brlen=treeStruct->BrlenMutateSubset( subtreeMemberNodes );
		if(mutated_brlen > 0){
			mutation_type |= brlen;
			dirty=true;
			}
		}

	if(r < effectiveTopoProb){
		r = rnd.uniform();
		int cut;
		if(r<adap->randNNIprob){
		  	//the node passed to the nni function can only be an internal node, so
		  	//pick from the first part of the list which contains the internals
		    cut = subtreeMemberNodes[(int)(rnd.uniform()*(subtreeMemberNodes.size()/2-1))];
		    int branch = rnd.uniform() < .5;
		    treeStruct->NNIMutate(cut, branch, optPrecision, subdomain);
		    mutation_type |= randNNI;
		    if(treeStruct->lnL !=-1.0){
			    fitness=treeStruct->lnL;
			    dirty=false;
		    	}
	   		else dirty=true;
		 	}
	  
		  else if(r < adap->randNNIprob + adap->randSPRprob){
			int broken;
		 	
		 	//the nodes passed to the spr function can be internals or terminals, so 
		 	//choose anywhere in the list
			do{
		    	cut=subtreeMemberNodes[(int)(rnd.uniform()*subtreeMemberNodes.size())];
		    
			    vector<int> SPRList;
				SPRList.reserve(subtreeMemberNodes.size());
			    treeStruct->allNodes[subdomain]->right->getSPRList(cut,SPRList);
			    treeStruct->allNodes[subdomain]->left->getSPRList(cut,SPRList);
			    
			    broken=SPRList[(int)(rnd.uniform()*SPRList.size())];
			    }while(treeStruct->allNodes[broken]->next==treeStruct->allNodes[cut] || 
		    	           treeStruct->allNodes[broken]->prev==treeStruct->allNodes[cut]);
		    	           //reattaching to cut's sib recreates the same tree, so avoid
    
		    treeStruct->SPRMutate(cut, broken, optPrecision, subdomain, 0);
		    mutation_type |= randSPR;
		    if(treeStruct->lnL !=-1.0){
			    fitness=treeStruct->lnL;
			    dirty=false;
			    }
	   		else dirty=true;
		  }
		  else{//limited spr
		 	//the nodes passed to the spr function can be internals or terminals, so 
		 	//choose anywhere in the list
		 	TreeNode *sib;
		 	do{
		    	cut=subtreeMemberNodes[(int)(rnd.uniform()*subtreeMemberNodes.size())];
		    	if(treeStruct->allNodes[cut]->next != NULL) sib=treeStruct->allNodes[cut]->next;
		    	else sib=treeStruct->allNodes[cut]->prev;
		    	}while(treeStruct->allNodes[cut]->anc->nodeNum == subdomain && sib->left==NULL);
		    			
		    treeStruct->SPRMutate(cut, -1, optPrecision, subdomain, adap->limSPRrange);
		    mutation_type |= limSPR;
		    if(treeStruct->lnL !=-1.0){
			    fitness=treeStruct->lnL;
			    dirty=false;
			    }
	   		else dirty=true;
		  }
		}
/*  
  else{
    assert(TaxonSwapList.size>0);
    double s2, s1 = params->rnd.uniform();
    int randint2, randint1 = TaxonSwapList.size * s1 + 1;
    if(randint1>TaxonSwapList.size) randint1 = TaxonSwapList.size;
    do{
      s2 = params->rnd.uniform();
      randint2 = TaxonSwapList.size * s2 + 1;
    }while(randint2==randint1);

    if(randint2>TaxonSwapList.size) randint2 = TaxonSwapList.size;

    treeStruct->TaxonSwap(randint1, randint2, optPrecision);
    mutation_type |= taxonSwap;
  } 
*/
	CalcFitness(subdomain);
	treeStruct->calcs=calcCount;
	calcCount=0;
	}

#ifdef VARIABLE_OPTIMIZATION

#define SPRMutate VariableSPRMutate
#define NNIMutate VariableNNIMutate
#endif

void Individual::Mutate(double optPrecision, Adaptation *adap){
	//this is the original version of mutate, and will be called by both 
	//master and remote when they are mutating a tree that does not have
	//its subtrees properly defined.

	double r = rnd.uniform();
	//DJZ 1-5-05 Moving branch length mutation to be before topo, so that if both are performed
	//the upward sweep needed for blen optimization in the topo mutation will automatically recalc
	//nodes that were dirtied by the blen mutation, and the score of the tree can be finalized at
	//an internal node after the last branch is optimized, rather than waiting until CalcAverageFitness
	//when it will require a sweep down to the root		
#ifndef MUTUALLY_EXCLUSIVE_MUTS
	if(adap->branchOptPrecision != adap->minOptPrecision || r > adap->modelMutateProb + adap->topoMutateProb){
#else
	if(r >= adap->modelMutateProb + adap->topoMutateProb){
#endif
		mutated_brlen=treeStruct->BrlenMutate();
		if(mutated_brlen > 0){
			mutation_type |= brlen;
			dirty=true;
			}
		}

	if(r <= adap->topoMutateProb){
	  r = rnd.uniform();
	  if(r<adap->limSPRprob){
	    treeStruct->SPRMutate(adap->limSPRrange, optPrecision);
	    mutation_type |= limSPR;
	    if(treeStruct->lnL !=-1.0){
		    fitness=treeStruct->lnL;
		    dirty=false;
		    }
   		else dirty=true;
	  }
	  else if (r< adap->randSPRprob + adap->limSPRprob){
	    treeStruct->SPRMutate(0, optPrecision);
	    mutation_type |= randSPR;
	    if(treeStruct->lnL !=-1.0){
		    fitness=treeStruct->lnL;
		    dirty=false;
		    }
		else dirty=true;
	  } 
	  
	  else if (r< adap->taxonSwapprob + adap->randSPRprob + adap->limSPRprob){
	    treeStruct->TaxonSwap(6, optPrecision);
	    mutation_type |= taxonSwap;
	    
	  } 
#ifdef GANESH 
     else if (r < adap->randPECRprob +  adap->taxonSwapprob + adap->randSPRprob + adap->limSPRprob){
        if  (Tree::random_p == true) {
            Tree::p_value = RandomInt(1, 9);
            Tree::ComputeRealCatalan();
        }
		treeStruct->lnL=-1;
        treeStruct->PECRMutate(rnd, optPrecision);
	    mutation_type |= randPECR;
	    if(treeStruct->lnL!=-1){
		    fitness=treeStruct->lnL;
		    dirty=false;
		    }
		else
		    dirty=true;
      }
#endif
	  else {
  		int node=treeStruct->GetRandomInternalNode();
	  	int branch=rnd.uniform() < .5;
	    treeStruct->NNIMutate(node, branch, optPrecision, 0);
	    mutation_type |= randNNI;
   	    if(treeStruct->lnL !=-1.0){
		    fitness=treeStruct->lnL;
		    dirty=false;
		    }
		else dirty=true;
	  }
	} // end if of topomutation
	
	//DJZ making model mutations mutually exclusive with topo mutations
	else if( r < adap->modelMutateProb + adap->topoMutateProb){
	  r = rnd.uniform();
	  double p;
	  if(mod->useFlexRates==true)
		p=1.0/13.0;
	  else if(mod->NoPinvInModel() == false)
	  	p=1.0/11.0;
	  else p=1.0/10.0;

	  if(r<p){
		treeStruct->ScaleWholeTree();
		mutation_type |= muScale;
		}
	  else if(r < p*6.0){
	    mod->MutateRates();
	    mutation_type |= rates;
		}
	  else if(r < p*9.0 && mod->useFlexRates==true){
	    mod->MutateRateMults();
	    mutation_type |= alpha;
	    }
	  else if(r < p*7 && mod->useFlexRates==false){
	    mod->MutateAlpha();
	    mutation_type |= alpha;	
		}
	  else if(r < p*10.0){
	    mod->MutatePis();
	    mutation_type |= pi;
		}
	  else {
		if(mod->useFlexRates==true){
			mod->MutateRateProbs();
			mutation_type |= pinv;
			}
		else{
			mod->MutatePropInvar();
			mutation_type |= pinv;	
			}
		}
	  treeStruct->MakeAllNodesDirty();
	  dirty = true;
	}
	//be sure that we have an accurate score before any CLAs get invalidated
	CalcFitness(0);
	treeStruct->calcs=calcCount;
	calcCount=0;
}

void Individual::CalcFitness(int subtreeNode){
	if(dirty){
		if(subtreeNode>0 && accurateSubtrees==true){
			treeStruct->Score( subtreeNode );
			}
		else treeStruct->Score( );
		
		fitness = treeStruct->lnL;
		dirty = 0;
		}
	
	if(memLevel > 0)
		treeStruct->RemoveTempClaReservations();
	}

void Individual::MakeRandomTree(){
	cout << "creating random starting tree..." << endl;
	treeStruct=new Tree();

	int n = params->data->NTax();
	Set taxset(n);
	for( int i = 1; i <= n; i++ )
		taxset += i;
		
	int placeInAllNodes=n+1;
	
	// add nodes randomly
	for( int i = 0; i < n; i++ ) {
		int pos = rnd.random_int( taxset.Size() );
		int k = taxset[pos];
		treeStruct->AddRandomNode(k, placeInAllNodes  );
		taxset -= k;
		}
	treeStruct->AssignCLAsFromMaster();
	}

void Individual::Randomize(char* fname, int rank){
	assert(params&&params->data);
	
	double ebf[4];
	params->data->CalcEmpiricalFreqs( ebf );
	mod->SetPis(ebf);

	if(mod->NoPinvInModel() || mod->useFlexRates==true){
		mod->SetPinv(0.0);
		mod->SetMaxPinv(0.0);		
		}
	else{
		mod->SetPinv(0.25 * ((double)params->data->NConstant()/(params->data->NConstant()+params->data->NInformative()+params->data->NAutapomorphic())));
		mod->SetMaxPinv((double)params->data->NConstant()/(params->data->NConstant()+params->data->NInformative()+params->data->NAutapomorphic()));
		}
	
	if(strcmp(fname, "random") == 0) MakeRandomTree();

	else { //using a startfile for the initial conditions
	  //12-28-05 This part used to check whether a tree had
	  //previously been read in before going into this loop.
	  //now it goes in regardlesss, since it needs to for
	  // bootstrapping from a starting tree
	  if(1){
	  //if( params->starting_tree.length() == 0 )	{
			if (!FileExists(fname))	{
				throw ErrorException("starting model/tree file \"%s\" does not exist!", fname);
				}
			int ntax=params->data->NTax();
			ifstream stf( fname, ios::in );
			
			if (!stf)	{
				throw ErrorException("starting model/tree file \"%s\" could not be opened!", fname);
				}
			
			//first we need to determine whether there is a model and/or a treestring and
			//check if the taxon numbers or names are present in the tree string
			char c=' ';
			c=stf.get();
			if(c=='#')
					throw ErrorException("Sorry, GARLI does not yet read Nexus tree files.  See manual for starting tree format.");
			int strlen=1;
			bool foundModel=false;
			bool foundTree=false;
			bool numericalTaxa=true;
			while(c!='\n' && c!='\r' && c!=';' && stf.eof()==false){
				if(foundModel==false && foundTree==false){
					if(isalpha(c)){
						if(c=='r'||c=='R'||c=='b'||c=='B'||c=='a'||c=='A'||c=='p'||c=='P') foundModel=true;
						else throw ErrorException("Unknown model parameter specification! \"%c\"", c);
						}
					}
				if(foundTree==false && c=='('){
					foundTree=true;
					}
				if(foundTree==true){
					if(isalpha(c) && c!='e' && c!='E'){//for scientific notation
						numericalTaxa=false;
						}
					}
				strlen++;
				c=stf.get();
				}

			//reopen the file
			stf.close();
			stf.clear();
			stf.open( fname, ios::in );

			char *temp=new char[strlen + 100];

			//if this is a remote population in a parallel run, find the proper tree (ie line number)
			int effectiveRank=rank;
			for(int r=0;r<effectiveRank;r++){
				c=stf.get();
				do{
					c=stf.get();	
					}while(c!='\r' && c!='\n');
				while(stf.peek()=='\r' || stf.peek()=='\n') c=stf.get();
				if(stf.eof() || stf.peek()==EOF){//we hit the end of the file, so we'll just start over.  Figure which tree we want
					effectiveRank=rank%(r+1);
					r=-1; //this is necessary so that when the loop above increments r it will =0
					stf.close();
					stf.clear();
					stf.open( fname, ios::in );
					}
				}
		
			if(foundModel == true){
				c=stf.get();
				do{
					if(c == 'R' || c == 'r'){//rate parameters
						double r[5];
						for(int i=0;i<5;i++){
							stf >> temp;
							r[i]=atof(temp);
							}
						mod->SetRmat(r);
						do{c=stf.get();}while(c==' ');
						if(isdigit(c) || c=='.'){
							stf >> temp;//this is necessary incase GT is included
							c=stf.get();
							}
						foundModel=true;
						}
					else if(c == 'B' || c == 'b'){
						double b[3];
						for(int i=0;i<3;i++){
							stf >> temp;
							b[i]=atof(temp);
							}
						mod->SetPis(b);
						do{c=stf.get();}while(c==' ');
						if(isdigit(c) || c=='.'){
							stf >> temp;
							c=stf.get();							
							}
						foundModel=true;
						}
					else if(c == 'A' || c == 'a'){
#ifdef FLEX_RATES
						assert(0);
#else
						stf >> temp;
						mod->SetAlpha(atof(temp));
						foundModel=true;
						c=stf.get();
#endif
						}				
					else if(c == 'P' || c == 'p'){
#ifdef FLEX_RATES
						assert(0);
#else
						stf >> temp;
						double p=atof(temp);
						if(mod->NoPinvInModel() == true && p > 0.0) throw ErrorException("Error: Value for proportion of invariable sites\nspecified in %s, but not allowed in model\n(see dontinferproportioninvariant in conf file).", fname);
						mod->SetPinv(p);
						foundModel=true;
						c=stf.get();
#endif
						}
					else if(c == 'F' || c == 'f'){
						stf >> temp;
						if(!(isdigit(*temp))) throw ErrorException("Error: expecting number of flex rate categories after \'f\' in starting condition file");
						int n=atoi(temp);
						if(n > 10) throw ErrorException("Error: %d rate categories exceeds maximum number categories (10)", n);
						mod->SetNRateCats(n);
						double rates[10];
						double probs[10];
						for(int i=0;i<n;i++){
							stf >> temp;
							rates[i]=atof(temp);
							stf >> temp;
							probs[i]=atof(temp);
							}		
						mod->SetFlexRates(rates, probs);					
						foundModel=true;
						c=stf.get();						
						}

					else if(isalpha(c)) throw ErrorException("Unknown model parameter specification! \"%c\"", c);
					else if(c != '(') c=stf.get();
					}while(c != '(' && c != '\r' && c != '\n' && !stf.eof());
				if(foundTree == true) stf.putback(c);
				}
			if(foundTree==true){
				stf >> temp;
				params->starting_tree=temp;
				treeStruct=new Tree( params->starting_tree.c_str(), numericalTaxa);
				treeStruct->AssignCLAsFromMaster();
				}
			else MakeRandomTree();

			if(foundTree==true){
				cout << "Obtained starting tree from file " << fname << endl;
				}
			if(foundModel==true){
				cout << "Obtained starting model from file " << fname << endl;
				}
			else{
				cout << "No starting model found in " << fname << endl << "Using default parameter values." << endl;
				}
			mod->OutputGamlFormattedModel(cout);
			cout << endl << endl;

			mod->UpdateQMat();
			stf.close();
			delete []temp;
			}
		assert(treeStruct->root->left->next!=treeStruct->root->right);
		}
	treeStruct->mod=mod;
	treeStruct->CheckBalance();
	dirty=true;
	CalcFitness(0);
	}

void Individual::RefineStartingConditions(bool optModel, double branchPrec){
	if(optModel) cout << "optimizing starting branch lengths and model..." << endl;
	else cout << "optimizing starting branch lengths..." << endl;
	double improve=999.9;
	CalcFitness(0);

/*
ofstream poo("alphatree.tre");
ofstream foo("models.log");
char str[8000000];
mod->OutputPaupBlockForModel(foo, "foo");
treeStruct->root->MakeNewick(str, false);
poo << endl << str << endl;
*/

	for(int i=1;improve > branchPrec;i++){
		double alphaImprove=0.0, pinvImprove=0.0, optImprove=0.0, scaleImprove=0.0;
		
		CalcFitness(0);
		double passStart=Fitness();
		
		optImprove=treeStruct->OptimizeAllBranches(branchPrec);

		SetDirty();
		CalcFitness(0);
		double trueImprove= Fitness() - passStart;
		assert(trueImprove >= 0.0);
		scaleImprove=treeStruct->OptimizeTreeScale();
		SetDirty();
		if(optModel==true){

				//DEBUG
/*						ofstream outf;
						outf.open( "temp.tre" );
						outf.precision(8);
						outf << "#nexus" << endl << endl;
						char treeString[100000];

						//rewritting this to output standard nexus tree files, not gamlviewer stuff
						int ntaxa = params->data->NTax();
						outf << "begin trees;\ntranslate\n";
						for(int k=0;k<ntaxa;k++){
							outf << "  " << (k+1);
							NxsString tnstr = params->data->TaxonLabel(k);
							tnstr.blanks_to_underscores();
							outf << "  " << tnstr.c_str();
							if(k < ntaxa-1) 
								outf << ",\n";
							}		

						outf << ";\n";
						
						outf << "tree best = [&U][" << Fitness() << "][";
						mod->OutputGamlFormattedModel(outf);
						outf << "]";

						outf.setf( ios::floatfield, ios::fixed );
						outf.setf( ios::showpoint );
						treeStruct->root->MakeNewick(treeString, false);
						outf << treeString << ";\n";
						outf << "end;\n";
						
						//add a paup block setting the model params
						mod->OutputPaupBlockForModel(outf, "temp.tre");
						outf.close();
*/				//

			//DEBUG
//			if(mod->useFlexRates==false)
//				alphaImprove=treeStruct->OptimizeAlpha();

			SetDirty();
			}
		improve=scaleImprove + trueImprove + alphaImprove + pinvImprove;
		cout.precision(8);
		cout << "pass " << setw(2) << i << ": +" << setw(7) << improve << "\t(branch=" << trueImprove <<  " scale=" << scaleImprove;
		if(optModel==true) cout << " alpha=" << alphaImprove << ")" << endl;
		else cout << ")" << endl;
		}

	treeStruct->MakeAllNodesDirty();
	treeStruct->nodeOptVector.clear();
	treeStruct->calcs=calcCount;
	calcCount=0;
	dirty=true;
	}

void Individual::ReadTreeFromFile(istream & inf)
{	char tmp[256];
	char ch = ' ';
	NxsString s;

	while( inf )
	{
		inf.get( tmp, 255, '\n' );
		inf.get(ch);
		tmp[255] = '\0';
		s += tmp;
		if( ch == '\n' ) 
			break;
		else
			s += ch;
	}
	treeStruct=new Tree(s.c_str());
	}

void Individual::CopyNonTreeFields(const Individual* ind ){
	fitness = ind->fitness;
	params = ind->params;
	accurateSubtrees=ind->accurateSubtrees;
	mod->CopyModel(ind->mod);
	
	dirty = ind->dirty;
	topo=ind->topo;
	}




