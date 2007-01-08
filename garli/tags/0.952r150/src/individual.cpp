// GARLI version 0.951 source code
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

#include <iosfwd>
#include <iomanip>

using namespace std;

#include "memchk.h"
#include "set.h"
#include "adaptation.h"
#include "model.h"
#include "tree.h"
#include "population.h"
#include "condlike.h"
#include "mlhky.h"
#include "treenode.h"
#include "individual.h"
#include "funcs.h"
#include "errorexception.h"
#include "outputman.h"

extern int memLevel;
extern int calcCount;
extern OutputManager outman;

#define MUTUALLY_EXCLUSIVE_MUTS

#undef VARIABLE_OPTIMIZATION

//
//
// Methods for class Individual
//
//
Individual::Individual() : dirty(1), fitness(0.0), 
	reproduced(false), willreproduce(false), parent(-1),
	willrecombine(false), recombinewith(-1), topo(-1), mutated_brlen(0), 
	mutation_type(0), accurateSubtrees(0){
	 
 	treeStruct=NULL;
	mod=new Model();
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
	    int reconDist = treeStruct->TopologyMutator(optPrecision, adap->limSPRrange, 0);
		if(reconDist == 1 || reconDist == -1) mutation_type |= randNNI;
	    else if(reconDist < 0) mutation_type |= limSPRCon;
		else  mutation_type |= limSPR;
	    if(treeStruct->lnL !=-1.0){
		    fitness=treeStruct->lnL;
		    dirty=false;
		    }
   		else dirty=true;
	  }
	  else if (r< adap->randSPRprob + adap->limSPRprob){
	    int reconDist = treeStruct->TopologyMutator(optPrecision, -1, 0);
		if(reconDist < 0){
			if(reconDist == -1) mutation_type |= randNNI;
			else if(reconDist < -1 * (int)adap->limSPRrange) mutation_type |= randSPRCon;
			else mutation_type |= limSPRCon;
			}
		else {
			if(reconDist == 1) mutation_type |= randNNI;
			else if(reconDist >  (int) adap->limSPRrange) mutation_type |= randSPR;
			else mutation_type |= limSPR;
			}
	    if(treeStruct->lnL !=-1.0){
		    fitness=treeStruct->lnL;
		    dirty=false;
		    }
		else dirty=true;
	  } 
	  else {
		treeStruct->TopologyMutator(optPrecision, 1, 0);
		mutation_type |= randNNI;
   	    if(treeStruct->lnL !=-1.0){
		    fitness=treeStruct->lnL;
		    dirty=false;
		    }
		else dirty=true;
	  }
	} // end if of topomutation
	
	//model mutations
	else if( r < adap->modelMutateProb + adap->topoMutateProb){
		mutation_type |= mod->PerformModelMutation();
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

void Individual::MakeRandomTree(int nTax){
	treeStruct=new Tree();

	int n = nTax;
	Set taxset(n);
	for( int i = 1; i <= n; i++ )
		taxset += i;
		
	int placeInAllNodes=n+1;
	
	if(treeStruct->constraints.empty() == true){
		// add nodes randomly
		for( int i = 0; i < n; i++ ) {
			int pos = rnd.random_int( taxset.Size() );
			int k = taxset[pos];
			treeStruct->AddRandomNode(k, placeInAllNodes  );
			taxset -= k;
			}
		}
	else{
		// add nodes randomly, ensuring that the resulting partial tree is compatible with constraints
		Bipartition mask;
		for( int i = 0; i < n; i++ ) {
			int pos = rnd.random_int( taxset.Size() );
			int k = taxset[pos];
			treeStruct->AddRandomNodeWithConstraints(k, placeInAllNodes, mask );
			taxset -= k;
			}
#ifndef NDEBUG
		treeStruct->CalcBipartitions();
		for(vector<Constraint>::iterator conit=treeStruct->constraints.begin();conit!=treeStruct->constraints.end();conit++){
			if((*conit).IsPositive() == true)
				assert(treeStruct->ContainsBipartitionOrComplement((*conit).GetBipartition()) != NULL);
			else 
				assert(treeStruct->ContainsBipartitionOrComplement((*conit).GetBipartition()) == NULL);
			}
#endif
		}
	treeStruct->AssignCLAsFromMaster();
	}

void Individual::GetStartingConditionsFromFile(const char* fname, int rank, int nTax, bool restart /*=false*/){
	//using a startfile for the initial conditions
	//12-28-05 This part used to check whether a tree had previously been read in before going into
	//this loop. Now it goes in regardlesss, since it needs to for bootstrapping from a starting tree
	
	if (!FileExists(fname))	
		throw ErrorException("starting model/tree file \"%s\" does not exist!", fname);
	ifstream stf( fname, ios::in );
	if (!stf)
		throw ErrorException("starting model/tree file \"%s\" could not be opened!", fname);
	
	bool foundModel, foundTree, numericalTaxa;
	int strlen;
	char c;

	if(restart == false){
		//first we need to determine whether there is a model and/or a treestring and
		//check if the taxon numbers or names are present in the tree string
		c=' ';
		c=stf.get();
		if(c=='#')
				throw ErrorException("Sorry, GARLI does not yet read Nexus tree files.  See manual for starting tree/model format.");
		strlen=1;
		foundModel=false;
		foundTree=false;
		numericalTaxa=true;
		while(c!='\n' && c!='\r' && c!=';' && stf.eof()==false){
			if(foundModel==false && foundTree==false){
				if(isalpha(c)){
					if(c=='r'||c=='R'||c=='b'||c=='B'||c=='a'||c=='A'||c=='p'||c=='P'||c=='n'||c=='f') foundModel=true;
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
		}
	else{//if we are restarting, we can a few things for granted
		//also note the the rank will be incremented by 1, since
		//we want to skip the first line, which had non-tree info on it
		foundModel=foundTree=numericalTaxa=true;
		rank++;
		strlen = (int)((nTax*2)*(10+DEF_PRECISION)+ (double) log10((double) ((double)nTax)*nTax*2));
		}

	//we know what we need to, now reopen the file
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

	//bool foundRmat, foundStateFreqs, foundAlpha, foundPinv;
	//foundRmat=foundStateFreqs=foundAlpha=foundPinv = false;

	if(foundModel == true){//this is UGLY!
		c=stf.get();
		do{
			if(c == 'R' || c == 'r'){//rate parameters
				double r[5];
				for(int i=0;i<5;i++){
					stf >> temp;
					if(temp[0] != '.' && (!isdigit(temp[0]))) throw(ErrorException("Problem reading rate matrix parameters in file %s!\nExamine file and check manual for format.\n", fname));
					r[i]=atof(temp);
					}
				mod->SetRmat(r, restart==false);
				do{c=stf.get();}while(c==' ');
				if(isdigit(c) || c=='.'){
					stf >> temp;//this is necessary incase GT is included
					c=stf.get();
					}
				modSpec.gotRmatFromFile=true;
				}
			else if(c == 'B' || c == 'b'){//base freqs
				double b[3];
				for(int i=0;i<3;i++){
					stf >> temp;
					if(temp[0] != '.' && (!isdigit(temp[0]))) throw(ErrorException("Problem reading base frequency parameters in file %s!\nExamine file and check manual for format.\n", fname));
					b[i]=atof(temp);
					}
				mod->SetPis(b, restart==false);
				do{c=stf.get();}while(c==' ');
				if(isdigit(c) || c=='.'){
					stf >> temp;
					c=stf.get();							
					}
				modSpec.gotStateFreqsFromFile=true;
				}
			else if(c == 'A' || c == 'a'){//alpha shape
				if(modSpec.flexRates==true) throw(ErrorException("Config file specifies ratehetmodel = flex, but starting model contains alpha!\n"));
				stf >> temp;
				if(temp[0] != '.' && (!isdigit(temp[0]))) throw(ErrorException("Problem reading alpha parameter in file %s!\nExamine file and check manual for format.\n", fname));
				mod->SetAlpha(atof(temp), restart==false);
				c=stf.get();
				modSpec.gotAlphaFromFile=true;
				}				
			else if(c == 'P' || c == 'p'){//proportion invariant
				stf >> temp;
				if(temp[0] != '.' && (!isdigit(temp[0]))) throw(ErrorException("Problem reading proportion of invariant sites parameter in file %s!\nExamine file and check manual for format.\n", fname));
				double p=atof(temp);
				mod->SetPinv(p, restart==false);
				c=stf.get();
				modSpec.gotPinvFromFile=true;
				}
			else if(c == 'F' || c == 'f'){//flex rates
				double rates[20];
				double probs[20];
				for(int i=0;i<mod->NRateCats();i++){
					stf >> temp;
					if(isalpha(temp[0])) throw ErrorException("Problem with flex rates specification in starting condition file");
					rates[i]=atof(temp);
					stf >> temp;
					if(isalpha(temp[0])) throw ErrorException("Problem with flex rates specification in starting condition file");
					probs[i]=atof(temp);
					}		
				mod->SetFlexRates(rates, probs);					
				c=stf.get();
				modSpec.gotFlexFromFile=true;
				}
			else if(c == 'n'){
				//the number of cats should now be set in the config file
				c=stf.get();
				assert(0);
				}
			else if(isalpha(c)) throw(ErrorException("Unknown model parameter specification in file %s!\nExamine file and check manual for format.\n", fname));
			else if(c != '(') c=stf.get();
			}while(c != '(' && c != '\r' && c != '\n' && !stf.eof());

		if(foundTree == true) stf.putback(c);
		}//if(foundModel == true)

	//Here we'll error out if something was fixed but didn't appear
	if(modSpec.fixRelativeRates == true && modSpec.gotRmatFromFile == false) throw ErrorException("ratematrix = fixed in conf file, but parameter values not found in %s.", fname);
	if(modSpec.fixStateFreqs == true && modSpec.equalStateFreqs == false && modSpec.empiricalStateFreqs == false && modSpec.gotStateFreqsFromFile == false) throw ErrorException("statefrequencies = fixed in conf file, but parameter values not found in %s.", fname);
	if(modSpec.fixAlpha == true && modSpec.gotAlphaFromFile == false) throw ErrorException("ratehetmodel = gammafixed in conf file, but no parameter value for alpha found in %s.", fname);
	if(modSpec.fixInvariantSites == true && modSpec.gotPinvFromFile == false) throw ErrorException("invariantsites = fixed in conf file, but no parameter value found in %s.", fname);

	if(foundTree==true){
		stf >> temp;
		treeStruct=new Tree(temp, numericalTaxa);

		//check that any defined constraints are present in the starting tree
		treeStruct->CalcBipartitions();
		int conNum=1;
		for(vector<Constraint>::iterator conit=treeStruct->constraints.begin();conit!=treeStruct->constraints.end();conit++){
			if((*conit).IsPositive()){
				if(treeStruct->ContainsBipartitionOrComplement((*conit).GetBipartition()) == NULL) throw ErrorException("Starting tree not compatible with constraint number %d!!!", conNum);
			}
			else if(treeStruct->ContainsBipartitionOrComplement((*conit).GetBipartition()) != NULL) throw ErrorException("Starting tree not compatible with constraint number %d!!!", conNum);
			conNum++;
			}
		treeStruct->AssignCLAsFromMaster();
		}

	else MakeRandomTree(nTax);

	if(restart == false){
		if(foundTree==true)
			outman.UserMessage("Obtained starting tree from file %s", fname);
		else 
			outman.UserMessage("No starting tree found in file %s, creating random tree", fname);

		if(foundModel==true) outman.UserMessage("Obtained starting model from file %s:", fname);
				
		else outman.UserMessage("No starting model found in %s\nUsing default parameter values:", fname);
			
		mod->OutputGarliFormattedModel(cout);
		outman.UserMessage("\n");
		}

	mod->UpdateQMat();
	stf.close();
	delete []temp;
	}

void Individual::RefineStartingConditions(bool optModel, double branchPrec){
	if(optModel && mod->NRateCats() > 1 && modSpec.gotFlexFromFile == false) outman.UserMessage("optimizing starting branch lengths and rate heterogeneity parameters...");
	else outman.UserMessage("optimizing starting branch lengths...");
	double improve=999.9;
	CalcFitness(0);

	for(int i=1;improve > branchPrec;i++){
		double alphaImprove=0.0, optImprove=0.0, scaleImprove=0.0;
		
		CalcFitness(0);
		double passStart=Fitness();
		
		optImprove=treeStruct->OptimizeAllBranches(branchPrec);

		SetDirty();
		CalcFitness(0);
		double trueImprove= Fitness() - passStart;
//		assert(trueImprove >= -1.0);
		if(trueImprove < 0.0) trueImprove = 0.0;
		scaleImprove=treeStruct->OptimizeTreeScale(branchPrec);
		SetDirty();
		if(optModel==true && mod->NRateCats() > 1 && modSpec.fixAlpha == false && modSpec.gotFlexFromFile == false){
			alphaImprove=treeStruct->OptimizeAlpha(branchPrec);
			SetDirty();
			}
		improve=scaleImprove + trueImprove + alphaImprove;
		outman.precision(8);
		if(optModel==true && mod->NRateCats() > 1){
			if(modSpec.flexRates == false) outman.UserMessage("pass %-2d: +%10.4f (branch=%8.2f scale=%8.2f alpha=%8.2f)", i, improve, trueImprove, scaleImprove, alphaImprove);
			else if(modSpec.gotFlexFromFile == false) outman.UserMessage("pass %-2d: +%10.4f (branch=%8.2f scale=%8.2f flex rates=%8.2f)", i, improve, trueImprove, scaleImprove, alphaImprove);
			else outman.UserMessage("pass %-2d: +%10.4f (branch=%8.2f scale=%8.2f)", i, improve, trueImprove, scaleImprove);
			}
		else outman.UserMessage("pass %-2d: +%10.4f (branch=%8.2f scale=%8.2f)", i, improve, trueImprove, scaleImprove);
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
	treeStruct=new Tree(s.c_str(), true);
	}

void Individual::CopyNonTreeFields(const Individual* ind ){
	fitness = ind->fitness;
	accurateSubtrees=ind->accurateSubtrees;
	mod->CopyModel(ind->mod);
	
	dirty = ind->dirty;
	topo=ind->topo;
	}




/* 7/21/06 needs to be fixed to correspond to changes in tree for constraints
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
*//*
	CalcFitness(subdomain);
	treeStruct->calcs=calcCount;
	calcCount=0;
	}
*/

/*7/21/06 needs to be fixed to correspond to changes in tree for constraints
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
*/ /*  
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
/*	CalcFitness(0);

	treeStruct->calcs=calcCount;
	calcCount=0;
}
*/