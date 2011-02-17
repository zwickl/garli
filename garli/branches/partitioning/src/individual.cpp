// GARLI version 0.96b8 source code
// Copyright 2005-2008 Derrick J. Zwickl
// email: zwickl@nescent.org
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

#include <iosfwd>
#include <iomanip>
#include <sstream>

using namespace std;

#include "defs.h"
#include "set.h"
#include "funcs.h"
#include "adaptation.h"
#include "model.h"
#include "tree.h"
#include "population.h"
#include "condlike.h"
#include "sequencedata.h"
#include "treenode.h"
#include "individual.h"
#include "errorexception.h"
#include "outputman.h"
#include "reconnode.h"
#include "utility.h"

extern int memLevel;
extern int calcCount;
extern OutputManager outman;
extern FLOAT_TYPE globalBest;

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
//	mod=new Model();
	}

Individual::Individual(const Individual *other) : 
	dirty(1), fitness(0.0), 
	reproduced(false), willreproduce(false), parent(-1),
	willrecombine(false), recombinewith(-1), topo(-1), mutated_brlen(0), 
	mutation_type(0), accurateSubtrees(0){

	//mod=new Model();
	treeStruct=new Tree();

	CopyNonTreeFields(other);
	
	treeStruct->MimicTopo(other->treeStruct);
	dirty=false;
	treeStruct->lnL=other->fitness;
	treeStruct->modPart = &modPart;
	}

Individual::~Individual(){
	if(treeStruct!=NULL)
		delete treeStruct;
	//if(mod!=NULL) delete mod;
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
	modPart.CopyModelPartition(&sourceOfInformation->modPart);
	treeStruct->modPart = &modPart;
	}

void Individual::DuplicateIndivWithoutCLAs(const Individual *sourceOfInformation){
	CopyNonTreeFields(sourceOfInformation);
	if(treeStruct == NULL)
		treeStruct = new Tree;

	for(int i=treeStruct->getNumTipsTotal()+1;i<(2*treeStruct->getNumTipsTotal()-2);i++)
		treeStruct->allNodes[i]->attached=false;
		
	treeStruct->MimicTopo(sourceOfInformation->treeStruct);
	dirty=true;
	treeStruct->lnL=sourceOfInformation->fitness;
	modPart.CopyModelPartition(&sourceOfInformation->modPart);
	treeStruct->modPart = &modPart;
	}

void Individual::Mutate(FLOAT_TYPE optPrecision, Adaptation *adap){
	//this is the original version of mutate, and will be called by both 
	//master and remote when they are mutating a tree that does not have
	//its subtrees properly defined.

	FLOAT_TYPE r = rnd.uniform();
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

	try{
		if(r <= adap->topoMutateProb){
		  r = rnd.uniform();
		  if(r<adap->limSPRprob){
			int reconDist = treeStruct->TopologyMutator(optPrecision, adap->limSPRrange, 0);
			if(reconDist == 1 || reconDist == -1) mutation_type |= randNNI;
			else if(reconDist < 0) mutation_type |= limSPRCon;
			else  mutation_type |= limSPR;
			if(!FloatingPointEquals(treeStruct->lnL, -ONE_POINT_ZERO, 1.0e-8)){
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
			if(!FloatingPointEquals(treeStruct->lnL, -ONE_POINT_ZERO, 1.0e-8)){
				fitness=treeStruct->lnL;
				dirty=false;
				}
			else dirty=true;
		  } 
		  else {
			treeStruct->TopologyMutator(optPrecision, 1, 0);
			mutation_type |= randNNI;
   			if(!FloatingPointEquals(treeStruct->lnL, -ONE_POINT_ZERO, 1.0e-8)){
				fitness=treeStruct->lnL;
				dirty=false;
				}
			else dirty=true;
		  }
		} // end if of topomutation
		
		//model mutations
		else if( r < adap->modelMutateProb + adap->topoMutateProb){
			mutation_type |= modPart.PerformModelMutation();
			treeStruct->MakeAllNodesDirty();
			dirty = true;
			}

		//be sure that we have an accurate score before any CLAs get invalidated
		CalcFitness(0);
		}
	catch(UnscoreableException &ex){
		//in some situations the tree just underflows no matter what - I've only seen this and only
		//throw this from orientedGap models with very poor trees.
		outman.DebugMessage("WARNING - created individual deemed unscorable!");
		treeStruct->lnL = -FLT_MAX;
		SetFitness(-FLT_MAX);
		}

/*	FLOAT_TYPE lnL = fitness;
	dirty = true;
	treeStruct->MakeAllNodesDirty();
	CalcFitness(0);
	if(!FloatingPointEquals(lnL, fitness, 1e-3)){
		outman.UserMessage("DEBUG - scoring problem:%f vs %f", lnL, fitness);
		//throw ErrorException("DEBUG - scoring problem:%f vs %f", lnL, fitness);
		}
*/
//	treeStruct->calcs=calcCount;
//	calcCount=0;
}

void Individual::CalcFitness(int subtreeNode){
	if(dirty || FloatingPointEquals(treeStruct->lnL, ZERO_POINT_ZERO, max(1.0e-8, GARLI_FP_EPS * 2.0)) || FloatingPointEquals(treeStruct->lnL, -ONE_POINT_ZERO, max(1.0e-8, GARLI_FP_EPS * 2))){
		if(subtreeNode>0 && accurateSubtrees==true){
			treeStruct->Score( subtreeNode );
			}
		else treeStruct->Score( );
		
		fitness = treeStruct->lnL;
		dirty = 0;
		}
	else{
		assert(!FloatingPointEquals(treeStruct->lnL, -ONE_POINT_ZERO, 1e-8));
		fitness = treeStruct->lnL;
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
			treeStruct->RandomlyAttachTip(k, placeInAllNodes  );
			taxset -= k;
			}
		}
	else{
		// add nodes randomly, ensuring that the resulting partial tree is compatible with constraints
		Bipartition mask;
		for( int i = 0; i < n; i++ ) {
			int pos = rnd.random_int( taxset.Size() );
			int k = taxset[pos];
			treeStruct->RandomlyAttachTipWithConstraints(k, placeInAllNodes, &mask );
			taxset -= k;
			}
#ifndef NDEBUG
		for(vector<Constraint>::iterator conit=treeStruct->constraints.begin();conit!=treeStruct->constraints.end();conit++){
			TreeNode *check = NULL;
			if((*conit).IsBackbone())
				check = treeStruct->ContainsMaskedBipartitionOrComplement(*(*conit).GetBipartition(), *(*conit).GetBackboneMask());
			else
				check = treeStruct->ContainsBipartitionOrComplement(*(*conit).GetBipartition());
			if((*conit).IsPositive()) assert(check != NULL);
			else assert(check == NULL);
			}
#endif
		}
	if(treeStruct->dummyRootBranchMidpoint)
		treeStruct->MoveDummyRootToBranchMidpoint();

	treeStruct->AssignCLAsFromMaster();
	}

void Individual::MakeStepwiseTree(int nTax, int attachesPerTaxon, FLOAT_TYPE optPrecision ){
	treeStruct=new Tree();
	treeStruct->modPart = &modPart;
	treeStruct->AssignCLAsFromMaster();

	Individual scratchI;
	scratchI.treeStruct=new Tree();
	Tree *scratchT = scratchI.treeStruct;
	scratchT->modPart = &scratchI.modPart;
	scratchT->AssignCLAsFromMaster();
	scratchI.CopySecByRearrangingNodesOfFirst(scratchT, this, true);

	int n = nTax;
	Set taxset(n);
	for( int i = 1; i <= n; i++ )
		taxset += i;
		
	int placeInAllNodes=n+1;
//	ofstream stepout("stepwise.log");
	outman.UserMessage("number of taxa added:");

	Bipartition mask;//mask is used for constrained trees
	for(int i = 0;i<3;i++){//add the first 3
		int pos = rnd.random_int( taxset.Size() );
		int k = taxset[pos];
		if(treeStruct->constraints.empty())
			scratchT->RandomlyAttachTip(k, placeInAllNodes  );
		else
			scratchT->RandomlyAttachTipWithConstraints(k, placeInAllNodes, &mask );
		taxset -= k;
		}
	//use information on the similarity between sequences to choose first stepwise additions
/*	
	const SequenceData *dat = treeStruct->data;
	int nstates = mod->NStates();
	FLOAT_TYPE **pdist = New2DArray<FLOAT_TYPE>(dat->NTax(), dat->NTax());
	for(int i=0;i<nTax;i++){
		pdist[i][i] = 0.0;
		for(int j=i+1;j<nTax;j++){
			pdist[i][j] = CalculateHammingDistance((char*) dat->GetRow(i), (char*) dat->GetRow(j), dat->GetCounts(), dat->NChar(), nstates);
			pdist[j][i] = pdist[i][j];
			}
		}
	//add the first 3
	//be careful because the taxa are indexed from 1->ntax
	int pos = rnd.random_int( taxset.Size() );
	int first = (taxset[pos]);
	scratchT->RandomlyAttachTip(first, placeInAllNodes  );
	taxset -= first;
	
	//add the furthest taxon to that
	int sec = 1;
	FLOAT_TYPE maxDist = pdist[first-1][sec-1];
	for(int i=sec+1;i<=dat->NTax();i++){
		if(pdist[first-1][i-1] > maxDist){
			sec = i; 
			maxDist = pdist[first-1][sec-1];
			}
		}
	scratchT->RandomlyAttachTip(sec, placeInAllNodes  );
	taxset -= sec;
	//add the furthest taxon to that (which may in fact be close to first, but should not have a pdist = 0 to it)
	int third = (first == 1 ? 2 : 1);
	maxDist = pdist[sec-1][third-1];
	for(int i=third+1;i<=dat->NTax();i++){
		if(pdist[sec-1][i] > maxDist && i != first && pdist[first-1][third-1] > ZERO_POINT_ZERO){
			third = i; 
			maxDist = pdist[sec-1][third-1];
			}
		}
	scratchT->RandomlyAttachTip(third, placeInAllNodes  );
	taxset -= third;
*/
	CopySecByRearrangingNodesOfFirst(treeStruct, &scratchI, true);

	for( int i = 3; i < n; i++ ) {
		//select a random node
		int pos = rnd.random_int( taxset.Size() );
		int k = taxset[pos];
		taxset -= k;
		//add the node randomly - this is a little odd, but for the existing swap collecting machinery
		//to work right, the taxon to be added needs to already be in the tree
		if(treeStruct->constraints.empty())
			scratchT->RandomlyAttachTip(k, placeInAllNodes  );
		else
			scratchT->RandomlyAttachTipWithConstraints(k, placeInAllNodes, &mask );
		TreeNode *added = scratchT->allNodes[k];

		scratchT->SweepDirtynessOverTree(added);
		scratchT->OptimizeBranchesWithinRadius(added->anc, optPrecision, 0, NULL);

		//backup what we have now
		CopySecByRearrangingNodesOfFirst(treeStruct, &scratchI, true);
		FLOAT_TYPE bestScore = scratchT->lnL;
		
		//collect reconnection points - this will automatically filter for constraints
		scratchT->GatherValidReconnectionNodes(scratchT->NTax()*2, added, NULL, &mask);
		
//			stepout << i << "\t" << k << "\t" << bestScore << "\t";

		//start swappin
		int num=0;
		//for(list<ReconNode>::iterator b = scratchT->sprRang.begin();b != scratchT->sprRang.end();b++){
		ReconList attempted;
		while(num < attachesPerTaxon && scratchT->sprRang.size() > 0){
			int connectNum = rnd.random_int(scratchT->sprRang.size());
			listIt broken = scratchT->sprRang.NthElement(connectNum);
			//try a reattachment point
			scratchT->SPRMutate(added->nodeNum, &(*broken), optPrecision, 0);
			//record the score
			broken->chooseProb = scratchT->lnL;
			attempted.AddNode(*broken);
			scratchT->sprRang.RemoveNthElement(connectNum);
//			stepout << scratchT->lnL << "\t";
			//restore the tree
			scratchI.CopySecByRearrangingNodesOfFirst(scratchT, this, true);
			num++;
			}
		//now find the best score
		ReconNode *best = NULL;
		
		//For debugging, add to random place, to check correct filtering of attachment points for constraints
/*
		if(attempted.size() != 0)
			best = attempted.RandomReconNode();
*/
		for(list<ReconNode>::iterator b = attempted.begin();b != attempted.end();b++){
			if((*b).chooseProb > bestScore){
				best = &(*b);
				bestScore = (*b).chooseProb;
				}
			}

		//if we didn't find anything better than the initial random attachment we don't need to do anything
		if(best != NULL){
			scratchT->SPRMutate(added->nodeNum, best, optPrecision, 0);
			}
		else scratchT->Score();
		scratchI.CalcFitness(0);

//		stepout << scratchT->lnL << endl;
		CopySecByRearrangingNodesOfFirst(treeStruct, &scratchI, true);

		//outman.UserMessage(" %d %f", i+1, scratchT->lnL);
		outman.UserMessageNoCR(" %d ", i+1);
		outman.flush();
		//when we've added half the taxa optimize alpha, flex or omega 
		if(i == (n/2)){
			FLOAT_TYPE improve = 0.0;
			for(int modnum = 0;modnum < modPart.NumModels();modnum++){
				Model *mod = scratchI.modPart.GetModel(modnum);
				const ModelSpecification *modSpec = mod->GetCorrespondingSpec();
				if(modSpec->IsCodon())//optimize omega even if there is only 1
					improve += scratchT->OptimizeOmegaParameters(optPrecision, modnum);
				else if(mod->NRateCats() > 1){
					if(modSpec->IsFlexRateHet()){//Flex rates
						//no longer doing alpha first, it was too hard to know if the flex rates had been partially optimized
						//already during making of a stepwise tree
						improve += scratchT->OptimizeFlexRates(optPrecision, modnum);
						}
					else if(modSpec->fixAlpha == false){//normal gamma
						//do NOT let alpha go too low here - on bad or random starting trees the branch lengths get crazy long
						improve += scratchT->OptimizeBoundedParameter(modnum, optPrecision, mod->Alpha(), 0, 0.05, 999.9, &Model::SetAlpha);
						}
					}
				if(modSpec->includeInvariantSites && !modSpec->fixInvariantSites)
					improve += scratchT->OptimizeBoundedParameter(modnum, optPrecision, mod->PropInvar(), 0, 1.0e-8, mod->maxPropInvar, &Model::SetPinv);
				}
			if(modSpecSet.InferSubsetRates()){
				improve += scratchT->OptimizeSubsetRates(optPrecision);
				}
			outman.UserMessageNoCR("\nOptimizing parameters... improved %.3f lnL", improve);
			scratchT->Score();
			FLOAT_TYPE start=scratchT->lnL;
			scratchT->OptimizeAllBranches(optPrecision);
			FLOAT_TYPE bimprove = max(scratchT->lnL - start, 0.0);
			outman.UserMessage("\nOptimizing branchlengths... improved %.3f lnL", bimprove);
			}
		}		

//	stepout.close();
	outman.UserMessage("");
	scratchI.treeStruct->RemoveTreeFromAllClas();
	delete scratchI.treeStruct;
	scratchI.treeStruct=NULL;
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
		if(c=='#'){//nexus tree files should now be going through NCL elsewhere, so we shouldn't be here
				assert(0);
				throw ErrorException("Sorry, GARLI does not yet read Nexus tree files.  See manual for starting tree/model format.");
				}
		strlen=1;
		foundModel=false;
		foundTree=false;
		numericalTaxa=true;
		while(c!='\n' && c!='\r' && c!=';' && stf.eof()==false){
			if(foundModel==false && foundTree==false){
				if(isalpha(c)){
					//changing from b for base freqs to e, for equilibrium freqs
					if(c=='r'||c=='R'||c=='b'||c=='B'||c=='e'||c=='E'||c=='a'||c=='A'||c=='p'||c=='P'||c=='i'||c=='I'||c=='f'||c=='o'||c=='O'||c=='M'||c=='m'||c=='S'||c=='s') 
						foundModel=true;
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
		assert(0);
		foundModel=foundTree=numericalTaxa=true;
		rank++;
		strlen = (int)((nTax*2)*(10+DEF_PRECISION)+ (FLOAT_TYPE) log10((FLOAT_TYPE) ((FLOAT_TYPE)nTax)*nTax*2));
		}

	//we know what we need to, now reopen the file
	stf.close();
	stf.clear();
	stf.open( fname, ios::in );

	char *temp=new char[strlen + 100];
	
	//if this is a remote population in a parallel run or a multirep run, find the proper tree (ie line number)
	int effectiveRank=rank;
	for(int r=0;r<effectiveRank;r++){
		c=stf.get();
		do{
			c=stf.get();	
			}while(c!='\r' && c!='\n' && !stf.eof());
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
	
	if(foundModel == true){
//		if(modPart.NumModels() > 1)
//			throw ErrorException("Specification of model parameter values is not yet supported with partitioned models");
		string modString;
		do{
			c=stf.get();
			modString += c;
			//}while(c != '(' && c != '\r' && c != '\n' && !stf.eof());
			}while(stf.peek() != '(' && stf.peek() != '\r' && stf.peek() != '\n' && !stf.eof());
		while((stf.peek() == '\n' || stf.peek() == '\r') && stf.eof() == false) 
			stf.get(c);
		modPart.ReadGarliFormattedModelStrings(modString);
		}

	if(foundTree==true){
		string treeString;
		char c;
		stf.get(c);
		do{
			treeString += c;
			stf.get(c);
			}while(c != '\n' && c!= '\r' && stf.eof() == false);
		while((stf.peek() == '\n' || stf.peek() == '\r') && stf.eof() == false) stf.get(c);

		//the call to the tree constructor can change the seed because random branch lengths are generated when the tree doesn't
		//have them.  So, store and restore the seed, mainly for output purposes (the seed output to the screen happens after this
		//call
		int seed = rnd.seed();

		//now allowing polytomies, since they will be taken care of in Population::SeedPopulationWithStartingTree
		treeStruct=new Tree(treeString.c_str(), numericalTaxa, true);
		//treeStruct=new Tree(treeString.c_str(), numericalTaxa);

		//check that any defined constraints are present in the starting tree
		int conNum=1;
		for(vector<Constraint>::iterator conit=treeStruct->constraints.begin();conit!=treeStruct->constraints.end();conit++){
			TreeNode *check = NULL;
			if((*conit).IsBackbone())
				check = treeStruct->ContainsMaskedBipartitionOrComplement(*(*conit).GetBipartition(), *(*conit).GetBackboneMask());
			else
				check = treeStruct->ContainsBipartitionOrComplement(*(*conit).GetBipartition());
			if(((*conit).IsPositive() && check == NULL) || ((*conit).IsPositive() == false  && check != NULL))
				throw ErrorException("Starting tree not compatible with constraint number %d!!!", conNum);
			}
		treeStruct->AssignCLAsFromMaster();
		}

	//if no tree is found the making of the random tree will now be taken care of back in Population::SeedPopulationWithStartingTree
	//else MakeRandomTree(nTax);

	if(restart == false){
		if(!foundTree && !foundModel) 
			throw ErrorException("No starting tree or model was found in the specified starting conditions\n\tfile %s.\n\tIf it is a Nexus file it must start with #NEXUS\n\tOtherwise see manual for information on starting condition format.", fname);

		if(foundTree==true)
			outman.UserMessage("Obtained starting tree %d from file %s",  effectiveRank+1, fname);
		else{
			outman.UserMessage("No starting tree found in file %s", fname);
			}

		if(foundModel==true){
			outman.UserMessage("Obtained starting or fixed model parameter values from file %s", fname);
			string m;
			modPart.FillGarliFormattedModelStrings(m);
			outman.UserMessage("%s", m.c_str());
			}
		else{
			//this checks whether we have already gotten some parameter values from file, which might have come from a garli block in the datafile
			if(!(modSpecSet.GotAnyParametersFromFile())){
				outman.UserMessage("No starting model parameter values found in %s\nUsing default parameter values", fname);
				}
			}
			
		outman.UserMessage("");
		}

	for(int m=0;m < modPart.NumModels();m++){
		modPart.GetModel(m)->UpdateQMat();
		}
	stf.close();
	delete []temp;
	}

void Individual::GetStartingTreeFromNCL(const NxsTreesBlock *treesblock, int rank, int nTax, bool restart /*=false*/){
	assert(treeStruct == NULL);

	int totalTrees = treesblock->GetNumTrees();

	int effectiveRank = rank % totalTrees;
	
	//the call to the tree constructor can change the seed because random branch lengths are generated when the tree doesn't
	//have them.  So, store and restore the seed, mainly for output purposes (the seed output to the screen happens after this
	//call
	int seed = rnd.seed();

	//we will get the tree string from NCL with taxon numbers (starting at 1), regardless of how it was initially read in 
	const NxsFullTreeDescription &t = treesblock->GetFullTreeDescription(effectiveRank);
	if(t.AllTaxaAreIncluded() == false && !treeStruct->someOrientedGap)
		throw ErrorException("Starting tree description must contain all taxa.");
	string ts = t.GetNewick();
	ts += ";";
	treeStruct=new Tree(ts.c_str(), true, true);

	rnd.set_seed(seed);

	//check that any defined constraints are present in the starting tree
	int conNum=1;
	for(vector<Constraint>::iterator conit=treeStruct->constraints.begin();conit!=treeStruct->constraints.end();conit++){
		TreeNode *check = NULL;
		if((*conit).IsBackbone())
			check = treeStruct->ContainsMaskedBipartitionOrComplement(*(*conit).GetBipartition(), *(*conit).GetBackboneMask());
		else
			check = treeStruct->ContainsBipartitionOrComplement(*(*conit).GetBipartition());
		if(((*conit).IsPositive() && check == NULL) || ((*conit).IsPositive() == false  && check != NULL))
			throw ErrorException("Starting tree not compatible with constraint number %d!!!", conNum);
		}
	treeStruct->AssignCLAsFromMaster();

	for(int m=0;m < modPart.NumModels();m++){
		modPart.GetModel(m)->UpdateQMat();
		}
	}

void Individual::RefineStartingConditions(bool optModel, FLOAT_TYPE branchPrec){
	bool optOmega, optAlpha, optFlex, optPinv, optFreqs, optRelRates, optSubsetRates;
	optOmega = optAlpha = optFlex = optPinv = optFreqs = optRelRates = optSubsetRates = false;

	bool optInsDel = false;

	if(optModel){
		for(int modnum = 0;modnum < modPart.NumModels();modnum++){
			Model *mod = modPart.GetModel(modnum);
			const ModelSpecification *modSpec = mod->GetCorrespondingSpec();
			if(modSpec->numRateCats > 1 && modSpec->IsNonsynonymousRateHet() == false && modSpec->IsFlexRateHet() == false) optAlpha = true;
			if(modSpec->IsFlexRateHet()) optFlex = true;
			if(modSpec->includeInvariantSites && modSpec->fixInvariantSites == false) optPinv = true;
			if(modSpec->IsCodon() && !modSpec->fixOmega) optOmega = true;
			if(modSpec->IsOrientedGap()) optInsDel = true;

			if(modSpec->IsCodon() == false && modSpec->fixStateFreqs == false && modSpec->IsEqualStateFrequencies() == false && modSpec->IsEmpiricalStateFrequencies() == false)
				optFreqs = true;
			//this is the case of forced freq optimization with codon models.  For everything to work they must be set as both not fixed but empirical
			if(modSpec->IsCodon() && modSpec->fixStateFreqs == false && modSpec->IsEqualStateFrequencies() == false && modSpec->IsEmpiricalStateFrequencies() == true)
				optFreqs = true;
			if(modSpec->fixRelativeRates == false && (modSpec->Nst() > 1 || modSpec->IsEstimateAAMatrix() || modSpec->IsTwoSerineRateMatrix()))
				optRelRates = true;
			}
		if(modSpecSet.InferSubsetRates() && modSpecSet.NumSpecs() > 1)
			optSubsetRates = true;
		}

	outman.UserMessageNoCR("optimizing: starting branch lengths");
	if(optAlpha) outman.UserMessageNoCR(", alpha shape");
	if(optPinv) outman.UserMessageNoCR(", prop. invar");
	if(optRelRates) outman.UserMessageNoCR(", rel rates");
	if(optFreqs) outman.UserMessageNoCR(", eq freqs");
	if(optOmega) outman.UserMessageNoCR(", dN/dS (aka omega) parameters");
	if(optInsDel){
		outman.UserMessageNoCR(", ins rate");
		outman.UserMessageNoCR(", del rate");
		}
	if(optSubsetRates) outman.UserMessageNoCR(", subset rates");
	outman.UserMessage("...");
	FLOAT_TYPE improve=(FLOAT_TYPE)999.9;
	CalcFitness(0);

	for(int i=1;improve > branchPrec;i++){
		FLOAT_TYPE alphaOptImprove=0.0, pinvOptImprove = 0.0, omegaOptImprove = 0.0, flexOptImprove = 0.0, optImprove=0.0, scaleOptImprove=0.0, subsetRateImprove=0.0, rateOptImprove=0.0;
		FLOAT_TYPE freqOptImprove=0.0, insDelImprove = 0.0;
		
		CalcFitness(0);
		FLOAT_TYPE passStart=Fitness();
		
		optImprove=treeStruct->OptimizeAllBranches(branchPrec);
		CalcFitness(0);

		FLOAT_TYPE trueImprove= Fitness() - passStart;
		assert(trueImprove >= -1.0);
		if(trueImprove < ZERO_POINT_ZERO) trueImprove = ZERO_POINT_ZERO;

		vector<FLOAT_TYPE> blens;
		treeStruct->StoreBranchlengths(blens);
		scaleOptImprove=treeStruct->OptimizeTreeScale(branchPrec);
		CalcFitness(0);
		//if some of the branch lengths were at the minimum or maximum boundaries the scale optimization
		//can actually worsen the score.  If so, return them to their original lengths.
		if(scaleOptImprove < ZERO_POINT_ZERO){
			treeStruct->RestoreBranchlengths(blens);
			CalcFitness(0);
			scaleOptImprove = ZERO_POINT_ZERO;
			}

		CalcFitness(0);
		if(optModel){
			for(int modnum = 0;modnum < modPart.NumModels();modnum++){
				Model *mod = modPart.GetModel(modnum);
				const ModelSpecification *modSpec = mod->GetCorrespondingSpec();
				if(modSpec->IsCodon() && !modSpec->fixOmega)//optimize omega even if there is only 1
					omegaOptImprove += treeStruct->OptimizeOmegaParameters(branchPrec, modnum);
				else if(mod->NRateCats() > 1){
					if(modSpec->IsFlexRateHet()){//Flex rates
						//no longer doing alpha first, it was too hard to know if the flex rates had been partially optimized
						//already during making of a stepwise tree
						//if(i == 1) rateOptImprove = treeStruct->OptimizeAlpha(branchPrec);
						//if(i == 1 && modSpec.gotFlexFromFile==false) rateOptImprove = treeStruct->OptimizeBoundedParameter(branchPrec, mod->Alpha(), 0, 1.0e-8, 999.9, &Model::SetAlpha);
						flexOptImprove += treeStruct->OptimizeFlexRates(branchPrec, modnum);
						}
					else if(modSpec->fixAlpha == false){//normal gamma
						//rateOptImprove = treeStruct->OptimizeAlpha(branchPrec);
						//do NOT let alpha go too low here - on bad or random starting trees the branch lengths get crazy long
						//rateOptImprove = treeStruct->OptimizeBoundedParameter(branchPrec, mod->Alpha(), 0, 1.0e-8, 999.9, &Model::SetAlpha);
						//alphaOptImprove += treeStruct->OptimizeBoundedParameter(branchPrec, mod->Alpha(), 0, 0.05, 999.9, modnum, &Model::SetAlpha);
						alphaOptImprove += treeStruct->OptimizeBoundedParameter(modnum, branchPrec, mod->Alpha(), 0, 0.05, 999.9, &Model::SetAlpha);
						}
					}
				if(modSpec->includeInvariantSites && !modSpec->fixInvariantSites)
					pinvOptImprove += treeStruct->OptimizeBoundedParameter(modnum, branchPrec, mod->PropInvar(), 0, 1.0e-8, mod->maxPropInvar, &Model::SetPinv);
				if(modSpec->IsOrientedGap()){
					insDelImprove += treeStruct->OptimizeInsertDeleteRates(branchPrec, modnum);
					}
				if(modSpec->IsCodon() == false && modSpec->fixStateFreqs == false && modSpec->IsEqualStateFrequencies() == false && modSpec->IsEmpiricalStateFrequencies() == false)
					freqOptImprove += treeStruct->OptimizeEquilibriumFreqs(branchPrec, modnum);
				if(modSpec->fixRelativeRates == false && (modSpec->Nst() > 1 || modSpec->IsEstimateAAMatrix() || modSpec->IsTwoSerineRateMatrix()))
					rateOptImprove += treeStruct->OptimizeRelativeNucRates(branchPrec, modnum);
				}
			if(optSubsetRates){
				subsetRateImprove += treeStruct->OptimizeSubsetRates(branchPrec);
				}
			}
		improve=scaleOptImprove + trueImprove + alphaOptImprove + pinvOptImprove + flexOptImprove + omegaOptImprove + subsetRateImprove + insDelImprove;
		outman.precision(8);
		outman.UserMessageNoCR("pass%2d:+%9.3f (branch=%7.2f scale=%6.2f", i, improve, trueImprove, scaleOptImprove);
		if(optOmega) outman.UserMessageNoCR(" omega=%6.2f", omegaOptImprove);
		if(optAlpha) outman.UserMessageNoCR(" alpha=%6.2f", alphaOptImprove);

		if(optFreqs) outman.UserMessageNoCR(" freqs=%6.2f", freqOptImprove);
		if(optRelRates) outman.UserMessageNoCR(" rel rates=%6.2f", rateOptImprove);

		if(optFlex) outman.UserMessageNoCR(" flex=%6.2f", flexOptImprove);
		if(optPinv) outman.UserMessageNoCR(" pinv=%6.2f", pinvOptImprove);
		if(optInsDel){
			outman.UserMessageNoCR(" ins/del=%6.2f", insDelImprove);
			}
		if(optSubsetRates) outman.UserMessageNoCR(" subset rates=%6.2f", subsetRateImprove);
		outman.UserMessage(")");
		}

	treeStruct->nodeOptVector.clear();
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
	modPart.CopyModelPartition(&ind->modPart);
	
	dirty = ind->dirty;
	topo=ind->topo;
	}

/* 7/21/06 needs to be fixed to correspond to changes in tree for constraints
void Individual::SubtreeMutate(int subdomain, FLOAT_TYPE optPrecision, vector<int> const &subtreeMemberNodes, Adaptation *adap){
  //this version is used only by remotes when they have had a subtree defined for them
  //it will mutate only within that subtree, and because we know that the next mutation
  //will also be within that subtree we can get away without recalculating some likelihood
  //arrays

	//because we don't do model mutations during subtree mode, factor the modelMutateProb out
  	FLOAT_TYPE effectiveTopoProb=adap->topoMutateProb / (1.0/(1.0-adap->modelMutateProb));
	FLOAT_TYPE r = rnd.uniform();
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
    FLOAT_TYPE s2, s1 = params->rnd.uniform();
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
void Individual::NonSubtreeMutate(const ParallelManager *pMan, FLOAT_TYPE optPrecision, Adaptation *adap)
{//this version is used only by the master when subtree mode is active
//it will make a mutation on one of the nodes that are not contained within
//a subtree, which are in a vector that is passed in

	//because we don't do model mutations during subtree mode, factor the modelMutateProb out
  	FLOAT_TYPE effectiveTopoProb=adap->topoMutateProb / (1.0/(1.0-adap->modelMutateProb));
	FLOAT_TYPE r = rnd.uniform();

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
	  FLOAT_TYPE r = rnd.uniform();
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
    FLOAT_TYPE s2, s1 = params->rnd.uniform();
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
