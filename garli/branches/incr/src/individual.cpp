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
#include <set>
#include <stack>
using namespace std;

#include "defs.h"
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



unsigned PopRandom(std::set<unsigned> & indexSet, rng & rnd);

// returns a random element from the set and removes the element.
// 	Assumes that the set is not empty.
unsigned PopRandom(std::set<unsigned> & indexSet, rng & rnd) {
	assert (! indexSet.empty());
	const unsigned pos = rnd.random_int( indexSet.size() );
	unsigned i = 0;
	std::set<unsigned>::iterator sIt = indexSet.begin();
	for (; i != pos ; ++i, ++sIt)
		{
		assert(sIt != indexSet.end());
		}
	assert(sIt != indexSet.end());
	const unsigned toReturn = *sIt;
	indexSet.erase(sIt);
	return toReturn;
	}


#define MUTUALLY_EXCLUSIVE_MUTS

#undef VARIABLE_OPTIMIZATION

//
//
// Methods for class Individual
//
//
Individual::Individual()
  :fitness(0.0),
  dirty(1), 
  mutation_type(0), 
  mutated_brlen(0),
  accurateSubtrees(0),
  mod(0L),
  treeStruct(0L),
  reproduced(false),
  willreproduce(false), 
  willrecombine(false),
  recombinewith(-1), 
  parent(-1),
  topologyInt(-1) { 
 	treeStruct=NULL;
	mod=new Model();
	}

Individual::Individual(const Individual *other)
  :fitness(0.0),
  dirty(1), 
  mutation_type(0), 
  mutated_brlen(0),
  accurateSubtrees(0),
  mod(0L),
  treeStruct(0L),
  reproduced(false),
  willreproduce(false), 
  willrecombine(false),
  recombinewith(-1), 
  parent(-1),
  topologyInt(-1) {
 	treeStruct=NULL;
	mod=new Model();
	if (other) {
		*this = *other;
	}
}

Individual & Individual::operator=(const Individual &other) {
	if (this->treeStruct == 0L)
		this->treeStruct = new Tree();
	CopyNonTreeFields(&other);
	treeStruct->MimicTopo(other.treeStruct);
	this->dirty = other.dirty;
	treeStruct->lnL = other.fitness;
	treeStruct->SetModel(this->mod);
	return *this;
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
		treeStruct->allNodes[i]->SetAttached(false);
		
	//DZ 10-28 changing this	
	treeStruct->MimicTopo(sourceOfInformation->treeStruct);
	treeStruct->CopyClaIndeces(sourceOfInformation->treeStruct,CLAassigned);
	dirty=false;
	treeStruct->lnL=sourceOfInformation->fitness;
	treeStruct->SetModel(mod);
	}


void Individual::DoBrLenMutation() {
	mutated_brlen = treeStruct->BrlenMutate();
	if(mutated_brlen > 0) {
		mutation_type |= brlen;
		dirty = true;
		}
}


bool  Individual::DoTopoMutation(FLOAT_TYPE optPrecision, const Adaptation & adap) {
	  FLOAT_TYPE r = rnd.uniform();
	  if(r<adap.limSPRprob){
	    int reconDist = treeStruct->TopologyMutator(optPrecision, adap.limSPRrange, 0);
		if(reconDist == 1 || reconDist == -1)
			mutation_type |= randNNI;
	    else if (reconDist < 0)
	    	mutation_type |= limSPRCon;
		else
			mutation_type |= limSPR;
	  }
	  else if (r< adap.randSPRprob + adap.limSPRprob){
	    int reconDist = treeStruct->TopologyMutator(optPrecision, -1, 0);
		if (reconDist < 0) {
			if(reconDist == -1)
				mutation_type |= randNNI;
			else if(reconDist < -1 * (int)adap.limSPRrange)
				mutation_type |= randSPRCon;
			else
				mutation_type |= limSPRCon;
			}
		else {
			if (reconDist == 1)
				mutation_type |= randNNI;
			else if(reconDist >  (int) adap.limSPRrange)
				mutation_type |= randSPR;
			else
				mutation_type |= limSPR;
			}
	  } 
	  else {
		treeStruct->TopologyMutator(optPrecision, 1, 0);
		mutation_type |= randNNI;
	  }
	if (!FloatingPointEquals(treeStruct->lnL, -ONE_POINT_ZERO, 1.0e-8)) {
		fitness = treeStruct->lnL;
		dirty=false;
		}
	else 
		dirty=true;
}

void Individual::DoModelMutation()	{
	mutation_type |= mod->PerformModelMutation();
	treeStruct->MakeAllNodesDirty();
	dirty = true;
}

void Individual::Mutate(FLOAT_TYPE optPrecision, const Adaptation & adap){
	//this is the original version of mutate, and will be called by both 
	//master and remote when they are mutating a tree that does not have
	//its subtrees properly defined.

	FLOAT_TYPE r = rnd.uniform();
	//DJZ 1-5-05 Moving branch length mutation to be before topo, so that if both are performed
	//the upward sweep needed for blen optimization in the topo mutation will automatically recalc
	//nodes that were dirtied by the blen mutation, and the score of the tree can be finalized at
	//an internal node after the last branch is optimized, rather than waiting until CalcAverageFitness
	//when it will require a sweep down to the root		
#	ifndef MUTUALLY_EXCLUSIVE_MUTS
		const bool forceBrLenMut = adap.branchOptPrecision != adap.minOptPrecision;
#	else
		const bool forceBrLenMut = false;
#	endif
	if (forceBrLenMut || r >= adap.modelMutateProb + adap.topoMutateProb)
		DoBrLenMutation();

	if(r <= adap.topoMutateProb)
		DoTopoMutation(optPrecision, adap); 
	else if( r < adap.modelMutateProb + adap.topoMutateProb) 
		DoModelMutation();

	//be sure that we have an accurate score before any CLAs get invalidated
	CalcFitness(0);
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



void Individual::MakeRandomTree(unsigned nTax){
	treeStruct=new Tree();

	unsigned n = (unsigned) nTax;
	std::set<unsigned> taxset;
	for( unsigned i = 1; i <= n; i++ )
		taxset.insert(i);
		
	int placeInAllNodes = n + 1;
	
	if(!Tree::IsUsingConstraints()){
		// add nodes randomly
		for( unsigned i = 0; i < n; i++ ) {
			unsigned k = PopRandom(taxset, rnd);
			treeStruct->AddRandomNode(k, placeInAllNodes  );
			}
		}
	else{
		// add nodes randomly, ensuring that the resulting partial tree is compatible with constraints
		Bipartition mask;
		for( unsigned i = 0; i < n; i++ ) {
			unsigned k = PopRandom(taxset, rnd);
			treeStruct->AddRandomNodeWithConstraints(k, placeInAllNodes, &mask );
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
	treeStruct->AssignCLAsFromMaster();
	}


void Individual::FinishIncompleteTreeByStepwiseAddition(unsigned nTax, 
														unsigned attachesPerTaxon, 
														FLOAT_TYPE optPrecision , 
														Individual & scratchI) {
	if (treeStruct == 0L) {
		this->treeStruct = new Tree();
		this->treeStruct->AssignCLAsFromMaster();
	}
	
	std::string modelString;
	this->mod->FillGarliFormattedModelString(modelString);
	scratchI.mod->ReadGarliFormattedModelString(modelString);

	Tree *scratchT = scratchI.treeStruct;
	scratchT->SetModel(scratchI.mod);
	this->treeStruct->SetModel(this->mod);
	assert(scratchT != 0L);
	this->SetDirty();
	const unsigned n = nTax;
	std::set<unsigned> taxaIndicesToAdd;
	for( unsigned i = 1; i <= n; i++ )
		taxaIndicesToAdd.insert(i);
		
	int placeInAllNodes = n + 1;
	Bipartition mask; //mask is used for constrained trees
	Bipartition temp;
	
	// Now we walk through the tree that has been supplied (scratchT) and remove
	//	all of the added leaves from our taxaIndicesToAdd.
	
	const TreeNode * nd = scratchT->GetRootConst();
	assert(nd);
	assert(nd->left);
	nd = nd->left;
	const TreeNode * tmpNd;
	assert(nd != 0L);
	std::stack<const TreeNode *> ndStack;
	for (; nd != 0L;) {
		tmpNd = nd->left;
		if (tmpNd)
			{
			++placeInAllNodes;
			if (nd->next)
				ndStack.push(nd->next);
			nd = tmpNd;
			}
		else {
			int nn = (unsigned) nd->nodeNum;
			assert (nn >= 0 && nn <= (int)nTax);
			taxaIndicesToAdd.erase((unsigned) nn);
			if (Tree::IsUsingConstraints()) {
				mask += temp.TerminalBipart(nn);
				}
			if (nd->next) {
				nd = nd->next;
				}
			else {
				if (ndStack.empty())
					nd = 0L;
				else {
					nd = ndStack.top();
					ndStack.pop();
					}
				}
			
			}
		}
	assert(taxaIndicesToAdd.size() <= nTax);
	const int nAdded = (int)nTax - (int)taxaIndicesToAdd.size();
	outman.UserMessage("number of taxa present in starting tree: %d", nAdded);
	if (nAdded < 3)
		{
		scratchI.treeStruct->RemoveTreeFromAllClas();
		delete scratchI.treeStruct;
		scratchI.treeStruct=NULL;
		outman.UserMessage("Constructing the entire tree using stepwise addition.", nAdded);
		MakeStepwiseTree(nTax, attachesPerTaxon, optPrecision);
		}
	else
		{
		outman.UserMessage("number of taxa added:");

		// MTH lets get the branch lengths approximately correct
		//scratchI.treeStruct->Score();
		//scratchI.treeStruct->OptimizeAllBranches(optPrecision);
				
		this->ContinueBuildingStepwiseTree(nTax, attachesPerTaxon, optPrecision, scratchI, taxaIndicesToAdd, placeInAllNodes, mask);
		}
}



void Individual::MakeStepwiseTree(unsigned nTax, unsigned attachesPerTaxon, FLOAT_TYPE optPrecision ){
	treeStruct=new Tree();
	treeStruct->AssignCLAsFromMaster();

	Individual scratchI;
	scratchI.treeStruct=new Tree();
	Tree *scratchT = scratchI.treeStruct;
	scratchT->AssignCLAsFromMaster();
	scratchI.CopySecByRearrangingNodesOfFirst(scratchT, this, true);

	const unsigned n = nTax;
	std::set<unsigned> taxaIndicesToAdd;
	for( unsigned i = 1; i <= n; i++ )
		taxaIndicesToAdd.insert(i);
		
	int placeInAllNodes = n + 1;
//	ofstream stepout("stepwise.log");
	outman.UserMessage("number of taxa added:");

	Bipartition mask;//mask is used for constrained trees
	for(unsigned i = 0; i < 3; i++){//add the first 3
		unsigned k = PopRandom(taxaIndicesToAdd, rnd);
		if(!Tree::IsUsingConstraints())
			scratchT->AddRandomNode(k, placeInAllNodes  );
		else
			scratchT->AddRandomNodeWithConstraints(k, placeInAllNodes, &mask );
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
	int pos = rnd.random_int( taxaIndicesToAdd.Size() );
	int first = (taxaIndicesToAdd[pos]);
	scratchT->AddRandomNode(first, placeInAllNodes  );
	taxaIndicesToAdd -= first;
	
	//add the furthest taxon to that
	int sec = 1;
	FLOAT_TYPE maxDist = pdist[first-1][sec-1];
	for(int i=sec+1;i<=dat->NTax();i++){
		if(pdist[first-1][i-1] > maxDist){
			sec = i; 
			maxDist = pdist[first-1][sec-1];
			}
		}
	scratchT->AddRandomNode(sec, placeInAllNodes  );
	taxaIndicesToAdd -= sec;
	//add the furthest taxon to that (which may in fact be close to first, but should not have a pdist = 0 to it)
	int third = (first == 1 ? 2 : 1);
	maxDist = pdist[sec-1][third-1];
	for(int i=third+1;i<=dat->NTax();i++){
		if(pdist[sec-1][i] > maxDist && i != first && pdist[first-1][third-1] > ZERO_POINT_ZERO){
			third = i; 
			maxDist = pdist[sec-1][third-1];
			}
		}
	scratchT->AddRandomNode(third, placeInAllNodes  );
	taxaIndicesToAdd -= third;
*/

	this->ContinueBuildingStepwiseTree(nTax, attachesPerTaxon, optPrecision, scratchI, taxaIndicesToAdd, placeInAllNodes, mask);
	}

void Individual::ContinueBuildingStepwiseTree(unsigned nTax, 
											  unsigned attachesPerTaxon,
											  FLOAT_TYPE optPrecision,
											  Individual & scratchI,
											  std::set<unsigned> & taxaIndicesToAdd,
											  int & placeInAllNodes, 
											  Bipartition & mask) {
	assert(scratchI.treeStruct);
	Tree *scratchT = scratchI.treeStruct;
	assert (nTax >= taxaIndicesToAdd.size());

	CopySecByRearrangingNodesOfFirst(treeStruct, &scratchI, true);
	

	
	unsigned i = nTax - taxaIndicesToAdd.size();
	while (!taxaIndicesToAdd.empty()) {
		//select a random node
		unsigned k = PopRandom(taxaIndicesToAdd, rnd);

		//add the node randomly - this is a little odd, but for the existing swap collecting machinery
		//		to work right, the taxon to be added needs to already be in the tree
		if(!Tree::IsUsingConstraints())
			scratchT->AddRandomNode(k, placeInAllNodes );
		else
			scratchT->AddRandomNodeWithConstraints(k, placeInAllNodes, &mask );
		TreeNode *added = scratchT->allNodes[k];

		scratchT->SweepDirtynessOverTree(added);
		scratchT->OptimizeBranchesWithinRadius(added->anc, optPrecision, 0, NULL);

		//backup what we have now
		CopySecByRearrangingNodesOfFirst(treeStruct, &scratchI, true);

		FLOAT_TYPE bestBeforeOpt = scratchT->lnL;
		if (gOptimizePlausibleDuringStepwise)
			scratchT->OptimizeAllBranches(optPrecision);
		FLOAT_TYPE bestScore = scratchT->lnL;
		FLOAT_TYPE biggestOptDiff = abs(bestScore - bestBeforeOpt);
		
		//collect reconnection points - this will automatically filter for constraints
		scratchT->GatherValidReconnectionNodes(scratchT->data->NTax()*2, added, NULL, &mask);
		
//			stepout << i << "\t" << k << "\t" << bestScore << "\t";

		//start swappin
		unsigned num=0;
		//for(list<ReconNode>::iterator b = scratchT->sprRang.begin();b != scratchT->sprRang.end();b++){
		ReconList attempted;
		while(num < attachesPerTaxon && scratchT->sprRang.size() > 0){
			int connectNum = rnd.random_int(scratchT->sprRang.size());
			listIt broken = scratchT->sprRang.NthElement(connectNum);
			//try a reattachment point
			scratchT->SPRMutate(added->nodeNum, &(*broken), optPrecision, 0);
			
			
			FLOAT_TYPE currBeforeOpt = scratchT->lnL;
			FLOAT_TYPE beforeOptDiffFromBest = bestBeforeOpt - beforeOptDiffFromBest;
			
			if (gOptimizePlausibleDuringStepwise 
				&& (beforeOptDiffFromBest < 100.0 || beforeOptDiffFromBest < 2*biggestOptDiff) ) {
				scratchT->OptimizeAllBranches(optPrecision);
				FLOAT_TYPE scoreDiff = abs(currBeforeOpt - scratchT->lnL);
				if (scoreDiff > biggestOptDiff)
					biggestOptDiff = scoreDiff;
				if (bestBeforeOpt < currBeforeOpt)
					currBeforeOpt = bestBeforeOpt;
				}

			Population::RecordStepwiseAdditionTree(*scratchT);

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
		if(best != NULL) 
			scratchT->SPRMutate(added->nodeNum, best, optPrecision, 0);
		else
			scratchT->Score();
		scratchI.CalcFitness(0);

//		stepout << scratchT->lnL << endl;
		CopySecByRearrangingNodesOfFirst(treeStruct, &scratchI, true);

		//outman.UserMessage(" %d %f", i+1, scratchT->lnL);
		outman.UserMessageNoCR(" %d ", i+1);
		outman.flush();
		//when we've added half the taxa optimize alpha, flex or omega 
		if (i == (nTax/2) && (modSpec.IsCodon() || mod->NRateCats() > 1) && modSpec.fixAlpha == false){
			FLOAT_TYPE rateOptImprove = 0.0;
			if(modSpec.IsCodon())//optimize omega even if there is only 1
				rateOptImprove = scratchT->OptimizeOmegaParameters(optPrecision);
			else if(mod->NRateCats() > 1){
				if(modSpec.IsFlexRateHet()){//Flex rates
					rateOptImprove = ZERO_POINT_ZERO;
					//for the first pass, use gamma to get in the right ballpark
					//rateOptImprove = scratchT->OptimizeAlpha(optPrecision);
					rateOptImprove = treeStruct->OptimizeBoundedParameter(optPrecision, mod->Alpha(), 0, 0.1, 999.9, &Model::SetAlpha);
					rateOptImprove += scratchT->OptimizeFlexRates(optPrecision);
					}
				else if(modSpec.fixAlpha == false){//normal gamma
					//rateOptImprove = scratchT->OptimizeAlpha(optPrecision);
					rateOptImprove = treeStruct->OptimizeBoundedParameter(optPrecision, mod->Alpha(), 0, 0.05, 999.9, &Model::SetAlpha);
					}
				}
			outman.UserMessage("\nOptimizing parameters... improved %f lnL", rateOptImprove);
		//	this used to depend on param improvement - not sure why
		//	if(rateOptImprove > 0.0){
				scratchT->Score();
				FLOAT_TYPE start=scratchT->lnL;
				scratchT->OptimizeAllBranches(optPrecision);
				FLOAT_TYPE bimprove = scratchT->lnL - start;
				outman.UserMessage("\nOptimizing branchlengths... improved %f lnL", bimprove);
	//			}
			}
		++i;
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
					if(c=='r'||c=='R'||c=='b'||c=='B'||c=='e'||c=='E'||c=='a'||c=='A'||c=='p'||c=='P'||c=='i'||c=='I'||c=='f'||c=='o'||c=='O') foundModel=true;
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
		string modString;
		do{
			c=stf.get();
			modString += c;
			//}while(c != '(' && c != '\r' && c != '\n' && !stf.eof());
			}while(stf.peek() != '(' && stf.peek() != '\r' && stf.peek() != '\n' && !stf.eof());
		while((stf.peek() == '\n' || stf.peek() == '\r') && stf.eof() == false) stf.get(c);
		mod->ReadGarliFormattedModelString(modString);
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

		//now allowing polytomies, since they will be taken care of in Population::SeedPopulationWithStartingTree
		treeStruct=new Tree(treeString.c_str(), numericalTaxa, true);

		std::string violCon;
		if (!treeStruct->ObeysConstraints(&violCon)) {
			std::string msg = "Starting tree is not compatible with the constraint:\n";
			msg.append(violCon);
			msg.append(1, '\n');
			throw ErrorException(msg.c_str());
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
/*			if(treeStruct->constraints.size() == 0)
				outman.UserMessage("No starting tree found in file %s, creating random tree", fname);
			else 
				outman.UserMessage("No starting tree found in file %s, creating random tree (compatible with constraint)", fname);
*/			}

		if(foundModel==true){
			outman.UserMessage("Obtained starting or fixed model parameter values from file %s", fname);
/*			string m;
			mod->FillGarliFormattedModelString(m);
			outman.UserMessage("%s", m.c_str());
*/			}
		else{
			//this checks whether we have already gotten some parameter values from file, which might have come from a garli block in the datafile
			if(!(modSpec.GotAnyParametersFromFile())){
				outman.UserMessage("No starting model parameter values found in %s\nUsing default parameter values", fname);
/*				string m;
				mod->FillGarliFormattedModelString(m);
				outman.UserMessage("%s", m.c_str());
*/				}
			}
			
		outman.UserMessage("");
		}

	mod->UpdateQMat();
	stf.close();
	delete []temp;
	}

void Individual::GetStartingTreeFromNCL(const NxsTreesBlock *treesblock, 
										int rank, 
										int , //nTax, 
										bool , //restart, // = false
										bool demandAllTaxa /* = true */) {
	assert(treeStruct == NULL);

	int totalTrees = treesblock->GetNumTrees();

	int effectiveRank = rank % totalTrees;
	
	//we will get the tree string from NCL with taxon numbers (starting at 1), regardless of how it was initially read in 
	const NxsFullTreeDescription &t = treesblock->GetFullTreeDescription(effectiveRank);
	this->ReadNxsFullTreeDescription(t, demandAllTaxa);
	}

void Individual::ReadNxsFullTreeDescription(const NxsFullTreeDescription &t, bool demandAllTaxa) {
	if(demandAllTaxa && t.AllTaxaAreIncluded() == false)
		throw ErrorException("Starting tree description must contain all taxa.");
	string ts = t.GetNewick();
	ts += ";";
	treeStruct=new Tree(ts.c_str(), true, true, !demandAllTaxa);

	std::string violCon;
	if (!treeStruct->ObeysConstraints(&violCon)) {
		std::string msg = "Starting tree is not compatible with the constraint:\n";
		msg.append(violCon);
		msg.append(1, '\n');
		throw ErrorException(msg.c_str());
	}

	treeStruct->AssignCLAsFromMaster();

	mod->UpdateQMat();
	}


void Individual::RefineStartingConditions(bool optModel, FLOAT_TYPE branchPrec){
	if(optModel && mod->NRateCats() > 1 && modSpec.IsNonsynonymousRateHet() == false && modSpec.gotFlexFromFile == false) outman.UserMessage("optimizing starting branch lengths and rate heterogeneity parameters...");
	else if(modSpec.IsCodon()) outman.UserMessage("optimizing starting branch lengths and dN/dS (aka omega) parameters...");
	else outman.UserMessage("optimizing starting branch lengths...");
	FLOAT_TYPE improve=(FLOAT_TYPE)999.9;
	CalcFitness(0);

	for(int i=1;improve > branchPrec;i++){
		FLOAT_TYPE rateOptImprove=0.0, pinvOptImprove = 0.0, optImprove=0.0, scaleImprove=0.0;
		
		CalcFitness(0);
		FLOAT_TYPE passStart=Fitness();
		
		optImprove=treeStruct->OptimizeAllBranches(branchPrec);
		CalcFitness(0);

		FLOAT_TYPE trueImprove= Fitness() - passStart;
		assert(trueImprove >= -1.0);
		if(trueImprove < ZERO_POINT_ZERO) trueImprove = ZERO_POINT_ZERO;

		vector<FLOAT_TYPE> blens;
		treeStruct->StoreBranchlengths(blens);
		scaleImprove=treeStruct->OptimizeTreeScale(branchPrec);
		CalcFitness(0);
		//if some of the branch lengths were at the minimum or maximum boundaries the scale optimization
		//can actually worsen the score.  If so, return them to their original lengths.
		if(scaleImprove < ZERO_POINT_ZERO){
			treeStruct->RestoreBranchlengths(blens);
			CalcFitness(0);
			scaleImprove = ZERO_POINT_ZERO;
			}

		CalcFitness(0);
		if(optModel){
			if(modSpec.IsCodon())//optimize omega even if there is only 1
				rateOptImprove = treeStruct->OptimizeOmegaParameters(branchPrec);
			else if(mod->NRateCats() > 1){
				if(modSpec.IsFlexRateHet()){//Flex rates
					rateOptImprove = ZERO_POINT_ZERO;
					//no longer doing alpha first, it was too hard to know if the flex rates had been partially optimized
					//already during making of a stepwise tree
					//if(i == 1) rateOptImprove = treeStruct->OptimizeAlpha(branchPrec);
					//if(i == 1 && modSpec.gotFlexFromFile==false) rateOptImprove = treeStruct->OptimizeBoundedParameter(branchPrec, mod->Alpha(), 0, 1.0e-8, 999.9, &Model::SetAlpha);
					rateOptImprove += treeStruct->OptimizeFlexRates(branchPrec);
					}
				else if(modSpec.fixAlpha == false){//normal gamma
					//rateOptImprove = treeStruct->OptimizeAlpha(branchPrec);
					//do NOT let alpha go too low here - on bad or random starting trees the branch lengths get crazy long
					//rateOptImprove = treeStruct->OptimizeBoundedParameter(branchPrec, mod->Alpha(), 0, 1.0e-8, 999.9, &Model::SetAlpha);
					rateOptImprove = treeStruct->OptimizeBoundedParameter(branchPrec, mod->Alpha(), 0, 0.05, 999.9, &Model::SetAlpha);
					}
				}
			if(modSpec.includeInvariantSites && !modSpec.fixInvariantSites)
				pinvOptImprove = treeStruct->OptimizeBoundedParameter(branchPrec, mod->PropInvar(), 0, 1.0e-8, mod->maxPropInvar, &Model::SetPinv);
			}
		improve=scaleImprove + trueImprove + rateOptImprove + pinvOptImprove;
		outman.precision(8);
		outman.UserMessageNoCR("pass %-2d: +%10.4f (branch=%7.2f scale=%6.2f", i, improve, trueImprove, scaleImprove);
		if(optModel && modSpec.IsCodon()) outman.UserMessageNoCR(" omega=%6.2f", rateOptImprove);
		else if(optModel && mod->NRateCats() > 1){
			if(modSpec.IsFlexRateHet() == false && modSpec.fixAlpha == false) outman.UserMessageNoCR(" alpha=%6.2f", rateOptImprove);
			else if(modSpec.IsFlexRateHet() && modSpec.gotFlexFromFile == false) outman.UserMessageNoCR(" flex rates=%6.2f", rateOptImprove);
			}
		if(optModel && modSpec.includeInvariantSites && !modSpec.fixInvariantSites) outman.UserMessageNoCR(" pinv=%6.2f", pinvOptImprove);
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
	mod->CopyModel(ind->mod);
	
	dirty = ind->dirty;
	SetTopo(ind->GetTopo());
	}
