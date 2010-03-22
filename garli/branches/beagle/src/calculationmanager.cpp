// GARLI version 1.00 source code
// Copyright 2005-2010 Derrick J. Zwickl
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


#include "defs.h"
#include "treenode.h"
#include "utility.h"
#include "calculationmanager.h"

//extern ClaManager claMan;

ClaManager *CalculationManager::claMan = NULL;
PmatManager *CalculationManager::pmatMan = NULL;
const SequenceData *CalculationManager::data = NULL;

ClaManager *NodeClaManager::claMan = NULL;
PmatManager *NodeClaManager::pmatMan = NULL;


const char *AdvanceDataPointer(const char *arr, int num);

#ifdef USE_BEAGLE

void InterpretBeagleResourceFlags(long flags, string &list){
    if (flags & BEAGLE_FLAG_PRECISION_DOUBLE)	list += " DOUBLE";
    if (flags & BEAGLE_FLAG_PRECISION_SINGLE)	list += " SINGLE";
    
	if (flags & BEAGLE_FLAG_COMPUTATION_ASYNCH) list += " ASYNCH";
    if (flags & BEAGLE_FLAG_COMPUTATION_SYNCH)  list += " SYNCH";
    
	if (flags & BEAGLE_FLAG_EIGEN_COMPLEX)		list += " COMPLEX";
	if (flags & BEAGLE_FLAG_EIGEN_REAL)			list += " REAL";

	if (flags & BEAGLE_FLAG_SCALING_MANUAL)		list += " MANUAL_SCALING";
	if (flags & BEAGLE_FLAG_SCALING_AUTO)		list += " AUTO_SCALING";
	if (flags & BEAGLE_FLAG_SCALING_ALWAYS)		list += " ALWAYS_SCALING";

	if (flags & BEAGLE_FLAG_SCALERS_RAW)		list += " RAW_SCALERS";
    if (flags & BEAGLE_FLAG_SCALERS_LOG)		list += " LOG_SCALERS";

    if (flags & BEAGLE_FLAG_VECTOR_SSE)			list += " SSE";
	if (flags & BEAGLE_FLAG_THREADING_OPENMP)	list += " OPENMP";
	
	if (flags & BEAGLE_FLAG_PROCESSOR_CPU)		list += " CPU";
    if (flags & BEAGLE_FLAG_PROCESSOR_GPU)		list += " GPU";
    if (flags & BEAGLE_FLAG_PROCESSOR_FPGA)		list += " FPGA";
    if (flags & BEAGLE_FLAG_PROCESSOR_CELL)		list += " CELL";	
	}

// print possible beagle resources
void OutputBeagleResources(){
    BeagleResourceList* rList;
    rList = beagleGetResourceList();
    outman.UserMessageNoCR("Available resources:\n");
    for (int i = 0; i < rList->length; i++) {
        outman.UserMessageNoCR("\tResource %i:\n\t\tName : %s\n", i, rList->list[i].name);
        outman.UserMessageNoCR("\t\tDesc : %s\n", rList->list[i].description);
        outman.UserMessageNoCR("\t\tFlags:");
		string flagList;
		InterpretBeagleResourceFlags(rList->list[i].supportFlags, flagList);
		outman.UserMessage("%s", flagList.c_str());
		}
    outman.UserMessageNoCR("\n");
	}

//print actual beagle resources used by specific instance
void CalculationManager::OutputInstanceDetails(const BeagleInstanceDetails *det) const{
    outman.UserMessageNoCR("Instance details:\n");
	outman.UserMessageNoCR("\t\tnumber: %d\n", det->resourceNumber);
	outman.UserMessageNoCR("\t\tresource name: %s\n", det->resourceName);
	outman.UserMessageNoCR("\t\timplementation name: %s\n", det->implName);
	outman.UserMessageNoCR("\t\tFlags: ");
	string flagList;
	InterpretBeagleResourceFlags(det->flags, flagList);
	outman.UserMessage("%s", flagList.c_str());
	}

void CalculationManager::CheckBeagleReturnValue(int err, const char *funcName) const{

	string mess;
	if(err == BEAGLE_SUCCESS || err >= 0)
		return;

	else if(err == BEAGLE_ERROR_GENERAL)
		mess = "General Beagle error in ";

	else if(err == BEAGLE_ERROR_OUT_OF_MEMORY)
		mess = "Beagle out-of-memory error in ";

	else if(err == BEAGLE_ERROR_UNIDENTIFIED_EXCEPTION)
		mess = "Beagle unidentified exception error in ";

	else if(err == BEAGLE_ERROR_UNINITIALIZED_INSTANCE)
		mess = "Beagle uninitialized instance error in ";

	else if(err == BEAGLE_ERROR_OUT_OF_RANGE)
		mess = "Beagle out-of-range error in ";

	else if(err == BEAGLE_ERROR_NO_RESOURCE)
		mess = "Beagle no-resource error in ";

	else if(err == BEAGLE_ERROR_NO_IMPLEMENTATION)
		mess = "Beagle no-implementation error in ";

	else if(err == BEAGLE_ERROR_FLOATING_POINT)
		mess = "Beagle floating point error in ";

	else
		mess = "Unknown Beagle error in ";

	mess += funcName;
	if(termOnBeagleError)
		throw ErrorException("%s", mess.c_str());

	else
		outman.UserMessage("%s", mess.c_str());

	return;
	}

void CalculationManager::InitializeBeagle(int nTips, int nClas, int nHolders, int nstates, int nchar, int nrates){
	assert(useBeagle);

	termOnBeagleError = false;
	long req_flag = 0;
	//the fact that GPU implies SP is taken care of elsewhere
	long pref_flag = (gpuBeagle ? BEAGLE_FLAG_PROCESSOR_GPU : BEAGLE_FLAG_PROCESSOR_CPU);
	if(singlePrecBeagle)
		req_flag |= BEAGLE_FLAG_PRECISION_SINGLE;
	//must accumulate log rescalers for my implementation
	if(rescaleBeagle)
		req_flag |= BEAGLE_FLAG_SCALERS_LOG;

	outman.UserMessage("BEAGLE INITIALIZING ...");
	OutputBeagleResources();

	outman.UserMessage("CREATING BEAGLE INSTANCE ...");
//without ambiguity
/*	
	int tipCount = nTips;
	int partialsCount = nClas;
	//I think that all compacts are tips, but not all tips must be compact. They always are in the nonambig case though.
	int compactCount = tipCount;
*/
	//to allow ambiguity we need to create partials for the tips with ambiguity, and normal tip states for the others
	int tipCount = nTips;
	int ambigTips = data->NumTaxaWithPartialAmbig();
	int normalTips = nTips - ambigTips;
	int partialsCount = nClas + ambigTips;
	int compactCount = normalTips;

	int eigCount = nClas;
	int matrixCount = (nClas * 3) * 2;//x3 for the pmats, d1mats and d2mats, x2 for both internals and terms
	//DEBUG - trying scaling
	//try one scaler per cla, as in normal garli.  These are doubles rather than ints though, so larger.
	//scaler for a given cla will share same index, and will be cumulative, as mine are now
	//add one for a destinationScaleWrite which will always be the last and will be used in all calls for scratch
	int scalerCount = (rescaleBeagle ? nClas + 1 : 0);
	int resourceList[1] = {NULL};
	int resourceListCount = 0;

	BeagleInstanceDetails det;

	//this returns either the instance number or a beagle error, which is somewhat annoying
   	beagleInst = beagleCreateInstance(tipCount, 
		partialsCount, 
		compactCount, 
		nstates, 
		nchar, 
		eigCount, 
		matrixCount, 
		nrates, 
		scalerCount, 
		resourceList, 
		resourceListCount, 
		pref_flag, 
		req_flag,
		&det);

	CheckBeagleReturnValue(
		beagleInst, 
		"beagleCreateInstance");

	outman.DebugMessage("INITING");

	outman.DebugMessage("BEAGLE ALLOCATIONS:");
	outman.DebugMessage("\tstates: %d char: %d rates: %d", nstates, nchar, nrates);
	outman.DebugMessage("\ttips: %d", tipCount);
	outman.DebugMessage("\tcompact tips arrays (no partial ambiguity): %d", compactCount);
	outman.DebugMessage("\tpartial tips arrays (some partial ambiguity): %d", ambigTips);
	outman.DebugMessage("\ttotal partial arrays: %d", partialsCount);
	outman.DebugMessage("\teigen solutions: %d", eigCount);
	outman.DebugMessage("\ttransition matrices (pmats and derivs): %d", matrixCount);
	outman.DebugMessage("\trescaling arrays %d", scalerCount);

	OutputInstanceDetails(&det);

	SendTipDataToBeagle();

	vector<double> counts;
	for(int pat = 0;pat < data->NChar();pat++)
		counts.push_back((double) data->Count(pat));

	beagleSetPatternWeights(beagleInst, &(counts[0]));

	outman.UserMessage("#######################################################");
	}

void CalculationManager::SendTipDataToBeagle(){
	outman.DebugMessage("SENDING DATA");
    char convert[16]={-1, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4};
	
	int nstates = data->NStates();
	for(int t = 0;t < data->NTax();t++){
		bool partialAmbig = data->TaxonHasPartialAmbig(t);
		outman.DebugMessageNoCR("tax %d ", t);

		const unsigned char *dataString = data->GetRow(t);
		if(partialAmbig){
			//ambiguity if currently for nuc only (all versions, not just beagle)
			assert(nstates == 4);
			outman.DebugMessage("(some ambiguity)");
			vector<double> tipPartial;
			for(int c = 0;c < data->NChar();c++){
				for(int s = 0;s < nstates;s++){
					((dataString[c] & (1 << s)) ? tipPartial.push_back(1.0) : tipPartial.push_back(0.0));
					}
				}
			CheckBeagleReturnValue(
				beagleSetTipPartials(beagleInst, t, &(tipPartial[0])),
				"beagleSetTipStates");
			}
		else{
			outman.DebugMessage("(no ambiguity)");
			vector<int> dat;
			
			for(int c = 0;c < data->NChar();c++){
				if(nstates == 4){
					//nucleotide data needs to be converted from the bitwise format to 0, 1, 2, 3, 4
					assert(dataString[c] < 16);
					assert(convert[dataString[c]] < 5);
					dat.push_back(convert[dataString[c]]);
					}
				else{
					//for non-nuc the data is already in the correct indexing scheme, but must be converted from
					//chars to ints
					dat.push_back(dataString[c]);
					}
				}
			CheckBeagleReturnValue(
				beagleSetTipStates(beagleInst, t, &(dat[0])),
				"beagleSetTipStates");
			}
		}
	outman.DebugMessage("DATA SENT");
	}
#endif

/*This is called by CalculateLikelihood and CalculateDerivatives and amasses all of the transmat, cla and scoring
operations necessary to get whatever is asked for.  When deriviatives is false, effectiveRoot will be treated as 
just that, the node at which all clas are pulled to.  It must be internal.  Which clas are combined is determined by this function.
In the derivative case it assumes that the effective root is the node a the top of the branch, i.e. the one that
corresponds to the branch for which the derivatives are requested.  It may be terminal.
Nothing is acutally calculated here - operations are just put into the queues*/
void CalculationManager::DetermineRequiredOperations(const TreeNode *effectiveRoot, bool derivatives){
	//outman.DebugMessage("**ENTERING DETERMINE REQUIRED OPS, ROOT = %d, derivatives = %d", effectiveRoot->nodeNum, derivatives);

	int dest, end1, end2, finalPmat;
	dest = end1 = end2 = finalPmat = -1;
	if(! derivatives){
		//here we are calculating the lnL only, and don't case which branch it is integrated on.  pick an internal one
		if(effectiveRoot->left->IsInternal()){
			end1 = effectiveRoot->myMan.ULHolderIndex;
			end2 = effectiveRoot->left->myMan.downHolderIndex;
			finalPmat = effectiveRoot->left->myMan.transMatIndex;
			}

		else if(effectiveRoot->right->IsInternal()){
			end1 = effectiveRoot->myMan.URHolderIndex;
			end2 = effectiveRoot->right->myMan.downHolderIndex;
			finalPmat = effectiveRoot->right->myMan.transMatIndex;
			}

		else if(effectiveRoot->anc){
			end1 = effectiveRoot->myMan.downHolderIndex;
			finalPmat = effectiveRoot->myMan.transMatIndex;
			if(effectiveRoot == effectiveRoot->anc->left)
				end2 = effectiveRoot->anc->myMan.ULHolderIndex;
			else if(effectiveRoot == effectiveRoot->anc->right)
				end2 = effectiveRoot->anc->myMan.URHolderIndex;
			else if(effectiveRoot == effectiveRoot->anc->left->next)//middle descendent of the effectiveRoot
				end2 = effectiveRoot->anc->myMan.downHolderIndex;
			else
				assert(0);
			}

		else{//if effectiveRoot is the true root of the tree
			end1 = effectiveRoot->myMan.downHolderIndex;
			end2 = effectiveRoot->left->next->myMan.downHolderIndex;
			finalPmat = effectiveRoot->left->next->myMan.transMatIndex;
			}
		}
	else{//derivatives
		//here we are specifically saying that the effectiveRoot is the "upper" of the two nodes, making it's anc
		//automatically the other node.  Figure out specifically which clas those are
		assert(! effectiveRoot->IsRoot());
		if(effectiveRoot->IsTerminal())
			end1 = -(effectiveRoot->nodeNum);
		else
			end1 = effectiveRoot->myMan.downHolderIndex;
		finalPmat = effectiveRoot->myMan.transMatIndex;
		if(effectiveRoot->anc->left == effectiveRoot)
			end2 = effectiveRoot->anc->myMan.ULHolderIndex;
		else if(effectiveRoot->anc->right == effectiveRoot)
			end2 = effectiveRoot->anc->myMan.URHolderIndex;
		else if(effectiveRoot->anc->left->next == effectiveRoot){//middle descendent of the effectiveRoot
			assert(effectiveRoot->anc->IsRoot());
			end2 = effectiveRoot->anc->myMan.downHolderIndex;
			}
		else 
			assert(0);
		}

	//this will recurse and grab all of the operations needed to get the end1 and end2 partials, and put them into the operationSetQueue
	//or, the new way they will all be put into nodeOps passed in, and then sorted out below
	list<NodeOperation> nodeOps;
	AccumulateOpsOnPath(end1, nodeOps);
	AccumulateOpsOnPath(end2, nodeOps);

	//this sorts the operations by dep level (with operator <), then graps all of a particular rank into a blocking set
	nodeOps.sort();

	operationSetQueue.clear();
	BlockingOperationsSet thisSet;

	int lvl = 1;
	list<NodeOperation>::iterator start = nodeOps.begin();
	list<NodeOperation>::iterator end = nodeOps.begin();
	while(end != nodeOps.end()){
		int num = 0;
		while(end != nodeOps.end() && end->claOp.depLevel == lvl){
			end++;
			num++;
			}
		//outman.UserMessage("level %d = %d nodes", lvl, num);
		for(list<NodeOperation>::iterator nit = start;nit != end;nit++){
			thisSet.claOps.push_back((*nit).claOp);
			thisSet.pmatOps.push_back((*nit).transOp1);
			thisSet.pmatOps.push_back((*nit).transOp2);
			}
		operationSetQueue.push_back(thisSet);
		start = end;
		lvl++;
		thisSet.claOps.clear();
		thisSet.pmatOps.clear();
		}

	//the matrix for the branch between end1 and end2 for the final operation 
	//there is no cla opt in this final operation, the calculation will instead 
	//be stored in the scoreOps list
	TransMatOperation final(finalPmat, 0, pmatMan->GetMutableHolder(finalPmat)->GetEdgelen(), derivatives);
	BlockingOperationsSet fin;
	fin.pmatOps.push_back(final);
	operationSetQueue.push_back(fin);
	scoreOps.push_back(ScoringOperation(-1, end1, end2, finalPmat, derivatives));
	}

void CalculationManager::OutputOperationsSummary() const{
	int numPmats = 0, numClas = 0, numSets = 0;
	for(list<BlockingOperationsSet>::const_iterator it = operationSetQueue.begin();it != operationSetQueue.end();it++){
		numSets++;
		for(list<TransMatOperation>::const_iterator pit = (*it).pmatOps.begin();pit != (*it).pmatOps.end();pit++)
			numPmats++;
		for(list<ClaOperation>::const_iterator cit = (*it).claOps.begin();cit != (*it).claOps.end();cit++)
			numClas++;
		}
	outman.DebugMessage("REQ OPS : %d op sets, %d pmats, %d clas", numSets, numPmats, numClas);
	}

ScoreSet CalculationManager::CalculateDerivatives(const TreeNode *effectiveRoot){
	outman.DebugMessage("#########################\nENTERING CALC DERIVATIVES, ROOT = %d", effectiveRoot->nodeNum);

	DetermineRequiredOperations(effectiveRoot, true);

	UpdateAllConditionals();

	//This assumes that there is only one score op in scoreOps
	assert(scoreOps.size() == 1);
	ScoreSet values;
	for(list<ScoringOperation>::iterator sit = scoreOps.begin() ; sit != scoreOps.end() ; sit++)
		values = PerformScoringOperation(&(*sit));

	operationSetQueue.clear();
	scoreOps.clear();
	return values;
	}	

FLOAT_TYPE CalculationManager::CalculateLikelihood(const TreeNode *effectiveRoot){
	outman.DebugMessage("#########################\nENTERING CALC LIKELIHOOD, ROOT = %d", effectiveRoot->nodeNum);

	DetermineRequiredOperations(effectiveRoot, false);

	OutputOperationsSummary();

	UpdateAllConditionals();

	//This assumes that there is only one score op in scoreOps
	assert(scoreOps.size() == 1);
	//only the lnL field of the values struct will be filled here
	ScoreSet values;
	for(list<ScoringOperation>::iterator sit = scoreOps.begin() ; sit != scoreOps.end() ; sit++)
		values = PerformScoringOperation(&(*sit));

	operationSetQueue.clear();
	scoreOps.clear();
	return values.lnL;
	}	

//this just performs all of the transmat and cla ops that were queued
//will generally be followed by a call to PerformScoring to get lnL or derivs
void CalculationManager::UpdateAllConditionals(){	
	//DEBUG
	OutputOperationsSummary();

#define BATCHED_CALLS
#ifdef BATCHED_CALLS
	//this just yoinks the pmat ops from each BlockingOpSet into a new list
	//note that PerformTransMatOperationBatch assumes that all the ops either
	//need derivs or not.  So, avoid adding ops that need derivs here.
	list<TransMatOperation> tops;
	for(list<BlockingOperationsSet>::iterator it = operationSetQueue.begin();it != operationSetQueue.end();it++){
		if(! (*it).pmatOps.begin()->calcDerivs)
			tops.splice(tops.end(), (*it).pmatOps);
		}
	if(tops.size() > 0)
		PerformTransMatOperationBatch(tops);

	//Any transmats that need derivs will still be included in a blockingOpSet, so will get done here.  Otherwise
	//no transmat calls will be done here because the ops were pulled out above to be done in one big batch
	for(list<BlockingOperationsSet>::const_iterator it = operationSetQueue.begin();it != operationSetQueue.end();it++){
		for(list<TransMatOperation>::const_iterator pit = (*it).pmatOps.begin();pit != (*it).pmatOps.end();pit++)
			PerformTransMatOperation(&(*pit));
		if((*it).claOps.size() > 0)
			PerformClaOperationBatch((*it).claOps);
		}
#else
	//do the actual operations.  right now it is doing 2 pmat ops followed by the cla op that requires then and looping
	//since there are enough pmats it would make more sense to put them into a single block to be done together
	for(list<BlockingOperationsSet>::const_iterator it = operationSetQueue.begin();it != operationSetQueue.end();it++){
		for(list<TransMatOperation>::const_iterator pit = (*it).pmatOps.begin();pit != (*it).pmatOps.end();pit++)
			PerformTransMatOperation(&(*pit));
		for(list<ClaOperation>::const_iterator cit = (*it).claOps.begin();cit != (*it).claOps.end();cit++)
			PerformClaOperation(&(*cit));
		}
#endif
	}

int CalculationManager::AccumulateOpsOnPath(int holderInd, list<NodeOperation> &nodeOps){
	//it is ok to call this with a tip, because that tip could be a terminal branch that is having derivs
	//calculated.  accumulate nothing in that case
	if(holderInd < 0)
		return 0;
	CondLikeArrayHolder *holder = claMan->GetMutableHolder(holderInd);
	if(!claMan->IsDirty(holderInd)){
		holder->depLevel = 0;
		return 0;
		}

	int depLevel1, depLevel2;
	depLevel1 = depLevel2 = 0;
	if(holder->holderDep1 >= 0){
		if(claMan->IsDirty(holder->holderDep1))
			depLevel1 = AccumulateOpsOnPath(holder->holderDep1, nodeOps);
		else
			claMan->GetMutableHolder(holder->holderDep1)->depLevel = 0;
		}
	if(holder->holderDep2 >= 0){
		if(claMan->IsDirty(holder->holderDep2))
			depLevel2 = AccumulateOpsOnPath(holder->holderDep2, nodeOps);
		else
			claMan->GetMutableHolder(holder->holderDep2)->depLevel = 0;
		}
	
	holder->depLevel = max(depLevel1, depLevel2) + 1;

	BlockingOperationsSet opSet;
	opSet.claOps.push_back(ClaOperation(holderInd, holder->holderDep1, holder->holderDep2, holder->transMatDep1, holder->transMatDep2, holder->depLevel));
	//DEBUG - second arg here is model index.  Not sure how that will be dealt with
	opSet.pmatOps.push_back(TransMatOperation(holder->transMatDep1, 0, pmatMan->GetMutableHolder(holder->transMatDep1)->GetEdgelen(), false));
	opSet.pmatOps.push_back(TransMatOperation(holder->transMatDep2, 0, pmatMan->GetMutableHolder(holder->transMatDep2)->GetEdgelen(), false));

	operationSetQueue.push_back(opSet);

	nodeOps.push_back(NodeOperation(ClaOperation(holderInd, holder->holderDep1, holder->holderDep2, holder->transMatDep1, holder->transMatDep2, holder->depLevel),
						TransMatOperation(holder->transMatDep1, 0, pmatMan->GetMutableHolder(holder->transMatDep1)->GetEdgelen(), false),
						TransMatOperation(holder->transMatDep2, 0, pmatMan->GetMutableHolder(holder->transMatDep2)->GetEdgelen(), false)));


	return holder->depLevel;
	}

void CalculationManager::PerformClaOperation(const ClaOperation *theOp){

#ifdef USE_BEAGLE
	//DEBUG
	useBeagle = true;
#endif
//	outman.DebugMessage("**ENTERING PERFORM CLA");	

	if(! useBeagle){
		int term1 = 0;
		int term2 = 0;
		if(theOp->childCLAIndex1 < 0)
			term1 = -theOp->childCLAIndex1;
	
		if(theOp->childCLAIndex2 < 0)
			term2 = -theOp->childCLAIndex2;

		if( term1 && term2 ){
			CalcFullCLATerminalTerminal(
				claMan->GetCla(theOp->destCLAIndex),
				//DEBUG - should probably figure some way not to do these static cassts
				(char*) static_cast<const NucleotideData *>(data)->GetAmbigString(term1-1),
				(char*) static_cast<const NucleotideData *>(data)->GetAmbigString(term2-1),
				pmatMan->GetPmat(theOp->transMatIndex1),
				pmatMan->GetPmat(theOp->transMatIndex2), 
				pmatMan->GetCorrespondingModel(theOp->transMatIndex1));
			}
		else if( ! (term1 || term2)){
			CalcFullClaInternalInternal(
				claMan->GetCla(theOp->destCLAIndex),
				claMan->GetCla(theOp->childCLAIndex1),
				claMan->GetCla(theOp->childCLAIndex2),
				pmatMan->GetPmat(theOp->transMatIndex1),
				pmatMan->GetPmat(theOp->transMatIndex2),
				pmatMan->GetCorrespondingModel(theOp->transMatIndex1));
			}
		else if( term1 ){
			this->CalcFullCLAInternalTerminal(
				claMan->GetCla(theOp->destCLAIndex),
				claMan->GetCla(theOp->childCLAIndex2),
				(char*) static_cast<const NucleotideData *>(data)->GetAmbigString(term1-1),
				pmatMan->GetPmat(theOp->transMatIndex2),
				pmatMan->GetPmat(theOp->transMatIndex1),
				pmatMan->GetCorrespondingModel(theOp->transMatIndex1), 
				NULL);
			}
		else{
			this->CalcFullCLAInternalTerminal(
				claMan->GetCla(theOp->destCLAIndex),
				claMan->GetCla(theOp->childCLAIndex1),
				(char*) static_cast<const NucleotideData *>(data)->GetAmbigString(term2-1),
				pmatMan->GetPmat(theOp->transMatIndex1),
				pmatMan->GetPmat(theOp->transMatIndex2),
				pmatMan->GetCorrespondingModel(theOp->transMatIndex1),
				NULL);
			}
		}
	else{
#ifdef USE_BEAGLE

		//not sure if this is right - will always use a single scale array for destWrite (essentially
		//scratch space, I think) and then pass a cumulative scaler to actually keep track of the scaling
		int destinationScaleWrite = (rescaleBeagle ? claMan->NumClas() : BEAGLE_OP_NONE);
		int	destinationScaleRead = BEAGLE_OP_NONE;

		int operationTuple[7] = {PartialIndexForBeagle(theOp->destCLAIndex),
                                destinationScaleWrite,
                                destinationScaleRead,
								PartialIndexForBeagle(theOp->childCLAIndex1),
 								PmatIndexForBeagle(theOp->transMatIndex1),
 								PartialIndexForBeagle(theOp->childCLAIndex2),
								PmatIndexForBeagle(theOp->transMatIndex2)};

		outman.DebugMessageNoCR("PARTS\t");
		outman.DebugMessageNoCR("\tD\t%d (%d)", PartialIndexForBeagle(theOp->destCLAIndex), theOp->destCLAIndex);
		outman.DebugMessageNoCR("\tC\t%d (%d)\t%d (%d)", PartialIndexForBeagle(theOp->childCLAIndex1), theOp->childCLAIndex1, PartialIndexForBeagle(theOp->childCLAIndex2), theOp->childCLAIndex2);
		outman.DebugMessage("\tP\t%d (%d)\t%d (%d)", PmatIndexForBeagle(theOp->transMatIndex1), theOp->transMatIndex1, PmatIndexForBeagle(theOp->transMatIndex2), theOp->transMatIndex2);
		
		int instanceCount = 1;
		int operationCount = 1;

		//accumulate rescaling factors - For scale arrays my indexing scheme and Beagle's happen to be the same, 
		//and my negative (tip) corresponds to a NULL in beagle
		int cumulativeScaleIndex = BEAGLE_OP_NONE;
		if(rescaleBeagle){
			cumulativeScaleIndex = theOp->destCLAIndex;
			AccumulateRescalers(cumulativeScaleIndex, theOp->childCLAIndex1, theOp->childCLAIndex2);
			}

		CheckBeagleReturnValue(
			beagleUpdatePartials(
				beagleInst,
				operationTuple,
				operationCount,
				cumulativeScaleIndex),
			"beagleUpdatePartials");

		//to indicate that the partial is now clean, fill the holder that represents it, despite the fact that the memory there is never being used
		//DEBUG - need to figure out second argument here (direction) which I think determined how likely a holder is to be recycled.
		claMan->FillHolder(theOp->destCLAIndex, 0);
#endif
		}
	}

void CalculationManager::PerformClaOperationBatch(const list<ClaOperation> &theOps){
		//not sure if this is right - will always use a single scale array for destWrite (essentially
		//scratch space, I think) and then pass a cumulative scaler to actually keep track of the scaling
		int destinationScaleWrite = (rescaleBeagle ? claMan->NumClas() : BEAGLE_OP_NONE);
		//need to work out batch rescaling yet
		assert(!rescaleBeagle);
		int	destinationScaleRead = BEAGLE_OP_NONE;
/*
		int operationTuple[7] = {PartialIndexForBeagle(theOp->destCLAIndex),
                                destinationScaleWrite,
                                destinationScaleRead,
								PartialIndexForBeagle(theOp->childCLAIndex1),
 								PmatIndexForBeagle(theOp->transMatIndex1),
 								PartialIndexForBeagle(theOp->childCLAIndex2),
								PmatIndexForBeagle(theOp->transMatIndex2)};
*/
		vector<int> tupleList;
		for(list<ClaOperation>::const_iterator it = theOps.begin();it != theOps.end();it++){
			tupleList.push_back(PartialIndexForBeagle((*it).destCLAIndex));
            tupleList.push_back(destinationScaleWrite);
            tupleList.push_back(destinationScaleRead);
			tupleList.push_back(PartialIndexForBeagle((*it).childCLAIndex1));
 			tupleList.push_back(PmatIndexForBeagle((*it).transMatIndex1));
 			tupleList.push_back(PartialIndexForBeagle((*it).childCLAIndex2));
			tupleList.push_back(PmatIndexForBeagle((*it).transMatIndex2));
			}
/*
		outman.DebugMessageNoCR("PARTS\t");
		outman.DebugMessageNoCR("\tD\t%d (%d)", PartialIndexForBeagle(theOp->destCLAIndex), theOp->destCLAIndex);
		outman.DebugMessageNoCR("\tC\t%d (%d)\t%d (%d)", PartialIndexForBeagle(theOp->childCLAIndex1), theOp->childCLAIndex1, PartialIndexForBeagle(theOp->childCLAIndex2), theOp->childCLAIndex2);
		outman.DebugMessage("\tP\t%d (%d)\t%d (%d)", PmatIndexForBeagle(theOp->transMatIndex1), theOp->transMatIndex1, PmatIndexForBeagle(theOp->transMatIndex2), theOp->transMatIndex2);
		
		int instanceCount = 1;
		int operationCount = 1;
*/
		int instanceCount = 1;
		int operationCount = theOps.size();

		//accumulate rescaling factors - For scale arrays my indexing scheme and Beagle's happen to be the same, 
		//and my negative (tip) corresponds to a NULL in beagle
		int cumulativeScaleIndex = BEAGLE_OP_NONE;
/*		if(rescaleBeagle){
			cumulativeScaleIndex = theOp->destCLAIndex;
			AccumulateRescalers(cumulativeScaleIndex, theOp->childCLAIndex1, theOp->childCLAIndex2);
			}
*/
		CheckBeagleReturnValue(
			beagleUpdatePartials(
				beagleInst,
				&tupleList[0],
				operationCount,
				cumulativeScaleIndex),
			"beagleUpdatePartials");

		//to indicate that the partial is now clean, fill the holder that represents it, despite the fact that the memory there is never being used
		//DEBUG - need to figure out second argument here (direction) which I think determined how likely a holder is to be recycled.
		for(list<ClaOperation>::const_iterator it = theOps.begin();it != theOps.end();it++){
			claMan->FillHolder((*it).destCLAIndex, 0);
			}
	}

#ifdef USE_BEAGLE
//For scale arrays my indexing scheme and Beagle's happen to be the same
void CalculationManager::AccumulateRescalers(int destIndex, int childIndex1, int childIndex2){
	//further accumulate the cumulative rescalings from the two children
	
	//clear out the destination scaler first
	CheckBeagleReturnValue(
		beagleResetScaleFactors(beagleInst,
			destIndex), 
		"beagleResetScaleFactors");

	vector<int> childScalers;
	if( ! (childIndex1 < 0))
		childScalers.push_back(childIndex1);
	if( ! (childIndex2 < 0))
		childScalers.push_back(childIndex2);

	if(childScalers.size() > 0){
		CheckBeagleReturnValue(
			beagleAccumulateScaleFactors(beagleInst,
				&(childScalers[0]),
				childScalers.size(),
				destIndex),
			"beagleAccumulateScaleFactors");
		}
	}
#endif

//need to change this to allow multiple simultaneous ops
void CalculationManager::PerformTransMatOperation(const TransMatOperation *theOp){
#ifndef USE_BEAGLE
	//mod->AltCalcPmat(theOp->edgeLength, pmatMan->GetPmat(theOp->destTransMatIndex)->theMat);
	//this is a bit ridiculous, but these funcs won't really be getting used except with beagle
	pmatMan->GetMutableHolder(theOp->destTransMatIndex)->myMod->AltCalcPmat(theOp->edgeLength, pmatMan->GetPmat(theOp->destTransMatIndex)->theMat);
#else 

//	outman.DebugMessage("**ENTERING PERFORM PMAT");
//	outman.DebugMessage("SETTING EIGEN");
	//DEBUG - currently just using a single eigen index and sending it every time
	int eigenIndex = 0;
	
	//this call will calculate the eigen solution first if necessary
	ModelEigenSolution sol;
	pmatMan->GetHolder(theOp->destTransMatIndex)->GetEigenSolution(sol);
	//not sure why I had this as mutable
	//pmatMan->GetMutableHolder(theOp->destTransMatIndex)->GetEigenSolution(sol);

	if(sol.changed){
		CheckBeagleReturnValue(
			beagleSetEigenDecomposition(
				beagleInst,
				eigenIndex,
				sol.eigenVecs,
				sol.invEigenVecs,
				sol.eigenVals),
			"beagleSetEigenDecomposition");
		}

	vector<double> categRates;
	pmatMan->GetHolder(theOp->destTransMatIndex)->GetCategoryRatesForBeagle(categRates);

//	outman.DebugMessage("SETTING CATEGORY RATES");
	if(categRates.size() > 0){
		CheckBeagleReturnValue(
			beagleSetCategoryRates(beagleInst,
				&(categRates[0])),
			"beagleSetCategoryRates");
		}

	//int pmatInd[2] = {PmatIndexForBeagle(theOp->transMatIndex1), PmatIndexForBeagle(theOp->transMatIndex2)};
	int pmatInd[1] = {PmatIndexForBeagle(theOp->destTransMatIndex)};
	int d1MatInd[1] = {D1MatIndexForBeagle(theOp->destTransMatIndex)};
	int d2MatInd[1] = {D2MatIndexForBeagle(theOp->destTransMatIndex)};

	//don't need to include the blen multiplier here, since it has already been used to scale the eigen values (not quite same as scaling blen)
	double edgeLens[1] = {pmatMan->GetHolder(theOp->destTransMatIndex)->edgeLen};
	int count = 1;

	//DEBUG
	//outman.UserMessage("edgelen = %e", edgeLens[0]);

	//outman.DebugMessageNoCR("UPDATING TRANS MAT ");
	//outman.DebugMessage("%d (%d), eigen %d, blen %f", PmatIndexForBeagle(theOp->destTransMatIndex), theOp->destTransMatIndex, eigenIndex, edgeLens[0]);
	CheckBeagleReturnValue(
		beagleUpdateTransitionMatrices(
			beagleInst,
			eigenIndex,
            pmatInd,
			(theOp->calcDerivs ? d1MatInd : NULL),
			(theOp->calcDerivs ? d2MatInd : NULL),
			edgeLens,
            count),
		"beagleUpdateTransitionMatrices");
#endif

#ifdef OUTPUT_PMATS
if(theOp->calcDerivs){
	int nrates = pmatMan->GetCorrespondingModel(theOp->destTransMatIndex)->NumRateCatsForBeagle();
	int nstates = data->NStates();
	vector<double> outMat(nstates * nstates * nrates);
	beagleGetTransitionMatrix(beagleInst, D1MatIndexForBeagle(theOp->destTransMatIndex), &(outMat[0]));
	vector<double>::iterator it = outMat.begin();
	for(int r = 0;r < nrates;r++){
		for(int i = 0;i < nstates;i++){
			for(int j=0;j < nstates;j++){
				outman.UserMessageNoCR("%.6f\t", *it++);
				}
			outman.UserMessage("");
			}
		}
	}		

#endif
	}

void CalculationManager::PerformTransMatOperationBatch(const list<TransMatOperation> &theOps){
#ifndef USE_BEAGLE
	//mod->AltCalcPmat(theOp->edgeLength, pmatMan->GetPmat(theOp->destTransMatIndex)->theMat);
	//this is a bit ridiculous, but these funcs won't really be getting used except with beagle
	pmatMan->GetMutableHolder(theOp->destTransMatIndex)->myMod->AltCalcPmat(theOp->edgeLength, pmatMan->GetPmat(theOp->destTransMatIndex)->theMat);
#else 
	//for now assuming that all transmat ops in a set have same eigen solution and rate multipliers

//	outman.DebugMessage("**ENTERING PERFORM PMAT");
//	outman.DebugMessage("SETTING EIGEN");
	//DEBUG - currently just using a single eigen index and sending it every time
	int eigenIndex = 0;
	
	//this call will calculate the eigen solution first if necessary
	ModelEigenSolution sol;
	pmatMan->GetHolder(theOps.begin()->destTransMatIndex)->GetEigenSolution(sol);
	//not sure why I had this as mutable
	//pmatMan->GetMutableHolder(theOp->destTransMatIndex)->GetEigenSolution(sol);

	if(sol.changed){
		CheckBeagleReturnValue(
			beagleSetEigenDecomposition(
				beagleInst,
				eigenIndex,
				sol.eigenVecs,
				sol.invEigenVecs,
				sol.eigenVals),
			"beagleSetEigenDecomposition");
		}

	vector<FLOAT_TYPE> categRates;
	pmatMan->GetHolder(theOps.begin()->destTransMatIndex)->GetCategoryRatesForBeagle(categRates);

//	outman.DebugMessage("SETTING CATEGORY RATES");
	if(categRates.size() > 0){
		CheckBeagleReturnValue(
			beagleSetCategoryRates(beagleInst,
				&(categRates[0])),
			"beagleSetCategoryRates");
		}


	vector<int> pmatInd;
	vector<int> d1MatInd;
	vector<int> d2MatInd;
	vector<double> edgeLens;

	for(list<TransMatOperation>::const_iterator it = theOps.begin();it != theOps.end();it++){
		pmatInd.push_back(PmatIndexForBeagle((*it).destTransMatIndex));
		d1MatInd.push_back(D1MatIndexForBeagle((*it).destTransMatIndex));
		d2MatInd.push_back(D2MatIndexForBeagle((*it).destTransMatIndex));
		edgeLens.push_back(pmatMan->GetHolder((*it).destTransMatIndex)->edgeLen);
		}

	//don't need to include the blen multiplier here, since it has already been used to scale the eigen values (not quite same as scaling blen)
//	double edgeLens[1] = {pmatMan->GetHolder(theOp->destTransMatIndex)->edgeLen};
	//int count = 1;
	int count = pmatInd.size();
	bool calcDerivs = theOps.begin()->calcDerivs;

	//DEBUG
	//outman.UserMessage("edgelen = %e", edgeLens[0]);

	//outman.DebugMessageNoCR("UPDATING TRANS MAT ");
	//outman.DebugMessage("%d (%d), eigen %d, blen %f", PmatIndexForBeagle(theOp->destTransMatIndex), theOp->destTransMatIndex, eigenIndex, edgeLens[0]);
	CheckBeagleReturnValue(
		beagleUpdateTransitionMatrices(
			beagleInst,
			eigenIndex,
            &pmatInd[0],
			(calcDerivs ? &d1MatInd[0] : NULL),
			(calcDerivs ? &d2MatInd[0] : NULL),
			&edgeLens[0],
            count),
		"beagleUpdateTransitionMatrices");
#endif

#ifdef OUTPUT_PMATS
	int nstates = data->NStates();
	
	int nrates = pmatMan->GetCorrespondingModel(theOps.begin()->destTransMatIndex)->NumRateCatsForBeagle();

	if(calcDerivs){

	vector<double> outMat(nstates * nstates * nrates);
	for(list<TransMatOperation>::const_iterator pit = theOps.begin();pit != theOps.end();pit++){				
		outMat.clear();
		for(int p=0;p<nstates*nstates*nrates;p++) outMat.push_back(p);
		//beagleGetTransitionMatrix(beagleInst, PmatIndexForBeagle((*pit).destTransMatIndex), &(outMat[0]));
		beagleGetTransitionMatrix(beagleInst, D1MatIndexForBeagle((*pit).destTransMatIndex), &(outMat[0]));
		vector<double>::iterator it = outMat.begin();
		outman.UserMessage("%d", (*pit).destTransMatIndex);
		for(int r = 0;r < nrates;r++){
			for(int i = 0;i < nstates;i++){
				for(int j=0;j < nstates;j++){
					outman.UserMessageNoCR("%.6f\t", *it++);
					}
				outman.UserMessage("");
				}
			}
		}
		}
#endif
	}

ScoreSet CalculationManager::PerformScoringOperation(const ScoringOperation *theOp){
	ScoreSet results = {0.0, 0.0, 0.0};
	//outman.DebugMessage("*ENTERING PERFORM SCORING");

	if(!useBeagle){
		results.lnL = GetScorePartialInternalRateHet(claMan->GetCla(theOp->childClaIndex1), 
			claMan->GetCla(theOp->childClaIndex2), 
			pmatMan->GetPmat(theOp->transMatIndex),
			pmatMan->GetCorrespondingModel(theOp->transMatIndex));
		}

	else{
#ifdef USE_BEAGLE

		int buffer1[1] = {PartialIndexForBeagle(theOp->childClaIndex1)};
		int buffer2[1] = {PartialIndexForBeagle(theOp->childClaIndex2)};
		int numInput = 2;

		//arrays to hold site lnLs and derivs returned from Beagle
		vector<double> siteLikesOut(data->NChar());
		vector<double> siteD1Out(data->NChar());
		vector<double> siteD2Out(data->NChar());

		//state freqs
		vector<FLOAT_TYPE> freqs(pmatMan->GetNumStates());
		pmatMan->GetCorrespondingModel(theOp->transMatIndex)->GetStateFreqs(&(freqs[0]));
		int freqIndex = 0;
		CheckBeagleReturnValue(
			beagleSetStateFrequencies(
				beagleInst,
				freqIndex,
				&(freqs[0])),
			"beagleSetStateFrequencies");

		//category weights (or probs)
		vector<FLOAT_TYPE> inWeights;
		pmatMan->GetHolder(theOp->transMatIndex)->GetCategoryWeightsForBeagle(inWeights);
		int weightIndex = 0;
		if(inWeights.size() > 0){
			CheckBeagleReturnValue(
				beagleSetCategoryWeights(
					beagleInst,
					weightIndex,
					&(inWeights[0])),
				"beagleSetCategoryWeights");
			}

		int pmatIndeces[1] = {PmatIndexForBeagle(theOp->transMatIndex)};
		int	d1MatIndeces[1] = {D1MatIndexForBeagle(theOp->transMatIndex)};
		int	d2MatIndeces[1] = {D2MatIndexForBeagle(theOp->transMatIndex)};

		outman.DebugMessageNoCR("Scr:\t");
		outman.DebugMessageNoCR("\tC\t%d (%d)\t%d (%d)", buffer2[0], theOp->childClaIndex2, buffer1[0], theOp->childClaIndex1);
		outman.DebugMessage("\tP\t%d (%d)", pmatIndeces[0], theOp->transMatIndex);

		int count = 1;

		//accumulate rescaling factors of the two clas that are being combined (one might be a tip)
		//For scale arrays my indexing scheme and Beagle's happen to be the same, and my negative (tip) corresponds to a NULL in beagle
		//use this scratch index to hold the final accumulated scalers 
		int cumulativeScaleIndex = BEAGLE_OP_NONE;
		if(rescaleBeagle){
			cumulativeScaleIndex = claMan->NumClas();
			AccumulateRescalers(cumulativeScaleIndex, theOp->childClaIndex1, theOp->childClaIndex2);
			}

		CheckBeagleReturnValue(
			beagleCalculateEdgeLogLikelihoods(
				beagleInst,
				buffer2,
				buffer1,
				pmatIndeces,
				(theOp->derivatives ? d1MatIndeces : NULL),
				(theOp->derivatives ? d2MatIndeces : NULL),
				&weightIndex,
				&freqIndex,
				/*
				&(inWeights[0]),
				&(freqs[0]),
				*/
				&cumulativeScaleIndex,
				count,
				&siteLikesOut[0],
				(theOp->derivatives ? &siteD1Out[0] : NULL),
				(theOp->derivatives ? &siteD2Out[0] : NULL)),
			"beagleCalculateEdgeLogLikelihoods");

		bool beagleReturnsSums = true;

		if(beagleReturnsSums){
/*			beagleGetSiteLogLikelihoods(beagleInst, &(siteLikesOut[0]));
			if(theOp->derivatives)
				beagleGetSiteDerivatives(beagleInst, &(siteD1Out[0]), &(siteD2Out[0]));
			results = SumSiteValues(&siteLikesOut[0], (theOp->derivatives ? &siteD1Out[0] : NULL), (theOp->derivatives ? &siteD2Out[0] : NULL));
*/

			results.lnL = siteLikesOut[0];
			results.d1 = siteD1Out[0];
			results.d2 = siteD2Out[0];
			}
		else
			results = SumSiteValues(&siteLikesOut[0], (theOp->derivatives ? &siteD1Out[0] : NULL), (theOp->derivatives ? &siteD2Out[0] : NULL));
		assert(results.lnL < 0.0 && results.lnL > -10.0e10);
		assert(results.d1 < 10.0e15 && results.d1 > -10.0e15);
		assert(results.d2 < 10.0e25 && results.d2 > -10.0e25);
		//DEBUG
		outman.DebugMessage("L %f D1 %f D2 %f", results.lnL, results.d1, results.d2);
#endif
		}
	return results;
	}

//these are my old functions, not used when in beagle mode
void CalculationManager::CalcFullClaInternalInternal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const TransMat *Lpmat, const TransMat *Rpmat, const Model *mod){

	FLOAT_TYPE *Lpr = **Lpmat->theMat;
	FLOAT_TYPE *Rpr = **Rpmat->theMat;

	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	FLOAT_TYPE *dest=destCLA->arr;
	const FLOAT_TYPE *LCL=LCLA->arr;
	const FLOAT_TYPE *RCL=RCLA->arr;
	FLOAT_TYPE L1, L2, L3, L4, R1, R2, R3, R4;

	const int nRateCats = mod->NRateCats();
	const int nchar = data->NChar();
	const int *counts = data->GetCounts();

#ifdef CUDA_GPU
	if (cudaman->GetGPUCLAEnabled()) {
		cudaman->ComputeGPUCLA(Lpr, Rpr, LCL, RCL, dest);
	} else {
#endif

#ifdef UNIX
	madvise(dest, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise((void *)LCL, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise((void *)RCL, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	if(nRateCats == 4){//the unrolled 4 rate version
#ifdef OMP_INTINTCLA
		#pragma omp parallel for private(dest, LCL, RCL, L1, L2, L3, L4, R1, R2, R3, R4)
		for(int i=0;i<nchar;i++){
			int index=4*4*i;
			dest = &(destCLA->arr[index]);
			LCL = &(LCLA->arr[index]);
			RCL= &(RCLA->arr[index]);
#else
		for(int i=0;i<nchar;i++){
#endif
#ifdef USE_COUNTS_IN_BOOT
			if(counts[i]> 0){
#else
			if(1){
#endif
				L1=((Lpr[0]*LCL[0])+(Lpr[1]*LCL[1]))+((Lpr[2]*LCL[2])+(Lpr[3]*LCL[3]));
				L2=((Lpr[4]*LCL[0])+(Lpr[5]*LCL[1]))+((Lpr[6]*LCL[2])+(Lpr[7]*LCL[3]));
				L3=((Lpr[8]*LCL[0])+(Lpr[9]*LCL[1]))+((Lpr[10]*LCL[2])+(Lpr[11]*LCL[3]));
				L4=((Lpr[12]*LCL[0])+(Lpr[13]*LCL[1]))+((Lpr[14]*LCL[2])+(Lpr[15]*LCL[3]));

				R1=((Rpr[0]*RCL[0])+(Rpr[1]*RCL[1]))+((Rpr[2]*RCL[2])+(Rpr[3]*RCL[3]));
				R2=((Rpr[4]*RCL[0])+(Rpr[5]*RCL[1]))+((Rpr[6]*RCL[2])+(Rpr[7]*RCL[3]));
				R3=((Rpr[8]*RCL[0])+(Rpr[9]*RCL[1]))+((Rpr[10]*RCL[2])+(Rpr[11]*RCL[3]));
				R4=((Rpr[12]*RCL[0])+(Rpr[13]*RCL[1]))+((Rpr[14]*RCL[2])+(Rpr[15]*RCL[3]));

				dest[0] = L1 * R1;
				dest[1] = L2 * R2;
				dest[2] = L3 * R3;
				dest[3] = L4 * R4;

				dest+=4;
				LCL+=4;
				RCL+=4;

				L1=(Lpr[16+0]*LCL[0]+Lpr[16+1]*LCL[1])+(Lpr[16+2]*LCL[2]+Lpr[16+3]*LCL[3]);
				L2=(Lpr[16+4]*LCL[0]+Lpr[16+5]*LCL[1])+(Lpr[16+6]*LCL[2]+Lpr[16+7]*LCL[3]);
				L3=(Lpr[16+8]*LCL[0]+Lpr[16+9]*LCL[1])+(Lpr[16+10]*LCL[2]+Lpr[16+11]*LCL[3]);
				L4=(Lpr[16+12]*LCL[0]+Lpr[16+13]*LCL[1])+(Lpr[16+14]*LCL[2]+Lpr[16+15]*LCL[3]);

				R1=(Rpr[16+0]*RCL[0]+Rpr[16+1]*RCL[1])+(Rpr[16+2]*RCL[2]+Rpr[16+3]*RCL[3]);
				R2=(Rpr[16+4]*RCL[0]+Rpr[16+5]*RCL[1])+(Rpr[16+6]*RCL[2]+Rpr[16+7]*RCL[3]);
				R3=(Rpr[16+8]*RCL[0]+Rpr[16+9]*RCL[1])+(Rpr[16+10]*RCL[2]+Rpr[16+11]*RCL[3]);
				R4=(Rpr[16+12]*RCL[0]+Rpr[16+13]*RCL[1])+(Rpr[16+14]*RCL[2]+Rpr[16+15]*RCL[3]);

				dest[0] = L1 * R1;
				dest[1] = L2 * R2;
				dest[2] = L3 * R3;
				dest[3] = L4 * R4;

				dest+=4;
				LCL+=4;
				RCL+=4;

				L1=(Lpr[32+0]*LCL[0]+Lpr[32+1]*LCL[1])+(Lpr[32+2]*LCL[2]+Lpr[32+3]*LCL[3]);
				L2=(Lpr[32+4]*LCL[0]+Lpr[32+5]*LCL[1])+(Lpr[32+6]*LCL[2]+Lpr[32+7]*LCL[3]);
				L3=(Lpr[32+8]*LCL[0]+Lpr[32+9]*LCL[1])+(Lpr[32+10]*LCL[2]+Lpr[32+11]*LCL[3]);
				L4=(Lpr[32+12]*LCL[0]+Lpr[32+13]*LCL[1])+(Lpr[32+14]*LCL[2]+Lpr[32+15]*LCL[3]);

				R1=(Rpr[32+0]*RCL[0]+Rpr[32+1]*RCL[1])+(Rpr[32+2]*RCL[2]+Rpr[32+3]*RCL[3]);
				R2=(Rpr[32+4]*RCL[0]+Rpr[32+5]*RCL[1])+(Rpr[32+6]*RCL[2]+Rpr[32+7]*RCL[3]);
				R3=(Rpr[32+8]*RCL[0]+Rpr[32+9]*RCL[1])+(Rpr[32+10]*RCL[2]+Rpr[32+11]*RCL[3]);
				R4=(Rpr[32+12]*RCL[0]+Rpr[32+13]*RCL[1])+(Rpr[32+14]*RCL[2]+Rpr[32+15]*RCL[3]);

				dest[0] = L1 * R1;
				dest[1] = L2 * R2;
				dest[2] = L3 * R3;
				dest[3] = L4 * R4;

				dest+=4;
				LCL+=4;
				RCL+=4;

				L1=(Lpr[48+0]*LCL[0]+Lpr[48+1]*LCL[1])+(Lpr[48+2]*LCL[2]+Lpr[48+3]*LCL[3]);
				L2=(Lpr[48+4]*LCL[0]+Lpr[48+5]*LCL[1])+(Lpr[48+6]*LCL[2]+Lpr[48+7]*LCL[3]);
				L3=(Lpr[48+8]*LCL[0]+Lpr[48+9]*LCL[1])+(Lpr[48+10]*LCL[2]+Lpr[48+11]*LCL[3]);
				L4=(Lpr[48+12]*LCL[0]+Lpr[48+13]*LCL[1])+(Lpr[48+14]*LCL[2]+Lpr[48+15]*LCL[3]);

				R1=(Rpr[48+0]*RCL[0]+Rpr[48+1]*RCL[1])+(Rpr[48+2]*RCL[2]+Rpr[48+3]*RCL[3]);
				R2=(Rpr[48+4]*RCL[0]+Rpr[48+5]*RCL[1])+(Rpr[48+6]*RCL[2]+Rpr[48+7]*RCL[3]);
				R3=(Rpr[48+8]*RCL[0]+Rpr[48+9]*RCL[1])+(Rpr[48+10]*RCL[2]+Rpr[48+11]*RCL[3]);
				R4=(Rpr[48+12]*RCL[0]+Rpr[48+13]*RCL[1])+(Rpr[48+14]*RCL[2]+Rpr[48+15]*RCL[3]);

				dest[0] = L1 * R1;
				dest[1] = L2 * R2;
				dest[2] = L3 * R3;
				dest[3] = L4 * R4;

				dest+=4;
				LCL+=4;
				RCL+=4;

#ifdef ALLOW_SINGLE_SITE
				if(siteToScore > -1) break;
#endif
				}
			}
		}

	else{//the general N rate version
		int r;
#ifdef OMP_INTINTCLA
		int index;
		#pragma omp parallel for private(r, index, dest, LCL, RCL, L1, L2, L3, L4, R1, R2, R3, R4)
		for(int i=0;i<nchar;i++) {
			index=4*nRateCats*i;
			dest = &(destCLA->arr[index]);
			LCL = &(LCLA->arr[index]);
			RCL= &(RCLA->arr[index]);
#else
		for(int i=0;i<nchar;i++) {
#endif
#ifdef USE_COUNTS_IN_BOOT
			if(counts[i] > 0){
#else
			if(1){
#endif
				for(r=0;r<nRateCats;r++){
					L1=( Lpr[16*r+0]*LCL[0]+Lpr[16*r+1]*LCL[1]+Lpr[16*r+2]*LCL[2]+Lpr[16*r+3]*LCL[3]);
					L2=( Lpr[16*r+4]*LCL[0]+Lpr[16*r+5]*LCL[1]+Lpr[16*r+6]*LCL[2]+Lpr[16*r+7]*LCL[3]);
					L3=( Lpr[16*r+8]*LCL[0]+Lpr[16*r+9]*LCL[1]+Lpr[16*r+10]*LCL[2]+Lpr[16*r+11]*LCL[3]);
					L4=( Lpr[16*r+12]*LCL[0]+Lpr[16*r+13]*LCL[1]+Lpr[16*r+14]*LCL[2]+Lpr[16*r+15]*LCL[3]);

					R1=(Rpr[16*r+0]*RCL[0]+Rpr[16*r+1]*RCL[1]+Rpr[16*r+2]*RCL[2]+Rpr[16*r+3]*RCL[3]);
					R2=(Rpr[16*r+4]*RCL[0]+Rpr[16*r+5]*RCL[1]+Rpr[16*r+6]*RCL[2]+Rpr[16*r+7]*RCL[3]);
					R3=(Rpr[16*r+8]*RCL[0]+Rpr[16*r+9]*RCL[1]+Rpr[16*r+10]*RCL[2]+Rpr[16*r+11]*RCL[3]);
					R4=(Rpr[16*r+12]*RCL[0]+Rpr[16*r+13]*RCL[1]+Rpr[16*r+14]*RCL[2]+Rpr[16*r+15]*RCL[3]);

					dest[0] = L1 * R1;
					dest[1] = L2 * R2;
					dest[2] = L3 * R3;
					dest[3] = L4 * R4;

					dest+=4;
					LCL+=4;
					RCL+=4;
					}
#ifdef ALLOW_SINGLE_SITE
				if(siteToScore > -1) break;
#endif
				}
			}
		}

#ifdef CUDA_GPU
	}
#endif

	const int *left_mult=LCLA->underflow_mult;
	const int *right_mult=RCLA->underflow_mult;
	int *undermult=destCLA->underflow_mult;

	for(int i=0;i<nchar;i++){
		undermult[i] = left_mult[i] + right_mult[i];
		}
	destCLA->rescaleRank = 2 + LCLA->rescaleRank + RCLA->rescaleRank;
	}

void CalculationManager::CalcFullCLATerminalTerminal(CondLikeArray *destCLA, const char *Ldata, const char *Rdata, const TransMat *Lpmat, const TransMat* Rpmat, const Model *mod){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.

	FLOAT_TYPE *Lpr = **Lpmat->theMat;
	FLOAT_TYPE *Rpr = **Rpmat->theMat;

	FLOAT_TYPE *dest=destCLA->arr;

	const int nRateCats = mod->NRateCats();
	const int nchar = data->NChar();
	const int *counts = data->GetCounts();

#ifdef UNIX
	madvise(dest, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

#ifdef ALLOW_SINGLE_SITE
	if(siteToScore > 0){
		Ldata = AdvanceDataPointer(Ldata, siteToScore);
		Rdata = AdvanceDataPointer(Rdata, siteToScore);
		}
#endif

	for(int i=0;i<nchar;i++){
#ifdef USE_COUNTS_IN_BOOT
		if(counts[i] > 0){
#else
		if(1){
#endif
			if(*Ldata > -1 && *Rdata > -1){
				for(int r=0;r<nRateCats;r++){
					*(dest++) = Lpr[(*Ldata)+16*r] * Rpr[(*Rdata)+16*r];
					*(dest++) = Lpr[(*Ldata+4)+16*r] * Rpr[(*Rdata+4)+16*r];
					*(dest++) = Lpr[(*Ldata+8)+16*r] * Rpr[(*Rdata+8)+16*r];
					*(dest++) = Lpr[(*Ldata+12)+16*r] * Rpr[(*Rdata+12)+16*r];
					}
				Ldata++;
				Rdata++;
				}

			else if((*Ldata == -4 && *Rdata == -4) || (*Ldata == -4 && *Rdata > -1) || (*Rdata == -4 && *Ldata > -1)){//total ambiguity of left, right or both

				if(*Ldata == -4 && *Rdata == -4) //total ambiguity of both
					for(int i=0;i< (4*nRateCats);i++) *(dest++) = ONE_POINT_ZERO;

				else if(*Ldata == -4){//total ambiguity of left
					for(int i=0;i<nRateCats;i++){
						*(dest++) = Rpr[(*Rdata)+16*i];
						*(dest++) = Rpr[(*Rdata+4)+16*i];
						*(dest++) = Rpr[(*Rdata+8)+16*i];
						*(dest++) = Rpr[(*Rdata+12)+16*i];
						assert(*(dest-4)>=ZERO_POINT_ZERO);
						}
					}
				else{//total ambiguity of right
					for(int i=0;i<nRateCats;i++){
						*(dest++) = Lpr[(*Ldata)+16*i];
						*(dest++) = Lpr[(*Ldata+4)+16*i];
						*(dest++) = Lpr[(*Ldata+8)+16*i];
						*(dest++) = Lpr[(*Ldata+12)+16*i];
						assert(*(dest-4)>=ZERO_POINT_ZERO);
						}
					}
				Ldata++;
				Rdata++;
				}
			else {//partial ambiguity of left, right or both
				if(*Ldata>-1){//unambiguous left
					for(int i=0;i<nRateCats;i++){
						*(dest+(i*4)) = Lpr[(*Ldata)+16*i];
						*(dest+(i*4)+1) = Lpr[(*Ldata+4)+16*i];
						*(dest+(i*4)+2) = Lpr[(*Ldata+8)+16*i];
						*(dest+(i*4)+3) = Lpr[(*Ldata+12)+16*i];
						assert(*(dest)>=ZERO_POINT_ZERO);
						}
					Ldata++;
					}
				else{
					if(*Ldata==-4){//fully ambiguous left
						for(int i=0;i< (4*nRateCats);i++){
							*(dest+i)=ONE_POINT_ZERO;
							}
						Ldata++;
						}

					else{//partially ambiguous left
						int nstates=-*(Ldata++);
						for(int q=0;q< (4*nRateCats);q++) dest[q]=0;
						for(int i=0;i<nstates;i++){
							for(int r=0;r<nRateCats;r++){
								*(dest+(r*4)) += Lpr[(*Ldata)+16*r];
								*(dest+(r*4)+1) += Lpr[(*Ldata+4)+16*r];
								*(dest+(r*4)+2) += Lpr[(*Ldata+8)+16*r];
								*(dest+(r*4)+3) += Lpr[(*Ldata+12)+16*r];
								}
							Ldata++;
							}
						}
					}
				if(*Rdata>-1){//unambiguous right
					for(int i=0;i<nRateCats;i++){
						*(dest++) *= Rpr[(*Rdata)+16*i];
						*(dest++) *= Rpr[(*Rdata+4)+16*i];
						*(dest++) *= Rpr[(*Rdata+8)+16*i];
						*(dest++) *= Rpr[(*Rdata+12)+16*i];
						}
					Rdata++;
					}
				else if(*Rdata != -4){//partially ambiguous right
					char nstates=-1 * *(Rdata++);
					//create a temporary cla to hold the results from the ambiguity of the right,
					//which need to be +'s
					//FLOAT_TYPE *tempcla=new FLOAT_TYPE[4*nRateCats];
					vector<FLOAT_TYPE> tempcla(4*nRateCats);
					for(int i=0;i<nstates;i++){
						for(int r=0;r<nRateCats;r++){
							tempcla[(r*4)]   += Rpr[(*Rdata)+16*r];
							tempcla[(r*4)+1] += Rpr[(*Rdata+4)+16*r];
							tempcla[(r*4)+2] += Rpr[(*Rdata+8)+16*r];
							tempcla[(r*4)+3] += Rpr[(*Rdata+12)+16*r];
							}
						Rdata++;
						}
					//Now multiply the temporary results against the already calced left
					for(int i=0;i<nRateCats;i++){
						*(dest++) *= tempcla[(i*4)];
						*(dest++) *= tempcla[(i*4)+1];
						*(dest++) *= tempcla[(i*4)+2];
						*(dest++) *= tempcla[(i*4)+3];
						}
					}
				else{//fully ambiguous right
					dest+=(4*nRateCats);
					Rdata++;
					}
				}
#ifdef ALLOW_SINGLE_SITE
			if(siteToScore > -1) break;
#endif
			}
		else{//if the count for this site is 0
#ifdef OPEN_MP
			//this is a little strange, but dest only needs to be advanced in the case of OMP
			//because sections of the CLAs corresponding to sites with count=0 are skipped
			//over in OMP instead of being eliminated
			dest += 4 * nRateCats;
#endif
			if(*Ldata > -1 || *Ldata == -4) Ldata++;
			else{
				int states = -1 * *Ldata;
				do{
					Ldata++;
					}while (states-- > 0);
				}
			if(*Rdata > -1 || *Rdata == -4) Rdata++;
			else{
				int states = -1 * *Rdata;
				do{
					Rdata++;
					}while (states-- > 0);
				}
			}
		}

		for(int site=0;site<nchar;site++){
			destCLA->underflow_mult[site]=0;
			}
	destCLA->rescaleRank=2;
	}

void CalculationManager::CalcFullCLAInternalTerminal(CondLikeArray *destCLA, const CondLikeArray *LCLA, char *dat2, const TransMat *Lpmat, const TransMat *Rpmat, const Model *mod, const unsigned *ambigMap){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.

	FLOAT_TYPE *pr1 = **Lpmat->theMat;
	FLOAT_TYPE *pr2 = **Rpmat->theMat;

	FLOAT_TYPE *des=destCLA->arr;
	FLOAT_TYPE *dest=des;
	const FLOAT_TYPE *CL=LCLA->arr;
	const FLOAT_TYPE *CL1=CL;
	const char *data2=dat2;

	const int nchar = data->NChar();
	const int nRateCats = mod->NRateCats();
	const int *counts = data->GetCounts();

#ifdef UNIX
	madvise(dest, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise((void*)CL1, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

#ifdef ALLOW_SINGLE_SITE
	if(siteToScore > 0) data2 = AdvanceDataPointer(data2, siteToScore);
#endif

	if(nRateCats==4){//unrolled 4 rate version
#ifdef OMP_INTTERMCLA
		#pragma omp parallel for private(dest, CL1, data2)
		for(int i=0;i<nchar;i++){
			dest=&des[4*4*i];
			CL1=&CL[4*4*i];
			data2=&dat2[ambigMap[i]];
#else
		for(int i=0;i<nchar;i++){
#endif
#ifdef USE_COUNTS_IN_BOOT
			if(counts[i] > 0){
#else
			if(1){
#endif
				if(*data2 > -1){ //no ambiguity
					dest[0] = ((pr1[0]*CL1[0]+pr1[1]*CL1[1])+(pr1[2]*CL1[2]+pr1[3]*CL1[3])) * pr2[*data2];
					dest[1] = ((pr1[4]*CL1[0]+pr1[5]*CL1[1])+(pr1[6]*CL1[2]+pr1[7]*CL1[3])) * pr2[*data2+4];
					dest[2] = ((pr1[8]*CL1[0]+pr1[9]*CL1[1])+(pr1[10]*CL1[2]+pr1[11]*CL1[3])) * pr2[*data2+8];
					dest[3] = ((pr1[12]*CL1[0]+pr1[13]*CL1[1])+(pr1[14]*CL1[2]+pr1[15]*CL1[3])) * pr2[*data2+12];

					dest[4] = ((pr1[16]*CL1[4]+pr1[17]*CL1[5])+(pr1[18]*CL1[6]+pr1[19]*CL1[7])) * pr2[*data2+16];
					dest[5] = ((pr1[20]*CL1[4]+pr1[21]*CL1[5])+(pr1[22]*CL1[6]+pr1[23]*CL1[7])) * pr2[*data2+4+16];
					dest[6] = ((pr1[24]*CL1[4]+pr1[25]*CL1[5])+(pr1[26]*CL1[6]+pr1[27]*CL1[7])) * pr2[*data2+8+16];
					dest[7] = ((pr1[28]*CL1[4]+pr1[29]*CL1[5])+(pr1[30]*CL1[6]+pr1[31]*CL1[7])) * pr2[*data2+12+16];

					dest[8] = ((pr1[32]*CL1[8]+pr1[33]*CL1[9])+(pr1[34]*CL1[10]+pr1[35]*CL1[11])) * pr2[*data2+32];
					dest[9] = ((pr1[36]*CL1[8]+pr1[37]*CL1[9])+(pr1[38]*CL1[10]+pr1[39]*CL1[11])) * pr2[*data2+4+32];
					dest[10] = ((pr1[40]*CL1[8]+pr1[41]*CL1[9])+(pr1[42]*CL1[10]+pr1[43]*CL1[11])) * pr2[*data2+8+32];
					dest[11] = ((pr1[44]*CL1[8]+pr1[45]*CL1[9])+(pr1[46]*CL1[10]+pr1[47]*CL1[11])) * pr2[*data2+12+32];

					dest[12] = ((pr1[48]*CL1[12]+pr1[49]*CL1[13])+(pr1[50]*CL1[14]+pr1[51]*CL1[15])) * pr2[*data2+48];
					dest[13] = ((pr1[52]*CL1[12]+pr1[53]*CL1[13])+(pr1[54]*CL1[14]+pr1[55]*CL1[15])) * pr2[*data2+4+48];
					dest[14] = ((pr1[56]*CL1[12]+pr1[57]*CL1[13])+(pr1[58]*CL1[14]+pr1[59]*CL1[15])) * pr2[*data2+8+48];
					dest[15] = ((pr1[60]*CL1[12]+pr1[61]*CL1[13])+(pr1[62]*CL1[14]+pr1[63]*CL1[15])) * pr2[*data2+12+48];

					dest+=16;
					data2++;
					}
				else if(*data2 == -4){//total ambiguity
					dest[0] = ( pr1[0]*CL1[0]+pr1[1]*CL1[1]+pr1[2]*CL1[2]+pr1[3]*CL1[3]);
					dest[1] = ( pr1[4]*CL1[0]+pr1[5]*CL1[1]+pr1[6]*CL1[2]+pr1[7]*CL1[3]);
					dest[2] = ( pr1[8]*CL1[0]+pr1[9]*CL1[1]+pr1[10]*CL1[2]+pr1[11]*CL1[3]);
					dest[3] = ( pr1[12]*CL1[0]+pr1[13]*CL1[1]+pr1[14]*CL1[2]+pr1[15]*CL1[3]);

					dest[4] = ( pr1[16]*CL1[4]+pr1[17]*CL1[5]+pr1[18]*CL1[6]+pr1[19]*CL1[7]);
					dest[5] = ( pr1[20]*CL1[4]+pr1[21]*CL1[5]+pr1[22]*CL1[6]+pr1[23]*CL1[7]);
					dest[6] = ( pr1[24]*CL1[4]+pr1[25]*CL1[5]+pr1[26]*CL1[6]+pr1[27]*CL1[7]);
					dest[7] = ( pr1[28]*CL1[4]+pr1[29]*CL1[5]+pr1[30]*CL1[6]+pr1[31]*CL1[7]);

					dest[8] = ( pr1[32]*CL1[8]+pr1[33]*CL1[9]+pr1[34]*CL1[10]+pr1[35]*CL1[11]);
					dest[9] = ( pr1[36]*CL1[8]+pr1[37]*CL1[9]+pr1[38]*CL1[10]+pr1[39]*CL1[11]);
					dest[10] = ( pr1[40]*CL1[8]+pr1[41]*CL1[9]+pr1[42]*CL1[10]+pr1[43]*CL1[11]);
					dest[11] = ( pr1[44]*CL1[8]+pr1[45]*CL1[9]+pr1[46]*CL1[10]+pr1[47]*CL1[11]);

					dest[12] = ( pr1[48]*CL1[12]+pr1[49]*CL1[13]+pr1[50]*CL1[14]+pr1[51]*CL1[15]);
					dest[13] = ( pr1[52]*CL1[12]+pr1[53]*CL1[13]+pr1[54]*CL1[14]+pr1[55]*CL1[15]);
					dest[14] = ( pr1[56]*CL1[12]+pr1[57]*CL1[13]+pr1[58]*CL1[14]+pr1[59]*CL1[15]);
					dest[15] = ( pr1[60]*CL1[12]+pr1[61]*CL1[13]+pr1[62]*CL1[14]+pr1[63]*CL1[15]);

					dest+=16;
					data2++;
					}
				else {//partial ambiguity
					//first figure in the ambiguous terminal
					int nstates=-1 * *(data2++);
					for(int j=0;j<16;j++) dest[j]=ZERO_POINT_ZERO;
					for(int s=0;s<nstates;s++){
						for(int r=0;r<4;r++){
							*(dest+(r*4)) += pr2[(*data2)+16*r];
							*(dest+(r*4)+1) += pr2[(*data2+4)+16*r];
							*(dest+(r*4)+2) += pr2[(*data2+8)+16*r];
							*(dest+(r*4)+3) += pr2[(*data2+12)+16*r];
							}
						data2++;
						}

					//now add the internal child
					*(dest++) *= ( pr1[0]*CL1[0]+pr1[1]*CL1[1]+pr1[2]*CL1[2]+pr1[3]*CL1[3]);
					*(dest++) *= ( pr1[4]*CL1[0]+pr1[5]*CL1[1]+pr1[6]*CL1[2]+pr1[7]*CL1[3]);
					*(dest++) *= ( pr1[8]*CL1[0]+pr1[9]*CL1[1]+pr1[10]*CL1[2]+pr1[11]*CL1[3]);
					*(dest++) *= ( pr1[12]*CL1[0]+pr1[13]*CL1[1]+pr1[14]*CL1[2]+pr1[15]*CL1[3]);

					*(dest++) *= ( pr1[16]*CL1[4]+pr1[17]*CL1[5]+pr1[18]*CL1[6]+pr1[19]*CL1[7]);
					*(dest++) *= ( pr1[20]*CL1[4]+pr1[21]*CL1[5]+pr1[22]*CL1[6]+pr1[23]*CL1[7]);
					*(dest++) *= ( pr1[24]*CL1[4]+pr1[25]*CL1[5]+pr1[26]*CL1[6]+pr1[27]*CL1[7]);
					*(dest++) *= ( pr1[28]*CL1[4]+pr1[29]*CL1[5]+pr1[30]*CL1[6]+pr1[31]*CL1[7]);

					*(dest++) *= ( pr1[32]*CL1[8]+pr1[33]*CL1[9]+pr1[34]*CL1[10]+pr1[35]*CL1[11]);
					*(dest++) *= ( pr1[36]*CL1[8]+pr1[37]*CL1[9]+pr1[38]*CL1[10]+pr1[39]*CL1[11]);
					*(dest++) *= ( pr1[40]*CL1[8]+pr1[41]*CL1[9]+pr1[42]*CL1[10]+pr1[43]*CL1[11]);
					*(dest++) *= ( pr1[44]*CL1[8]+pr1[45]*CL1[9]+pr1[46]*CL1[10]+pr1[47]*CL1[11]);

					*(dest++) *= ( pr1[48]*CL1[12]+pr1[49]*CL1[13]+pr1[50]*CL1[14]+pr1[51]*CL1[15]);
					*(dest++) *= ( pr1[52]*CL1[12]+pr1[53]*CL1[13]+pr1[54]*CL1[14]+pr1[55]*CL1[15]);
					*(dest++) *= ( pr1[56]*CL1[12]+pr1[57]*CL1[13]+pr1[58]*CL1[14]+pr1[59]*CL1[15]);
					*(dest++) *= ( pr1[60]*CL1[12]+pr1[61]*CL1[13]+pr1[62]*CL1[14]+pr1[63]*CL1[15]);
					}
				CL1+=16;
#ifdef ALLOW_SINGLE_SITE
				if(siteToScore > -1) break;
#endif
				}
			else{
				data2 = AdvanceDataPointer(data2, 1);
				}
			}
		}
	else{//general N rate version
#ifdef OMP_INTTERMCLA
		#pragma omp parallel for private(dest, CL1, data2)
		for(int i=0;i<nchar;i++){
			dest=&des[4*nRateCats*i];
			CL1=&CL[4*nRateCats*i];
			data2=&dat2[ambigMap[i]];
#else
		for(int i=0;i<nchar;i++){
#endif
#ifdef USE_COUNTS_IN_BOOT
			if(counts[i] > 0){
#else
			if(1){
#endif
				if(*data2 > -1){ //no ambiguity
					for(int r=0;r<nRateCats;r++){
						dest[0] = ( pr1[16*r+0]*CL1[4*r+0]+pr1[16*r+1]*CL1[4*r+1]+pr1[16*r+2]*CL1[4*r+2]+pr1[16*r+3]*CL1[4*r+3]) * pr2[(*data2)+16*r];
						dest[1] = ( pr1[16*r+4]*CL1[4*r+0]+pr1[16*r+5]*CL1[4*r+1]+pr1[16*r+6]*CL1[4*r+2]+pr1[16*r+7]*CL1[4*r+3]) * pr2[(*data2+4)+16*r];
						dest[2] = ( pr1[16*r+8]*CL1[4*r+0]+pr1[16*r+9]*CL1[4*r+1]+pr1[16*r+10]*CL1[4*r+2]+pr1[16*r+11]*CL1[4*r+3]) * pr2[(*data2+8)+16*r];
						dest[3] = ( pr1[16*r+12]*CL1[4*r+0]+pr1[16*r+13]*CL1[4*r+1]+pr1[16*r+14]*CL1[4*r+2]+pr1[16*r+15]*CL1[4*r+3]) * pr2[(*data2+12)+16*r];
						dest+=4;
						}
					data2++;
					}
				else if(*data2 == -4){//total ambiguity
					for(int r=0;r<nRateCats;r++){
						dest[0] = ( pr1[16*r+0]*CL1[4*r+0]+pr1[16*r+1]*CL1[4*r+1]+pr1[16*r+2]*CL1[4*r+2]+pr1[16*r+3]*CL1[4*r+3]);
						dest[1] = ( pr1[16*r+4]*CL1[4*r+0]+pr1[16*r+5]*CL1[4*r+1]+pr1[16*r+6]*CL1[4*r+2]+pr1[16*r+7]*CL1[4*r+3]);
						dest[2] = ( pr1[16*r+8]*CL1[4*r+0]+pr1[16*r+9]*CL1[4*r+1]+pr1[16*r+10]*CL1[4*r+2]+pr1[16*r+11]*CL1[4*r+3]);
						dest[3] = ( pr1[16*r+12]*CL1[4*r+0]+pr1[16*r+13]*CL1[4*r+1]+pr1[16*r+14]*CL1[4*r+2]+pr1[16*r+15]*CL1[4*r+3]);
						dest+=4;
						}
					data2++;
					}
				else {//partial ambiguity
					//first figure in the ambiguous terminal
					int nstates=-1 * *(data2++);
					for(int q=0;q<4*nRateCats;q++) dest[q]=0;
					for(int s=0;s<nstates;s++){
						for(int r=0;r<nRateCats;r++){
							*(dest+(r*4)) += pr2[(*data2)+16*r];
							*(dest+(r*4)+1) += pr2[(*data2+4)+16*r];
							*(dest+(r*4)+2) += pr2[(*data2+8)+16*r];
							*(dest+(r*4)+3) += pr2[(*data2+12)+16*r];
							}
						data2++;
						}
					//now add the internal child
					for(int r=0;r<nRateCats;r++){
						*(dest++) *= ( pr1[16*r+0]*CL1[4*r+0]+pr1[16*r+1]*CL1[4*r+1]+pr1[16*r+2]*CL1[4*r+2]+pr1[16*r+3]*CL1[4*r+3]);
						*(dest++) *= ( pr1[16*r+4]*CL1[4*r+0]+pr1[16*r+5]*CL1[4*r+1]+pr1[16*r+6]*CL1[4*r+2]+pr1[16*r+7]*CL1[4*r+3]);
						*(dest++) *= ( pr1[16*r+8]*CL1[4*r+0]+pr1[16*r+9]*CL1[4*r+1]+pr1[16*r+10]*CL1[4*r+2]+pr1[16*r+11]*CL1[4*r+3]);
						*(dest++) *= ( pr1[16*r+12]*CL1[4*r+0]+pr1[16*r+13]*CL1[4*r+1]+pr1[16*r+14]*CL1[4*r+2]+pr1[16*r+15]*CL1[4*r+3]);
						}
					}
				CL1 += 4*nRateCats;
#ifdef ALLOW_SINGLE_SITE
				if(siteToScore > -1) break;
#endif
				}
			else{
				data2 = AdvanceDataPointer(data2, 1);
				}
			}
		}

	for(int i=0;i<nchar;i++)
		destCLA->underflow_mult[i]=LCLA->underflow_mult[i];

	destCLA->rescaleRank=LCLA->rescaleRank+2;
	}


FLOAT_TYPE CalculationManager::GetScorePartialInternalRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const TransMat *pmat, const Model *mod){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.

	const FLOAT_TYPE *prmat = **pmat->theMat;

	const FLOAT_TYPE *CL1=childCLA->arr;
	const FLOAT_TYPE *partial=partialCLA->arr;
	const int *underflow_mult1=partialCLA->underflow_mult;
	const int *underflow_mult2=childCLA->underflow_mult;

	const int nchar=data->NChar();
	const int nRateCats=mod->NRateCats();

	const int *countit=data->GetCounts();

	const FLOAT_TYPE *rateProb=const_cast<Model *> (mod)->GetRateProbs();
	const FLOAT_TYPE prI=mod->PropInvar();
	const int lastConst=data->LastConstant();
	const int *conBases=data->GetConstStates();

	FLOAT_TYPE freqs[4];
	for(int i=0;i<4;i++) freqs[i]=mod->StateFreq(i);


#ifdef UNIX
	madvise((void*)partial, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise((void*)CL1, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	//DEBUG
	//FLOAT_TYPE siteL, totallnL=ZERO_POINT_ZERO, grandSumlnL=ZERO_POINT_ZERO, unscaledlnL;
	double siteL, totallnL=ZERO_POINT_ZERO, grandSumlnL=ZERO_POINT_ZERO, unscaledlnL;
	FLOAT_TYPE La, Lc, Lg, Lt;

	vector<double> siteLikes(nchar);

	for(int i=0;i<nchar;i++){
#ifdef USE_COUNTS_IN_BOOT
		if(countit[i] > 0){
#else
		if(1){
#endif
			La=Lc=Lg=Lt=ZERO_POINT_ZERO;
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
				FLOAT_TYPE btot=ZERO_POINT_ZERO;
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

			unscaledlnL = (log(siteL) - underflow_mult1[i] - underflow_mult2[i]);
			totallnL += (countit[i] * unscaledlnL);

#ifdef ALLOW_SINGLE_SITE
			if(siteToScore > -1) break;
#endif
			}
		else{
#ifdef OPEN_MP
			//this is a little strange, but the arrays only needs to be advanced in the case of OMP
			//because sections of the CLAs corresponding to sites with count=0 are skipped
			//over in OMP instead of being eliminated
			partial+=4*nRateCats;
			CL1+=4*nRateCats;
#endif
			}
#ifdef LUMP_LIKES
		if((i + 1) % LUMP_FREQ == 0){
			grandSumlnL += totallnL;
			totallnL = ZERO_POINT_ZERO;
			}
#endif
//DEBUG
/*
		if(sitelikeLevel != 0)
			siteLikes[i] = unscaledlnL;
*/		}

#ifdef LUMP_LIKES
	totallnL += grandSumlnL;
#endif
//DEBUG
/*
	if(sitelikeLevel != 0){
		OutputSiteLikelihoods(siteLikes, underflow_mult1, underflow_mult2);
		}
*/
	return totallnL;
	}
FLOAT_TYPE CalculationManager::GetScorePartialTerminalRateHet(const CondLikeArray *partialCLA, const TransMat *pmat, const char *Ldata, const Model *mod){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	const FLOAT_TYPE *prmat = **pmat->theMat;
	const FLOAT_TYPE *partial=partialCLA->arr;
	const int *underflow_mult=partialCLA->underflow_mult;
	const int nRateCats=mod->NRateCats();
	const int nchar=data->NChar();
	const int *countit=data->GetCounts();
	const FLOAT_TYPE *rateProb=const_cast<Model *> (mod)->GetRateProbs();
	const int lastConst=data->LastConstant();
	const int *conBases=data->GetConstStates();
	const FLOAT_TYPE prI=mod->PropInvar();
	FLOAT_TYPE freqs[4];
	for(int i=0;i<4;i++) freqs[i]=mod->StateFreq(i);
#ifdef UNIX
	madvise((void*)partial, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

#ifdef ALLOW_SINGLE_SITE
	if(siteToScore > 0) Ldata = AdvanceDataPointer(Ldata, siteToScore);
#endif

	FLOAT_TYPE siteL, totallnL=ZERO_POINT_ZERO, grandSumlnL=ZERO_POINT_ZERO, unscaledlnL;
	FLOAT_TYPE La, Lc, Lg, Lt;

	vector<double> siteLikes(nchar);

	for(int i=0;i<nchar;i++){
#ifdef USE_COUNTS_IN_BOOT
		if(countit[i] > 0){
#else
		if(1){
#endif
			La=Lc=Lg=Lt=ZERO_POINT_ZERO;
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
			unscaledlnL = (log(siteL) - underflow_mult[i]);
			totallnL += (countit[i] * unscaledlnL);

#ifdef ALLOW_SINGLE_SITE
			if(siteToScore > -1) break;
#endif
			}
		else{
#ifdef OPEN_MP
			//this is a little strange, but partial only needs to be advanced in the case of OMP
			//because sections of the CLAs corresponding to sites with count=0 are skipped
			//over in OMP instead of being eliminated
			partial += 4*nRateCats;
#endif
			if(*Ldata > -1 || *Ldata == -4) Ldata++;
			else{
				int states = -1 * *Ldata;
				do{
					Ldata++;
					}while (states-- > 0);
				}
			}
#ifdef LUMP_LIKES
		if((i + 1) % LUMP_FREQ == 0){
			grandSumlnL += totallnL;
			totallnL = ZERO_POINT_ZERO;
			}
#endif
//DEBUG
/*
			if(sitelikeLevel != 0)
				siteLikes[i] = unscaledlnL;
				*/
		}
#ifdef LUMP_LIKES
	totallnL += grandSumlnL;
#endif
//DEBUG
/*
	if(sitelikeLevel != 0){
		OutputSiteLikelihoods(siteLikes, underflow_mult, NULL);
		}
*/
	return totallnL;
	}
