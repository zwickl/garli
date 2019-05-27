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

//#undef FULL_BEAGLE_DEBUG
//#undef OUTPUT_OTHER_BEAGLE
//#define OUTPUT_OTHER_BEAGLE
//#define OUTPUT_PMATS
//#undef OUTPUT_PARTIALS
//#define OUTPUT_BEAGLE_SITELIKES

#ifdef FULL_BEAGLE_DEBUG
	#define	OUTPUT_PMATS
	#define OUTPUT_BEAGLE_SITELIKES
	#define OUTPUT_PARTIALS
	#define OUTPUT_OTHER_BEAGLE
#endif

#ifdef USE_BEAGLE

// print possible beagle resources
void CalculationManager::OutputBeagleResources() const{
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
		if((rList->list[i].supportFlags & req_flags) == req_flags)
			outman.UserMessage("\t\t(Meets beagle instance requirements)");
		else
			outman.UserMessage("\t\t(Does not meet beagle instance requirements)");
		if((rList->list[i].supportFlags & pref_flags) == pref_flags)
			outman.UserMessage("\t\t(Meets beagle instance preferences)");
		else
			outman.UserMessage("\t\t(Does not meet beagle instance preferences)");
		}
    outman.UserMessageNoCR("\n");
	}

//print actual beagle resources used by specific instance
void CalculationManager::ParseInstanceDetails(const BeagleInstanceDetails *det){
    outman.UserMessageNoCR("Instance details:\n");
	outman.UserMessageNoCR("\t\tnumber: %d\n", det->resourceNumber);
	outman.UserMessageNoCR("\t\tresource name: %s\n", det->resourceName);
	outman.UserMessageNoCR("\t\timplementation name: %s\n", det->implName);
	outman.UserMessageNoCR("\t\tFlags: ");
	
	actual_flags = det->flags;
	InterpretBeagleResourceFlags(actual_flags, actualBeagleFlags);
	outman.UserMessage("%s", actualBeagleFlags.c_str());
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
		mess = "Beagle was unable to find a suitable resource in ";

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

void CalculationManager::FillBeagleOptionsMaps(){
	//define the mappings
	flagToName.insert(pair<long, string>(BEAGLE_FLAG_PROCESSOR_CPU, "CPU"));
	nameToFlag.insert(pair<string, long>("CPU", BEAGLE_FLAG_PROCESSOR_CPU));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_PROCESSOR_GPU, "GPU"));
	nameToFlag.insert(pair<string, long>("GPU", BEAGLE_FLAG_PROCESSOR_GPU));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_PROCESSOR_GPU, "CUDA"));
	nameToFlag.insert(pair<string, long>("GPU", BEAGLE_FLAG_FRAMEWORK_CUDA));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_PROCESSOR_GPU, "OPENCL"));
	nameToFlag.insert(pair<string, long>("GPU", BEAGLE_FLAG_FRAMEWORK_OPENCL));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_PRECISION_DOUBLE, "DOUBLE"));
	nameToFlag.insert(pair<string, long>("DOUBLE", BEAGLE_FLAG_PRECISION_DOUBLE));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_PRECISION_SINGLE, "SINGLE"));
	nameToFlag.insert(pair<string, long>("SINGLE", BEAGLE_FLAG_PRECISION_SINGLE));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_SCALERS_LOG, "RESCALE"));
	nameToFlag.insert(pair<string, long>("RESCALE", BEAGLE_FLAG_SCALERS_LOG));	

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_VECTOR_SSE, "SSE"));
	nameToFlag.insert(pair<string, long>("SSE", BEAGLE_FLAG_VECTOR_SSE));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_VECTOR_SSE, "AVX"));
	nameToFlag.insert(pair<string, long>("SSE", BEAGLE_FLAG_VECTOR_AVX));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_THREADING_OPENMP, "OPENMP"));
	nameToFlag.insert(pair<string, long>("OPENMP", BEAGLE_FLAG_THREADING_OPENMP));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_THREADING_OPENMP, "THREADING"));
	nameToFlag.insert(pair<string, long>("OPENMP", BEAGLE_FLAG_THREADING_CPP));

#ifdef OUTPUT_OTHER_BEAGLE
	//we don't usually care about these flags
	flagToName.insert(pair<long, string>(BEAGLE_FLAG_COMPUTATION_ASYNCH, "ASYNCH"));
	nameToFlag.insert(pair<string, long>("ASYNCH", BEAGLE_FLAG_COMPUTATION_ASYNCH));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_COMPUTATION_SYNCH, "SYNCH"));
	nameToFlag.insert(pair<string, long>("SYNCH", BEAGLE_FLAG_COMPUTATION_SYNCH));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_EIGEN_COMPLEX, "COMPLEX"));
	nameToFlag.insert(pair<string, long>("COMPLEX", BEAGLE_FLAG_EIGEN_COMPLEX));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_EIGEN_REAL, "REAL"));
	nameToFlag.insert(pair<string, long>("REAL", BEAGLE_FLAG_EIGEN_REAL));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_SCALING_MANUAL, "MANUAL_SCALING"));
	nameToFlag.insert(pair<string, long>("MANUAL_SCALING", BEAGLE_FLAG_SCALING_MANUAL));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_SCALING_AUTO, "AUTO_SCALING"));
	nameToFlag.insert(pair<string, long>("AUTO_SCALING", BEAGLE_FLAG_SCALING_AUTO));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_SCALING_ALWAYS, "ALWAYS_SCALING"));
	nameToFlag.insert(pair<string, long>("ALWAYS_SCALING", BEAGLE_FLAG_SCALING_ALWAYS));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_SCALERS_RAW, "RAW_SCALERS"));
	nameToFlag.insert(pair<string, long>("RAW_SCALERS", BEAGLE_FLAG_SCALERS_RAW));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_PROCESSOR_FPGA, "FPGA"));
	nameToFlag.insert(pair<string, long>("FPGA", BEAGLE_FLAG_PROCESSOR_FPGA));

	flagToName.insert(pair<long, string>(BEAGLE_FLAG_PROCESSOR_CELL, "CELL"));
	nameToFlag.insert(pair<string, long>("CELL", BEAGLE_FLAG_PROCESSOR_CELL));
#endif

	invalidFlags.push_back(BEAGLE_FLAG_PRECISION_DOUBLE | BEAGLE_FLAG_PRECISION_SINGLE);
	invalidFlags.push_back(BEAGLE_FLAG_PROCESSOR_CPU | BEAGLE_FLAG_PROCESSOR_GPU);
	}

//flag mask allows for ignoring bits that have already be interpreted, for example req_flags trump pref_flags, and thus
//redundant pref_flags don't need to be dealt with 
long CalculationManager::ParseBeagleFlagString(string flagsString, long flagMask /*=0*/) const{
	long ret = 0;

	NxsString::to_upper(flagsString);
	stringstream s(flagsString);

	map<string, long>::const_iterator itToFlag;

	string flag;
	while(!s.eof()){
		s >> flag;
		if (flag.size()) {
			itToFlag = nameToFlag.find(flag);
			if (itToFlag == nameToFlag.end()) {
				outman.UserMessage("Ignoring unknown beagle option: %s", flag.c_str());
			}
			else
				ret |= (*itToFlag).second & ~flagMask;
		}
		}
	
	for(vector<long>::const_iterator it = invalidFlags.begin(); it != invalidFlags.end(); it++){
		if(((*it) & ret) == (*it)){
			string str;
			InterpretBeagleResourceFlags((*it), str);
			throw 
				ErrorException("INVALID BEAGLE FLAG COMBINATION\nThese flags cannot be specified together: %s", str.c_str());
			}
		}
	return ret;
	}

//flag mask allows for ignoring bits that have already be interpreted, for example req_flags trump pref_flags, and thus
//redundant pref_flags don't need to be dealt with  
void CalculationManager::InterpretBeagleResourceFlags(long flags, string &list, long flagMask /*=0*/) const{
	map<long, string>::const_iterator itToName;
	long bit = 1;

	do{
		if(bit & flags & ~flagMask){
			itToName = flagToName.find(bit);
			if (itToName == flagToName.end()){
#ifdef OUTPUT_OTHER_BEAGLE
				outman.DebugMessage("Warning: unknown beagle flag number: %d", bit);
#endif
				}
			else
				list += (*itToName).second + " ";
			}
		bit = bit << 1;
		}while(bit <= (1 << 18));
	}

void CalculationManager::InitializeBeagle(int nTips, int nClas, int nHolders, int nstates, int nchar, int nrates){
	assert(useBeagle);

	termOnBeagleError = true;
	req_flags = ParseBeagleFlagString(requiredBeagleFlags);
	//dealing with rescaling depending on other assigned flags gets annoying, so just always demand it
	req_flags = req_flags | BEAGLE_FLAG_SCALERS_LOG;
	pref_flags = ParseBeagleFlagString(preferredBeagleFlags, req_flags);

	outman.UserMessage("BEAGLE INITIALIZING ...");

 	string str;
 	InterpretBeagleResourceFlags(req_flags, str);
	outman.UserMessage("Required beagle flags: %s", str.c_str());
	str.clear();
	InterpretBeagleResourceFlags(pref_flags, str, req_flags);
	outman.UserMessage("Preferred beagle flags: %s", str.c_str());

	if(beagleDeviceNum > -1){
		outman.UserMessage("Requested beagle device number: %d", beagleDeviceNum);
		}

	outman.UserMessage("");
	OutputBeagleResources();

	outman.UserMessage("\nCREATING BEAGLE INSTANCE ...");

	//to allow ambiguity we need to create partials for the tips with ambiguity, and normal tip states for the others
	int tipCount = nTips;
	int ambigTips = data->NumTaxaWithPartialAmbig();
	int normalTips = nTips - ambigTips;
	int partialsCount = nClas + ambigTips;
	int compactCount = normalTips;

#ifndef DONT_SEND_TRANSMATS
	//num freqs = numEigens in beagle, so need at least one
	int eigCount = 1;
#else
	int eigCount = nClas;
#endif
	//this is a gross overestimate of how many matrices are needed, but haven't worked out recycling quite yet for mats
	int matrixCount = (nHolders * 3) * 2;//x3 for the pmats, d1mats and d2mats, x2 for both internals and terms
	//int matrixCount = (nClas * 3) * 2;//x3 for the pmats, d1mats and d2mats, x2 for both internals and terms
	
	//try one scaler per cla, as in normal garli.  These are doubles rather than ints though, so larger.
	//scaler for a given cla will share same index, and will be cumulative, as mine are now
	//add one for a destinationScaleWrite which will always be the last and will be used in all calls for scratch
	//Note that although there is a 1to1 correspondence between the cla and scaler indeces, this is NOT true for what is
	//passed to beagle, so beagle partial index = my partial index + #tips, beagle scaler = my scaler

	//DEBUG - it is hard to tell whether we will later be forced in a SP instance, in which case we'd like to turn on rescaling
	//so, just allocate enough scalers to do it whether we will or not
	//int scalerCount = (IsRescaling() ? nClas + 1 : 0);
	int scalerCount = nClas + 1;
	int resourceList[1] = {0};
	int resourceListCount = 0;

	if(beagleDeviceNum > -1){
		resourceList[0] = beagleDeviceNum;
		resourceListCount = 1;
		}

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
		pref_flags, 
		req_flags,
		&det);

	CheckBeagleReturnValue(
		beagleInst, 
		"beagleCreateInstance");

	ParseInstanceDetails(&det);

	int fpSize = (IsSinglePrecision() ? 4 : 8);

	outman.DebugMessage("BEAGLE ALLOCATIONS:");
	outman.DebugMessage("\tstates: %d char: %d rates: %d", nstates, nchar, nrates);
	outman.DebugMessage("\ttips: %d", tipCount);
	outman.DebugMessage("\tcompact tips arrays (no partial ambiguity): %d (%.1f KB)", compactCount, compactCount * nchar * 4 / 1024.0);
	outman.DebugMessage("\tpartial tips arrays (some partial ambiguity): %d (%.1f KB)", ambigTips, ambigTips * nstates * nchar * fpSize / 1024.0);
	outman.DebugMessage("\ttotal partial arrays: %d (%.1f KB)", partialsCount, partialsCount * nstates * nchar * nrates * fpSize / 1024.0);
	outman.DebugMessage("\teigen solutions: %d (%.1f KB (?))", eigCount, eigCount * nstates * 2 * fpSize / 1024.0);
	outman.DebugMessage("\ttransition matrices (pmats and derivs): %d (%.1f KB)", matrixCount, matrixCount * nstates * nstates * nrates * fpSize / 1024.0);
	outman.DebugMessage("\trescaling arrays %d (%.1f KB)", scalerCount, scalerCount * nchar * fpSize / 1024.0);
	
	double totMem = (compactCount * nchar * 4 / 1024.0) + (ambigTips * nstates * nchar * fpSize / 1024.0) + (partialsCount * nstates * nchar * nrates * fpSize / 1024.0) +
			(eigCount * nstates * 2 * fpSize / 1024.0) + (eigCount, eigCount * nstates * 2 * fpSize / 1024.0) + (matrixCount * nstates * nstates * nrates * fpSize / 1024.0) +
			(scalerCount * nchar * fpSize / 1024.0);

	outman.DebugMessage("Total beagle memory allocations approx %.1f MB)", totMem / 1024.0);

	SendTipDataToBeagle();

	vector<double> counts;
	for(int pat = 0;pat < data->NChar();pat++)
		counts.push_back((double) data->Count(pat));
#ifdef OUTPUT_OTHER_BEAGLE
	outman.DebugMessageNoCR("counts ");
	for(vector<double>::iterator it = counts.begin();it != counts.end();it++){
		outman.DebugMessageNoCR("%d ", (int) *it);
		}
	outman.DebugMessage("");
#endif
#ifdef OUTPUT_BEAGLE_SITELIKES
	string name = ofprefix + ".Bsitelikes.log";
	ofstream file(name.c_str());	
	file.close();
#endif
#ifdef OUTPUT_BEAGLE_PARTIALS
	string name2 = ofprefix + ".partials.log";
	ofstream part(name2.c_str());	
	part.close();
#endif
	beagleSetPatternWeights(beagleInst, &(counts[0]));

	outman.UserMessage("#######################################################");
	}

void CalculationManager::SendTipDataToBeagle(){
	outman.DebugMessage("SENDING DATA");
    char convert[16]={-1, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4};
	
	int nstates = data->NStates();
	for(int t = 0;t < data->NTax();t++){
		bool partialAmbig = data->TaxonHasPartialAmbig(t);
		outman.DebugMessageNoCR("t %d ", t);

		const unsigned char *dataString = data->GetRow(t); 
		if(partialAmbig){
			//ambiguity if currently for nuc only (all versions, not just beagle)
			assert(nstates == 4);
			outman.DebugMessage("(some ambig)");
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
			outman.DebugMessage("(no ambig)");
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

ScoreSet CalculationManager::CalculateLikelihoodAndDerivatives(const TreeNode *effectiveRoot, bool calcDerivs){
#ifdef DEBUG_OPS
	if(calcDerivs)
		outman.DebugMessage("#\nCALC D RT=%d", effectiveRoot->nodeNum);
	else
		outman.DebugMessage("#\nCALC L RT=%d", effectiveRoot->nodeNum);
#endif

	//this collects all of the necessary computation details and stores them in 
	//the operationSetQueue and scorOps list
	//NOTE that any clas that are indentified as clean dependencies will be temp reserved here
	DetermineRequiredOperations(effectiveRoot, calcDerivs);	

	//list<BlockingOperationsSet> operationSetQueue;
		//list claOps
			//dest and child clas, transmat ind, deplevel
		//list pmatOps
			//destTransMat, edgelen, model index, calcDerivs
    
	//list<ScoringOperation> scoreOps;
		//dest and children, transmatInd, calcDerivs

	//do all of the ops in the operationsSetQueue (i.e., all transmat ops and updatePartials, 
	//INCLUDING the transmat needed for the scoring operation)
	UpdateAllConditionals();

	//Now do the score op, which will return a ScoreSet that has at least the lnL, and 
	//D1 and D2 if requested.  This assumes that there is only one score op in scoreOps
	assert(scoreOps.size() == 1);
	ScoreSet values;
	for(list<ScoringOperation>::iterator sit = scoreOps.begin() ; sit != scoreOps.end() ; sit++)
		values = PerformScoringOperation(&(*sit));

	ResetDepLevelsAndReservations();

	operationSetQueue.clear();
	scoreOps.clear();
	return values;
	}	

/*This is called by CalculateLikelihoodAndDerivatives and amasses all of the transmat, cla and scoring
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
		//here we are calculating the lnL only, and don't care which branch it is integrated on.  pick an internal one
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
		assert(start->claOp.opDepLevel == lvl);
		if(start->claOp.opDepLevel == lvl){
			while(end != nodeOps.end() && end->claOp.opDepLevel == lvl){
				end++;
				num++;
				}
			//outman.UserMessage("level %d = %d nodes", lvl, num);
			for(list<NodeOperation>::iterator nit = start;nit != end;nit++){
				thisSet.claOps.push_back((*nit).claOp);
				thisSet.pmatOps.push_back((*nit).transOp1);
				thisSet.pmatOps.push_back((*nit).transOp2);
				}
			thisSet.opSetDepLevel = lvl;
			operationSetQueue.push_back(thisSet);
			start = end;
			thisSet.claOps.clear();
			thisSet.pmatOps.clear();
			lvl++;
			}
		}

	//the matrix for the branch between end1 and end2 for the final operation 
	//there is no cla opt in this final operation, the calculation will instead 
	//be stored in the scoreOps list
	TransMatOperation final(finalPmat, 0, pmatMan->GetEdgelen(finalPmat), derivatives);
	BlockingOperationsSet fin;
	fin.pmatOps.push_back(final);
	operationSetQueue.push_back(fin);
	scoreOps.push_back(ScoringOperation(-1, end1, end2, finalPmat, derivatives));
	}

//this will recursively move over the tree from the root node, spreading until tips or clean clas
//are hit.  If the endpoints are clas, they will be reserved now
int CalculationManager::AccumulateOpsOnPath(int holderInd, list<NodeOperation> &nodeOps){
	//it is ok to call this with a tip, because that tip could be a terminal branch that is having derivs
	//calculated.  accumulate nothing in that case
	if(holderInd < 0)
		return 0;

	CondLikeArrayHolder *holder = claMan->GetMutableHolder(holderInd);
	if(! holder->IsDirty()){
		//this cla is clean, and now an initial dependency, so be sure that it doesn't get yoniked
		claMan->ClaimClaFillIfNecessary(holderInd, 0);
		return 0;
		}

	int depLevel1 = 0, depLevel2 = 0;
	
	if(holder->HolderDep1() >= 0){
		depLevel1 = AccumulateOpsOnPath(holder->HolderDep1(), nodeOps);
		}
	if(holder->HolderDep2() >= 0){
		depLevel2 = AccumulateOpsOnPath(holder->HolderDep2(), nodeOps);
		}
	
	claMan->SetDepLevel(holderInd, max(depLevel1, depLevel2) + 1);

/*
	BlockingOperationsSet opSet;
	//opSet.claOps.push_back(ClaOperation(holderInd, holder->HolderDep1(), holder->HolderDep2(), holder->TransMatDep1(), holder->TransMatDep2(), holder->DepLevel()));
	opSet.claOps.push_back(ClaOperation(holderInd, holder));
	//DEBUG - second arg here is model index.  Not sure how that will be dealt with
	opSet.pmatOps.push_back(TransMatOperation(holder->TransMatDep1(), 0, pmatMan->GetEdgelen(holder->TransMatDep1()), false));
	opSet.pmatOps.push_back(TransMatOperation(holder->TransMatDep2(), 0, pmatMan->GetEdgelen(holder->TransMatDep2()), false));

	operationSetQueue.push_back(opSet);
*/
	nodeOps.push_back(NodeOperation(ClaOperation(holderInd, holder->HolderDep1(), holder->HolderDep2(), holder->TransMatDep1(), holder->TransMatDep2(), holder->DepLevel()),
						TransMatOperation(holder->TransMatDep1(), 0, pmatMan->GetEdgelen(holder->TransMatDep1()), false),
						TransMatOperation(holder->TransMatDep2(), 0, pmatMan->GetEdgelen(holder->TransMatDep2()), false)));

	return holder->DepLevel();
	}

//this just performs all of the transmat and cla ops that were queued
//will generally be followed by a call to PerformScoring to get lnL or derivs
void CalculationManager::UpdateAllConditionals(){	
	//DEBUG
#ifdef DEBUG_OPS
	OutputOperationsSummary();
#endif

#define BATCHED_CALLS
#ifdef BATCHED_CALLS
/*
	//this just yoinks the pmat ops from each BlockingOpSet into a new list
	//note that PerformTransMatOperationBatch assumes that all the ops either
	//need derivs or not.  So, avoid adding ops that need derivs here.
	list<TransMatOperation> tops;
	for(list<BlockingOperationsSet>::iterator it = operationSetQueue.begin();it != operationSetQueue.end();it++){
		if(! (*it).pmatOps.begin()->calcDerivs)
			tops.splice(tops.end(), (*it).pmatOps);
		}

	if(tops.size() > 0){
#ifndef DONT_SEND_TRANSMATS
		SendTransMatsToBeagle(tops);
#else
		PerformTransMatOperationBatch(tops);
#endif
		}

	//Any transmats that need derivs will still be included in a blockingOpSet, so will get done here.  Otherwise
	//no transmat calls will be done here because the ops were pulled out above to be done in one big batch
*/
	/*Calculation of all necessary partials will go like this:
		0. Before getting here all cla ops will have their cla and pmat dependencies set, and
		   the dependencies of the first destination clas will have been reserved if clean clas
		1. (calc any transmats necessary - currently this will actually happen above)
		2. Claim the first set of destination clas (this could in theory cause a call to Recycle, but 
		   shouldn't be necessary unless something is wrong, I think)
		3. Calc the first level
		4. Unreserve deps of first level (which would be clas clean to begin with), don't add to freeable
		5. Claim second (could call recyc, but wouldn't touch anything used in this entire scoring/deriv call)
		6. Calc second 
		7. Mark dependencies of second (i.e., the first level) as freeable and push to claman freeable queue
		6. Claim third (could call recyc, AND would first reclaim the first level from claman freeable queue)
		7. Calc third
		6. Mark dependencies of third (i.e., the second level) as freeable and push to claman freeable queue
		7. Claim fourth (could call recyc, AND would first reclaim the first and second level from claman freeable queue)
		6. Calc fourth
		8. Mark deps of fourth as freeable and push to claman freeable queue
		etc...
		9. Unreserve anything remaining, clear freeable queue and op queue

	*/

	//BMERGE DEBUG
//#define DONT_SEND_TRANSMATS

	claMan->SetCurReclaimableLevel(0);
	list<BlockingOperationsSet>::iterator toFree = operationSetQueue.begin();

	for(list<BlockingOperationsSet>::iterator it = operationSetQueue.begin();it != operationSetQueue.end();){
		if((*it).pmatOps.size() > 0){
#ifndef DONT_SEND_TRANSMATS
			SendTransMatsToBeagle((*it).pmatOps);
#else
			PerformTransMatOperationBatch((*it).pmatOps);
#endif
			}
		if((*it).claOps.size() > 0){
			//this will reserve the DESTINATION clas. The child cla dependencies should already have
			//been reserved if they are clas (not tips).  If they were clean in the first place they
			//should have been reserved in AccumulateOps, otherwise when they were reserved when previous
			//destinations
			ReserveDestinationClas(*it);
			PerformClaOperationBatch((*it).claOps);
			
			//No need to do all this unreservation if the number of clas is large enough to simultaneously hold
			//all clas of the tree. We'll still unreserve them all at the end
			//DEBUG
			//if(memLevel > 0){
			if(1){
				//remove the temp reservations from the level 1 DEPENDENCIES (if they are clean clas and not tips)
				//but DO NOT add them to the freeable queue (since they must have been clean in the first place)
				if((*it).opSetDepLevel == 1){
					UnreserveDependencyClas(*it);
					claMan->SetCurReclaimableLevel((*it).opSetDepLevel - 1);
					it = operationSetQueue.erase(it);
					}
				//when we've finished level 2 or greater, we can now tell the claMan that the dependencies are ripe for the taking
				else if((*it).opSetDepLevel > 1){
					QueueDependencyClasForFreeing(*it);
					claMan->SetCurReclaimableLevel((*it).opSetDepLevel - 1);
					it = operationSetQueue.erase(it);
					}
				}
			else it++;
			}
		else 
			it++;
		}
	//TODO - what should happen with reclaim level after all the calcs?  
	//trying to allow freeing of first level dependencies in general
	claMan->SetCurReclaimableLevel(1);
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

void CalculationManager::PerformTransMatOperationBatch(const list<TransMatOperation> &theOps){
#ifndef USE_BEAGLE
	//mod->AltCalcPmat(theOp->edgeLength, pmatMan->GetPmat(theOp->destTransMatIndex)->theMat);
	//this is a bit ridiculous, but these funcs won't really be getting used except with beagle
	pmatMan->GetMutableHolder(theOp->destTransMatIndex)->myMod->AltCalcPmat(theOp->edgeLength, pmatMan->GetPmat(theOp->destTransMatIndex)->theMat);
#else 

	//for now assuming that all transmat ops in a set have same eigen solution and rate multipliers
	
	//currently just using a single eigen index and sending it every time
	int eigenIndex = 0;
	
	//this call will calculate the eigen solution first if necessary
	//currently it will always be resent (sol.changed will always be true) because
	//it is non-trivial to figure out when it has changed since it was last sent to beagle
	ModelEigenSolution sol;
	pmatMan->GetEigenSolution(theOps.begin()->destTransMatIndex, sol);

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
	pmatMan->GetCategoryRatesForBeagle(theOps.begin()->destTransMatIndex, categRates);

	//as above, always resend currently
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
		//DEBUG
		//outman.DebugMessage("transmat %d mod %d modadr %d edge %.4f", (*it).destTransMatIndex, (*it).modelIndex, pmatMan->GetCorrespondingModel((*it).destTransMatIndex), (*it).edgeLength);
		//pmatMan->GetCorrespondingModel(theOps.begin()->destTransMatIndex)->NumRateCatsForBeagle();
		pmatMan->ClaimTransmatsetFillIfNecessary((*it).destTransMatIndex);
		pmatInd.push_back(PmatIndexForBeagle((*it).destTransMatIndex));
		d1MatInd.push_back(D1MatIndexForBeagle((*it).destTransMatIndex));
		d2MatInd.push_back(D2MatIndexForBeagle((*it).destTransMatIndex));
		edgeLens.push_back(pmatMan->GetEdgelen((*it).destTransMatIndex));
		}

	int count = pmatInd.size();
	bool calcDerivs = theOps.begin()->calcDerivs;

#ifdef OUTPUT_OTHER_BEAGLE
	outman.DebugMessageNoCR("UPDATING TRANS MAT ");
//	outman.DebugMessage("%d (%d), eigen %d, blen %f", PmatIndexForBeagle(theOp->destTransMatIndex), theOp->destTransMatIndex, eigenIndex, edgeLens[0]);
#endif
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
	
	int nrates = pmatMan->GetModelForTransMatHolder(theOps.begin()->destTransMatIndex)->NumRateCatsForBeagle();

//	if(calcDerivs){

	/*wrong	*/
	//double inmat[16] = {0.9070614,	0.02927732,	0.02870608,	0.0349794, 0.03213807,	0.90418679,	0.02868981,	0.03497981, 0.03213807,	0.02926432,	0.90361226, 0.03497981, 0.03213831,	0.02927236,	0.02868767,	0.90990174};
	//double inmat[16] = {0.9070613980, 0.0292773210, 0.0287060775, 0.0349793993, 0.0321380682, 0.9041867852, 0.0286898147, 0.0349798091, 0.0321380682, 0.0292643178, 0.9036122561, 0.0349798091, 0.0321383066, 0.0292723589, 0.0286876746, 0.9099017382}; 
	//wrong, normalized rows sum to zero
	//double inmat[16] = {0.9070394515, 0.0292766126, 0.0287053830, 0.0349785530, 0.0321382457, 0.9041917789, 0.0286899731, 0.0349800023, 0.0321382465, 0.0292644802, 0.9036172701, 0.0349800032, 0.0321383041, 0.0292723566, 0.0286876724, 0.9099016670};
	//wrong, rescaled such that column totals are same, three lower elements in each column match
	//double inmat[16] = {0.9070613980, 0.0292713326, 0.0286945223, 0.0349796725, 0.0321381477, 0.9041867852, 0.0286945223, 0.0349796725, 0.0321381477, 0.0292713326, 0.9036122561, 0.0349796725, 0.0321381477, 0.0292713326, 0.0286945223, 0.9099017382};
	/*right*/
	//double inmat[16] = {0.9070611, 0.0292719, 0.02868779, 0.03497917, 0.03213868, 0.90419436, 0.02868777, 0.03497917, 0.03213868, 0.0292719, 0.90361023, 0.03497917, 0.03213868, 0.02927192, 0.02868777, 0.90990162};
	//double inmat[16] = {0.9070611396, 0.0292718949, 0.0286877828, 0.0349791827, 0.0321386719, 0.9041943625, 0.0286877828, 0.0349791827, 0.0321386719, 0.0292718949, 0.9036102504, 0.0349791827, 0.0321386719, 0.0292718949, 0.0286877828, 0.9099016503};
	//pmat from DP beagle
	double inmat[16] = {0.907060771,  0.02927234,  0.028687678,  0.034979238,  0.032138303,  0.904194759,  0.028687707,  0.034979216,  0.032138303,  0.029272277,  0.903610189,  0.034979216,  0.032138297,  0.029272304,  0.028687689,  0.909901709};

	vector<double> outMat(nstates * nstates * nrates);
	for(list<TransMatOperation>::const_iterator pit = theOps.begin();pit != theOps.end();pit++){				
//		outman.UserMessage("MANUALLY SETTING PMATS");
//		beagleSetTransitionMatrix(beagleInst, PmatIndexForBeagle((*pit).destTransMatIndex), inmat);

		outMat.clear();
		for(int p=0;p<nstates*nstates*nrates;p++) 
			outMat.push_back(p);
		beagleGetTransitionMatrix(beagleInst, PmatIndexForBeagle((*pit).destTransMatIndex), &(outMat[0]));
		//beagleGetTransitionMatrix(beagleInst, D1MatIndexForBeagle((*pit).destTransMatIndex), &(outMat[0]));
		vector<double>::iterator it = outMat.begin();
		outman.DebugMessage("%d", (*pit).destTransMatIndex);
		for(int r = 0;r < nrates;r++){
			for(int i = 0;i < nstates;i++){
				for(int j=0;j < nstates;j++){
					outman.DebugMessageNoCR("%.10f\t", *it++);
					}
				outman.DebugMessage("");
				}
			}
		}
//		}
#endif
	}

void CalculationManager::SendTransMatsToBeagle(const list<TransMatOperation> &theOps){
#if ! (defined(BATCHED_CALLS) && ! defined(DONT_SEND_TRANSMATS))
	assert(0);
#endif

#ifdef OUTPUT_PMATS
	ofstream *deb = outman.GetDebugStream();
#endif

	for(list<TransMatOperation>::const_iterator pit = theOps.begin();pit != theOps.end();pit++){

#define NEW
#ifdef NEW
		int theIndex = (*pit).destTransMatIndex;
		FLOAT_TYPE ***thePMat = NULL;
		FLOAT_TYPE ***theD1Mat = NULL;
		FLOAT_TYPE ***theD2Mat = NULL;
		pmatMan->ObtainAppropriateTransmats(*pit, thePMat, theD1Mat, theD2Mat);
#else
		double theBlen = (*pit).edgeLength;
		int theIndex = (*pit).destTransMatIndex;
		FLOAT_TYPE ***thePMat = pmatMan->GetPmatArray(theIndex);
		FLOAT_TYPE ***theD1Mat = NULL;
		FLOAT_TYPE ***theD2Mat = NULL;
		if(pit->calcDerivs){
			theD1Mat = pmatMan->GetD1MatArray(theIndex);
			theD2Mat = pmatMan->GetD2MatArray(theIndex);
			pmatMan->GetMutableModelForTransMatHolder(theIndex)->FillDerivativeMatrices(theBlen, thePMat, theD1Mat, theD2Mat);
			pmatMan->FillDerivativeMatrices(theIndex, thePMat, theD1Mat, theD2Mat);
			}
		else{
			pmatMan->GetMutableModelForTransMatHolder(theIndex)->AltCalcPmat(theBlen, thePMat);
			}
#endif
#ifdef OUTPUT_PMATS
		//output pmats to be sent
		*deb << "send p " << theIndex << "\n";
		pmatMan->GetMutableModelForTransMatHolder(theIndex)->OutputPmat(*deb, thePMat);
#endif
		beagleSetTransitionMatrix(beagleInst, PmatIndexForBeagle(theIndex), **thePMat, 1.0);

#ifdef OUTPUT_PMATS
		/*BMERGE this doesn't seem to have been working
		if(!gpuBeagle){
			*deb << "ret p GI=" << theIndex << " BI=" << PmatIndexForBeagle(theIndex)<< "\n";
			OutputBeagleTransMat(PmatIndexForBeagle(theIndex));
			}
		*/
#endif
		if(pit->calcDerivs){
#ifdef OUTPUT_PMATS
			*deb << "send D1 " << theIndex << "\n";
			pmatMan->GetMutableModelForTransMatHolder(theIndex)->OutputPmat(*deb, theD1Mat);
			*deb << "send D2 " << theIndex << "\n";
			pmatMan->GetMutableModelForTransMatHolder(theIndex)->OutputPmat(*deb, theD2Mat);
#endif
			beagleSetTransitionMatrix(beagleInst, D1MatIndexForBeagle(theIndex), **theD1Mat, 0.0);
			beagleSetTransitionMatrix(beagleInst, D2MatIndexForBeagle(theIndex), **theD2Mat, 0.0);

#ifdef OUTPUT_PMATS
			/*BMERGE this doesn't seem to be working
			if(!gpuBeagle){
				if((*pit).calcDerivs){
					*deb << "ret D1 GI=" << theIndex << " BI=" << D1MatIndexForBeagle(theIndex) << "\n";
					OutputBeagleTransMat(D1MatIndexForBeagle(theIndex));
					*deb << "ret D2 GI=" << theIndex << " BI=" << D2MatIndexForBeagle(theIndex) << "\n";
					OutputBeagleTransMat(D2MatIndexForBeagle(theIndex));
					}
				}
			*/
				
#endif
			}
		}
	}

void CalculationManager::PerformClaOperationBatch(const list<ClaOperation> &theOps){
	//not sure if this is right - will always use a single scale array for destWrite (essentially
	//scratch space, I think) and then pass a cumulative scaler to actually keep track of the scaling
	//int destinationScaleWrite = (rescaleBeagle ? claMan->NumClas() : BEAGLE_OP_NONE);

	//scale read is if you are telling it to use precomputed scaling factors
	int	destinationScaleRead = BEAGLE_OP_NONE;

	vector<int> tupleList;
	for(list<ClaOperation>::const_iterator it = theOps.begin();it != theOps.end();it++){
		//This is happening a level higher now
		//this will fill the destination holder with a cla index if necessary
		//claMan->ClaimClaFillIfNecessary((*it).destClaIndex, (*it).depLevel);


		tupleList.push_back(PartialIndexForBeagle((*it).destClaIndex));
		//if rescaling, for each partial op pass in an array to store the amount of rescaling
		//at this node, in my case stored as logs
		if(IsRescaling())
			tupleList.push_back(ScalerIndexForBeagle((*it).destClaIndex));
		else 
			tupleList.push_back(BEAGLE_OP_NONE);
        tupleList.push_back(destinationScaleRead);
		tupleList.push_back(PartialIndexForBeagle((*it).childClaIndex1));
		tupleList.push_back(PmatIndexForBeagle((*it).transMatIndex1));
		tupleList.push_back(PartialIndexForBeagle((*it).childClaIndex2));
		tupleList.push_back(PmatIndexForBeagle((*it).transMatIndex2));
		}

	int instanceCount = 1;
	int operationCount = (int) theOps.size();

	//accumulate rescaling factors - For scale arrays my indexing scheme and Beagle's happen to be the same, 
	//and my negative (tip) corresponds to a NULL in beagle
	int cumulativeScaleIndex = BEAGLE_OP_NONE;

	CheckBeagleReturnValue(
		beagleUpdatePartials(
			beagleInst,
			(BeagleOperation*)&tupleList[0],
			operationCount,
			cumulativeScaleIndex),
		"beagleUpdatePartials");

#ifdef OUTPUT_PARTIALS
	string name2 = ofprefix + ".partials.log";
	ofstream part(name2.c_str());
	part.precision(5);

	vector<double> outPartials(data->NStates() * data->NChar() * pmatMan->GetNumRates());
	for(list<ClaOperation>::const_iterator it = theOps.begin();it != theOps.end();it++){
		//could pass a scaler index here to have the values unscaled, I think
		beagleGetPartials(beagleInst, PartialIndexForBeagle((*it).destClaIndex), BEAGLE_OP_NONE, &outPartials[0]);
		for(int site = 0;site < data->NChar();site++){
			part << site << "\t";
			for(int rate = 0;rate < pmatMan->GetNumRates();rate++){
				for(int state = 0;state < 4;state++)
					part << outPartials[site * 4 + state] << "\t";
				}
			part << "\n";
			}
		}
	part.close();
#endif

	//this is a bit sneaky. each partial now has a rescaling array associated with it, holding
	//the amount done at that node alone.  Take the rescaling amounts from the two children,
	//and add then to the rescaling already done at this node.  In the next dependency level pass
	//the dest here will become the child there.

	if(IsRescaling()){
		for(list<ClaOperation>::const_iterator it = theOps.begin();it != theOps.end();it++){
			AccumulateRescalers(
				(*it).destClaIndex, 
				(*it).childClaIndex1, 
				(*it).childClaIndex2);
			}
		}
/*
	//to indicate that the partial is now clean, fill the holder that represents it, despite the fact that the memory there is never being used
	//DEBUG - need to figure out second argument here (direction) which I think determined how likely a holder is to be recycled.
	for(list<ClaOperation>::const_iterator it = theOps.begin();it != theOps.end();it++){
		claMan->FillHolder((*it).destClaIndex, 0);
		}
*/
	}

#ifdef USE_BEAGLE
//For scale arrays my indexing scheme and Beagle's happen to be the same
//further accumulate the cumulative rescalings from the two children.  In the batched case
//the dest array will already have the rescaling done at this node included, so don't reset it
void CalculationManager::AccumulateRescalers(int destIndex, int childIndex1, int childIndex2){

#ifdef OUTPUT_OTHER_BEAGLE
	//DEBUG
	//outman.DebugMessage("accum %d + %d = %d", childIndex1, childIndex2, destIndex);
#endif
#ifndef BATCHED_CALLS
	//clear out the destination scaler first
	CheckBeagleReturnValue(
		beagleResetScaleFactors(beagleInst,
			ScalerIndexForBeagle(destIndex)), 
		"beagleResetScaleFactors");
#endif

	vector<int> childScalers;
	if( ! (childIndex1 < 0))
		childScalers.push_back(ScalerIndexForBeagle(childIndex1));
	if( ! (childIndex2 < 0))
		childScalers.push_back(ScalerIndexForBeagle(childIndex2));

	if(childScalers.size() > 0){
		CheckBeagleReturnValue(
			beagleAccumulateScaleFactors(beagleInst,
				&(childScalers[0]),
				(int) childScalers.size(),
				ScalerIndexForBeagle(destIndex)),
			"beagleAccumulateScaleFactors");
		}
	}
#endif

void CalculationManager::OutputOperationsSummary() const{
	//DEBUG
	bool outputFullSummary = true;
	//bool outputFullSummary = false;

	int numPmats = 0, numClas = 0, numSets = 0;
	for(list<BlockingOperationsSet>::const_iterator it = operationSetQueue.begin();it != operationSetQueue.end();it++){
		numSets++;
		for(list<TransMatOperation>::const_iterator pit = (*it).pmatOps.begin();pit != (*it).pmatOps.end();pit++)
			numPmats++;
		for(list<ClaOperation>::const_iterator cit = (*it).claOps.begin();cit != (*it).claOps.end();cit++){
			numClas++;
			if(outputFullSummary){
				//DEBUG
				int D=(*cit).destClaIndex;
				int C1=(*cit).childClaIndex1;
				int C2=(*cit).childClaIndex2;
				outman.DebugMessage("D\t%d\t%d\tC\t%d\t%f\t%d\t%f M1 %d M2 %d", 
					(*cit).destClaIndex, 
					(*cit).opDepLevel, 
					(*cit).childClaIndex1, 
					pmatMan->GetEdgelen((*cit).transMatIndex1), 
					(*cit).childClaIndex2, 
					pmatMan->GetEdgelen((*cit).transMatIndex2), 
					(*cit).transMatIndex1,
					(*cit).transMatIndex2);
				}
			}
		}
	outman.DebugMessage("%d sets\t%d mats\t%d clas", numSets, numPmats, numClas);
	}

void CalculationManager::ReserveDestinationClas(BlockingOperationsSet &set){
	//outman.DebugMessage("reserving dests of level %d", set.opSetDepLevel);
	for(list<ClaOperation>::const_iterator op = set.claOps.begin(); op != set.claOps.end() ; op++) 
		claMan->ClaimClaFillIfNecessary((*op).destClaIndex, (*op).opDepLevel);
	}

void CalculationManager::UnreserveDependencyClas(const BlockingOperationsSet &set){
	//outman.DebugMessage("unreserving deps of level %d", set.opSetDepLevel);
	for(list<ClaOperation>::const_iterator op = set.claOps.begin(); op != set.claOps.end() ; op++){
		//This ONLY removes the temp reservations from the dependencies
		//It does NOT add them to the freeableQueue
		if((*op).childClaIndex1 > -1)
			claMan->RemoveTempReservation((*op).childClaIndex1);
		if((*op).childClaIndex2 > -1)
			claMan->RemoveTempReservation((*op).childClaIndex2);
		}
	}

void CalculationManager::QueueDependencyClasForFreeing(const BlockingOperationsSet &set){
	//outman.DebugMessage("queuing deps level %d for freeing", set.opSetDepLevel);
	for(list<ClaOperation>::const_iterator op = set.claOps.begin(); op != set.claOps.end() ; op++){
		//Adds DEPENDENCY clas to the freeable queue AND removes their temp reservations
		claMan->AddToFreeableQueue((*op).childClaIndex1);
		claMan->AddToFreeableQueue((*op).childClaIndex2);
		}
	}

void CalculationManager::QueueDestinationClasForFreeing(const BlockingOperationsSet &set){
	//outman.DebugMessage("queuing dests level %d for freeing", set.opSetDepLevel);
	int dontFree = scoreOps.begin()->childClaIndex1;
	int dontFree2 = scoreOps.begin()->childClaIndex2;
	for(list<ClaOperation>::const_iterator op = set.claOps.begin(); op != set.claOps.end() ; op++){
		int thisIndex = (*op).destClaIndex;
		//Adds destination clas to the freeable queue AND removes their temp reservations
		//this is annoying, but can't do with dests needed for final scoring.  They will
		//eventually get freed in ResetDepLevelsAndRes
		if(thisIndex != dontFree && thisIndex != dontFree2)
			claMan->AddToFreeableQueue(thisIndex);
		}
	}

//this version resets a single opset
void CalculationManager::ResetDepLevels(const BlockingOperationsSet &set){
	assert(0); //trying not resetting deps
	outman.DebugMessage("resetting deps of level %d", set.opSetDepLevel);
	for(list<ClaOperation>::const_iterator op = set.claOps.begin(); op != set.claOps.end() ; op++){
		claMan->SetDepLevel((*op).destClaIndex, 0);
		//unreseving the deps may be very redundant since they would be unreserved when they
		//are the destClaIndex in the next step.  However, the tip-most dependenices 
		//which might be clean clas cause reservations as well
		if((*op).childClaIndex1 > -1)
			claMan->SetDepLevel((*op).childClaIndex1, 0);
		if((*op).childClaIndex2 > -1) 
			claMan->SetDepLevel((*op).childClaIndex2, 0);
		}
	}

//this version resets all ops that remain in the queue
void CalculationManager::ResetDepLevelsAndReservations(){
	for(list<BlockingOperationsSet>::const_iterator set = operationSetQueue.begin(); set != operationSetQueue.end();set++){
		if((*set).claOps.size() > 0)
			outman.DebugMessage("resetting all - level %d", (*set).opSetDepLevel);
		for(list<ClaOperation>::const_iterator op = (*set).claOps.begin(); op != (*set).claOps.end() ; op++){
//			claMan->SetDepLevel((*op).destClaIndex, 0);
			claMan->RemoveTempReservation((*op).destClaIndex);
			if((*op).childClaIndex1 > -1){
//				claMan->SetDepLevel((*op).childClaIndex1, 0);
				claMan->RemoveTempReservation((*op).childClaIndex1);
				}
			if((*op).childClaIndex2 > -1){
//				claMan->SetDepLevel((*op).childClaIndex2, 0);
				claMan->RemoveTempReservation((*op).childClaIndex2);
				}
			}
		}
	for(list<ScoringOperation>::const_iterator sop = scoreOps.begin();sop != scoreOps.end();sop++){
		if((*sop).childClaIndex1 > -1){
//			claMan->SetDepLevel((*sop).childClaIndex1, 0);
			claMan->RemoveTempReservation((*sop).childClaIndex1);
			}
		if((*sop).childClaIndex2 > -1){
//			claMan->SetDepLevel((*sop).childClaIndex2, 0);
			claMan->RemoveTempReservation((*sop).childClaIndex2);
			}
		}
	claMan->ClearFreeableQueue();
	if(claMan->debug_clas)
		claMan->ReportClaTotals("reset all ", -1);
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
		if(theOp->childClaIndex1 < 0)
			term1 = -theOp->childClaIndex1;
	
		if(theOp->childClaIndex2 < 0)
			term2 = -theOp->childClaIndex2;

		if( term1 && term2 ){
			CalcFullCLATerminalTerminal(
				claMan->GetClaFillIfNecessary(theOp->destClaIndex),
				//DEBUG - should probably figure some way not to do these static cassts
				(char*) static_cast<const NucleotideData *>(data)->GetAmbigString(term1-1),
				(char*) static_cast<const NucleotideData *>(data)->GetAmbigString(term2-1),
				pmatMan->GetPmatArray(theOp->transMatIndex1),
				pmatMan->GetPmatArray(theOp->transMatIndex2), 
				pmatMan->GetModelForTransMatHolder(theOp->transMatIndex1));
			}
		else if( ! (term1 || term2)){
			CalcFullClaInternalInternal(
				claMan->GetClaFillIfNecessary(theOp->destClaIndex),
				claMan->GetClaFillIfNecessary(theOp->childClaIndex1),
				claMan->GetClaFillIfNecessary(theOp->childClaIndex2),
				pmatMan->GetPmatArray(theOp->transMatIndex1),
				pmatMan->GetPmatArray(theOp->transMatIndex2),
				pmatMan->GetModelForTransMatHolder(theOp->transMatIndex1));
			}
		else if( term1 ){
			this->CalcFullCLAInternalTerminal(
				claMan->GetClaFillIfNecessary(theOp->destClaIndex),
				claMan->GetClaFillIfNecessary(theOp->childClaIndex2),
				(char*) static_cast<const NucleotideData *>(data)->GetAmbigString(term1-1),
				pmatMan->GetPmatArray(theOp->transMatIndex2),
				pmatMan->GetPmatArray(theOp->transMatIndex1),
				pmatMan->GetModelForTransMatHolder(theOp->transMatIndex1), 
				NULL);
			}
		else{
			this->CalcFullCLAInternalTerminal(
				claMan->GetClaFillIfNecessary(theOp->destClaIndex),
				claMan->GetClaFillIfNecessary(theOp->childClaIndex1),
				(char*) static_cast<const NucleotideData *>(data)->GetAmbigString(term2-1),
				pmatMan->GetPmatArray(theOp->transMatIndex1),
				pmatMan->GetPmatArray(theOp->transMatIndex2),
				pmatMan->GetModelForTransMatHolder(theOp->transMatIndex1),
				NULL);
			}
		}
	else{
#ifdef USE_BEAGLE

		//this is deprecated in favor of batched calls
#ifdef BATCHED_CALLS
	assert(0);
#endif	
		//not sure if this is right - will always use a single scale array for destWrite (essentially
		//scratch space, I think) and then pass a cumulative scaler to actually keep track of the scaling
		int destinationScaleWrite = (IsRescaling() ? claMan->NumClas() : BEAGLE_OP_NONE);
		int	destinationScaleRead = BEAGLE_OP_NONE;

		int operationTuple[7] = {PartialIndexForBeagle(theOp->destClaIndex),
                                destinationScaleWrite,
                                destinationScaleRead,
								PartialIndexForBeagle(theOp->childClaIndex1),
 								PmatIndexForBeagle(theOp->transMatIndex1),
 								PartialIndexForBeagle(theOp->childClaIndex2),
								PmatIndexForBeagle(theOp->transMatIndex2)};

		outman.DebugMessageNoCR("PARTS\t");
		outman.DebugMessageNoCR("\tD\t%d (%d)", PartialIndexForBeagle(theOp->destClaIndex), theOp->destClaIndex);
		outman.DebugMessageNoCR("\tC\t%d (%d)\t%d (%d)", PartialIndexForBeagle(theOp->childClaIndex1), theOp->childClaIndex1, PartialIndexForBeagle(theOp->childClaIndex2), theOp->childClaIndex2);
		outman.DebugMessage("\tP\t%d (%d)\t%d (%d)", PmatIndexForBeagle(theOp->transMatIndex1), theOp->transMatIndex1, PmatIndexForBeagle(theOp->transMatIndex2), theOp->transMatIndex2);
		
		int instanceCount = 1;
		int operationCount = 1;

		//accumulate rescaling factors - For scale arrays my indexing scheme and Beagle's happen to be the same, 
		//and my negative (tip) corresponds to a NULL in beagle
		int cumulativeScaleIndex = BEAGLE_OP_NONE;
		if(IsRescaling()){
			cumulativeScaleIndex = theOp->destClaIndex;
			AccumulateRescalers(cumulativeScaleIndex, theOp->childClaIndex1, theOp->childClaIndex2);
			}

		CheckBeagleReturnValue(
			beagleUpdatePartials(
				beagleInst,
				(BeagleOperation*)operationTuple,
				operationCount,
				cumulativeScaleIndex),
			"beagleUpdatePartials");

		//to indicate that the partial is now clean, fill the holder that represents it, despite the fact that the memory there is never being used
		//DEBUG - need to figure out second argument here (direction) which I think determined how likely a holder is to be recycled.
		claMan->FillHolder(theOp->destClaIndex, 0);
#endif
		}
	}



//deprecated in favor of PerformTransMatOperationBatch
void CalculationManager::PerformTransMatOperation(const TransMatOperation *theOp){
#ifndef USE_BEAGLE
	//mod->AltCalcPmat(theOp->edgeLength, pmatMan->GetPmat(theOp->destTransMatIndex)->theMat);
	//this is a bit ridiculous, but these funcs won't really be getting used except with beagle
	pmatMan->GetMutableHolder(theOp->destTransMatIndex)->myMod->AltCalcPmat(theOp->edgeLength, pmatMan->GetPmat(theOp->destTransMatIndex)->theMat);
#else 

#ifdef BATCHED_CALLS
	assert(0);
#endif

#ifdef OUTPUT_OTHER_BEAGLE
	outman.DebugMessage("SETTING EIGEN");
#endif
	//DEBUG - currently just using a single eigen index and sending it every time
	int eigenIndex = 0;
	
	//this call will calculate the eigen solution first if necessary
	ModelEigenSolution sol;
	pmatMan->GetEigenSolution(theOp->destTransMatIndex, sol);

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
	pmatMan->GetCategoryRatesForBeagle(theOp->destTransMatIndex, categRates);

#ifdef OUTPUT_OTHER_BEAGLE
	outman.DebugMessage("SETTING CATEGORY RATES");
#endif
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
	double edgeLens[1] = {pmatMan->GetEdgelen(theOp->destTransMatIndex)};
	int count = 1;
	
#ifdef OUTPUT_OTHER_BEAGLE
	outman.DebugMessageNoCR("UPDATING MAT ");
	outman.DebugMessage("%d (%d), eigen %d, blen %f", PmatIndexForBeagle(theOp->destTransMatIndex), theOp->destTransMatIndex, eigenIndex, edgeLens[0]);
#endif

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
	int nrates = pmatMan->GetModelForTransMatHolder(theOp->destTransMatIndex)->NumRateCatsForBeagle();
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
void CalculationManager::OutputBeagleTransMat(int beagleIndex){
	int nstates = data->NStates();
	int nrates = pmatMan->NumRateCatsForBeagle(beagleIndex);
	vector<double> outMat(nstates * nstates * nrates);

	beagleGetTransitionMatrix(beagleInst, beagleIndex, &(outMat[0]));
	vector<double>::iterator it = outMat.begin();
	for(int r = 0;r < nrates;r++){
		outman.DebugMessage("r%d", r);
		for(int i = 0;i < nstates;i++){
			for(int j=0;j < nstates;j++){
				outman.DebugMessageNoCR("%.6f\t", *it++);
				}
			outman.DebugMessage("");
			}
		}
	}

//this does the last operation, beagleCalculateEdgeLogLikelihoods, and expects that the 
//required partials have been calculated as well as the transmats necessary for this operation
ScoreSet CalculationManager::PerformScoringOperation(const ScoringOperation *theOp){
	ScoreSet results = {0.0, 0.0, 0.0};
#ifdef OUTPUT_OTHER_BEAGLE
//	outman.DebugMessage("ENTER PERF SCR");
#endif

	if(!useBeagle){
		results.lnL = GetScorePartialInternalRateHet(claMan->GetClaFillIfNecessary(theOp->childClaIndex1), 
			claMan->GetClaFillIfNecessary(theOp->childClaIndex2), 
			pmatMan->GetPmatArray(theOp->transMatIndex),
			pmatMan->GetModelForTransMatHolder(theOp->transMatIndex));
		}

	else{
#ifdef USE_BEAGLE
		//state freqs - these are always being sent and stored in the same slot
		vector<FLOAT_TYPE> freqs(pmatMan->GetNumStates());
		//BEAGLE MERGE
		pmatMan->GetModelForTransMatHolder(theOp->transMatIndex)->GetStateFreqs(&(freqs[0]));
		//pmatMan->GetModelForTransMatHolder(theOp->transMatIndex)->stateFreqs;

#ifdef OUTPUT_OTHER_BEAGLE
/*		outman.DebugMessageNoCR("freqs sent\t");
		for(int i = 0; i < 4;i++)
			outman.DebugMessageNoCR("%.6f\t", freqs[i]);
		outman.DebugMessage("");
		*/
#endif
		int freqIndex = 0;
		CheckBeagleReturnValue(
			beagleSetStateFrequencies(
				beagleInst,
				freqIndex,
				&(freqs[0])),
			"beagleSetStateFrequencies");

		//category weights (or probs) - these are also always being sent.  annoying to figure out when they need to be updated
		vector<FLOAT_TYPE> inWeights;
		pmatMan->GetCategoryWeightsForBeagle(theOp->transMatIndex, inWeights);
#ifdef OUTPUT_OTHER_BEAGLE
		if(inWeights.size() > 1){
			outman.DebugMessageNoCR("rate props sent\t");
			for(int i = 0; i < inWeights.size();i++)
				outman.DebugMessageNoCR("%.6f\t", inWeights[i]);
			outman.DebugMessage("");
			}
#endif
		int weightIndex = 0;
		if(inWeights.size() > 0){
			CheckBeagleReturnValue(
				beagleSetCategoryWeights(
					beagleInst,
					weightIndex,
					&(inWeights[0])),
				"beagleSetCategoryWeights");
			}

		//the two children to be included- they should have been calculated and reserved
		//earlier in UpdateClas
		int buffer1[1] = {PartialIndexForBeagle(theOp->childClaIndex1)};
		int buffer2[1] = {PartialIndexForBeagle(theOp->childClaIndex2)};
		int numInput = 2;

		//arrays to hold site lnLs and derivs returned from Beagle
		//note that later beagle versions don't return the array anymore, just the sum
		vector<double> siteLikesOut(data->NChar());
		vector<double> siteD1Out(data->NChar());
		vector<double> siteD2Out(data->NChar());

		int pmatIndeces[1] = {PmatIndexForBeagle(theOp->transMatIndex)};
		int	d1MatIndeces[1] = {D1MatIndexForBeagle(theOp->transMatIndex)};
		int	d2MatIndeces[1] = {D2MatIndexForBeagle(theOp->transMatIndex)};

#ifdef OUTPUT_OTHER_BEAGLE
		outman.DebugMessageNoCR("Scr:\t");
		outman.DebugMessageNoCR("\tC\t%d (%d)\t%d (%d)", buffer2[0], theOp->childClaIndex2, buffer1[0], theOp->childClaIndex1);
		outman.DebugMessage("\tP\t%d (%d)", pmatIndeces[0], theOp->transMatIndex);
#endif

		int count = 1;

		//accumulate rescaling factors of the two clas that are being combined (one might be a tip)
		//For scale arrays my indexing scheme and Beagle's happen to be the same, and my negative (tip) corresponds to a NULL in beagle
		//use this scratch index to hold the final accumulated scalers 
		int cumulativeBeagleScaleIndex = BEAGLE_OP_NONE;
		if(IsRescaling()){
			cumulativeBeagleScaleIndex = ScalerIndexForBeagle(GARLI_FINAL_SCALER_INDEX);
			beagleResetScaleFactors(beagleInst, cumulativeBeagleScaleIndex);
			//add all of the rescaling up to the two involved children.  This will be passed to
			//calcEdgeLike, which will give the final scores
			AccumulateRescalers(GARLI_FINAL_SCALER_INDEX, theOp->childClaIndex1, theOp->childClaIndex2);
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
				&cumulativeBeagleScaleIndex,
				count,
				&siteLikesOut[0],
				(theOp->derivatives ? &siteD1Out[0] : NULL),
				(theOp->derivatives ? &siteD2Out[0] : NULL)),
			"beagleCalculateEdgeLogLikelihoods");

		bool beagleReturnsSums = true;

		if(beagleReturnsSums){
			results.lnL = siteLikesOut[0];
			results.d1 = siteD1Out[0];
			results.d2 = siteD2Out[0];

#ifdef OUTPUT_BEAGLE_SITELIKES
			string name = ofprefix + ".Bsitelikes.log";
			ofstream file(name.c_str(), ios::app);	
			file.precision(10);
			OutputBeagleSiteValues(file, theOp->derivatives);
			file.close();
			//mySummation = SumSiteValues(&siteLikesOut[0], (theOp->derivatives ? &siteD1Out[0] : NULL), (theOp->derivatives ? &siteD2Out[0] : NULL));
			//outman.UserMessage("mine = %.4f beag = %.4f", mySummation.lnL, results.lnL);
#endif
			}
		else
			results = SumSiteValues(&siteLikesOut[0], (theOp->derivatives ? &siteD1Out[0] : NULL), (theOp->derivatives ? &siteD2Out[0] : NULL));
		assert(results.lnL < 0.0 && results.lnL > -10.0e10);
		assert(results.d1 < 10.0e15 && results.d1 > -10.0e15);
		assert(results.d2 < 10.0e25 && results.d2 > -10.0e25);
#ifdef OUTPUT_OTHER_BEAGLE
		outman.DebugMessage("L\t%f\tD1\t%f\tD2\t%f", results.lnL, results.d1, results.d2);
#endif
#ifdef TERMINATE_AFTER_SCORE
		throw ErrorException("TERMINATING AFTER ONE SCORE");
#endif

#endif
		}
	return results;
	}

ScoreSet CalculationManager::SumSiteValues(const FLOAT_TYPE *sitelnL, const FLOAT_TYPE *siteD1, const FLOAT_TYPE *siteD2) const{
	FLOAT_TYPE lnL = 0.0, D1 = 0.0, D2 = 0.0;
	for(int i = 0;i < data->NChar();i++){
		assert(sitelnL[i] == sitelnL[i]);
		assert(sitelnL[i] <= 1e-6);

		if(sitelnL[i] > 0.0){
			outman.DebugMessage("BEAGLE-GARLI WARNING - site lnL > 0 : %e", sitelnL[i]);
			}

		assert(sitelnL[i] > -1e6);

		if(sitelnL[i] < -1e6){
			outman.DebugMessage("BEAGLE-GARLI WARNING - very small site lnL : %e", sitelnL[i]);
			}

		double actuallnL = sitelnL[i];
		if(actuallnL > 0.0)
			actuallnL = 0.0;
		if(actuallnL < -FLT_MAX)
			actuallnL = -FLT_MAX;			

		lnL += actuallnL * data->Count(i);
		
		if(siteD1 != NULL){
			assert(siteD1[i] == siteD1[i]);
//				assert(siteD1[i] > -1e15);
//				assert(siteD1[i] < 1e15);

			if(siteD1[i] > 1e15){
				outman.DebugMessage("BEAGLE-GARLI WARNING - very large site D1 : %e", siteD1[i]);
				}
			if(siteD1[i] < -1e15){
				outman.DebugMessage("BEAGLE-GARLI WARNING - very small site D1 : %e", siteD1[i]);
				}
			
			D1 += siteD1[i] * data->Count(i);
			}
		if(siteD2 != NULL){
			assert(siteD2[i] == siteD2[i]);
//				assert(siteD2[i] > -1e15);
//				assert(siteD2[i] < 1e15);
			
			if(siteD2[i] > 1e15){
				outman.DebugMessage("BEAGLE-GARLI WARNING - very large site D1 : %e", siteD2[i]);
				}
			if(siteD2[i] < -1e15){
				outman.DebugMessage("BEAGLE-GARLI WARNING - very small site D1 : %e", siteD2[i]);
				}

			D2 += siteD2[i] * data->Count(i);
			}
		}
	ScoreSet res = {lnL, D1, D2};
	return res;
	}

void CalculationManager::GetBeagleSiteLikelihoods(double *likes){
	beagleGetSiteLogLikelihoods(beagleInst, likes);
	}

void CalculationManager::OutputBeagleSiteValues(ofstream &out, bool derivs) const{
	//there is only one set of sitelikes stored by an instance, the most recent ones		
	vector<double> siteLikesOut(data->NChar());
	vector<double> siteD1Out(data->NChar());
	vector<double> siteD2Out(data->NChar());

	const int *counts = data->GetCounts();

	beagleGetSiteLogLikelihoods(beagleInst, &(siteLikesOut[0]));
	if(derivs)
		beagleGetSiteDerivatives(beagleInst, &(siteD1Out[0]), &(siteD2Out[0]));

	//ScoreSet mySummation = SumSiteValues(&siteLikesOut[0], (derivs ? &siteD1Out[0] : NULL), (derivs ? &siteD2Out[0] : NULL));
	//out << "mine\t" << mySummation.lnL << "\t" << mySummation.d1 << "\t" << mySummation.d2 << "\n";
	//outman.DebugMessage("myL\t%f\tD1\t%f\tD2\t%f", mySummation.lnL, mySummation.d1, mySummation.d2);
	//outman.UserMessage("mine = %.4f beag = %.4f", mySummation.lnL, results.lnL);

	for(int s = 0;s < data->GapsIncludedNChar();s++){
		int packed = data->Number(s);
		if(packed >= 0){
			out << s << "\t" << packed << "\t" << counts[packed] << "\t" << siteLikesOut[packed];
			if(derivs){
				out << "\t" << siteD1Out[packed] << "\t" << siteD2Out[packed];
				}
			}
		else 
			out << s << "\t" << packed << "\t" << counts[packed] << "\tNA";
		out << "\n";
		}
	}

//these are my old functions, not used when in beagle mode
void CalculationManager::CalcFullClaInternalInternal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, MODEL_FLOAT ***Lpmat,  MODEL_FLOAT ***Rpmat, const Model *mod){

	FLOAT_TYPE *Lpr = **Lpmat;
	FLOAT_TYPE *Rpr = **Rpmat;

	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	FLOAT_TYPE *dest=destCLA->arr;
	const FLOAT_TYPE *LCL=LCLA->arr;
	const FLOAT_TYPE *RCL=RCLA->arr;
	FLOAT_TYPE L1, L2, L3, L4, R1, R2, R3, R4;

	const int nRateCats = mod->NRateCats();
	const int nchar = data->NChar();
	const int *counts = data->GetCounts();

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

	const int *left_mult=LCLA->underflow_mult;
	const int *right_mult=RCLA->underflow_mult;
	int *undermult=destCLA->underflow_mult;

	for(int i=0;i<nchar;i++){
		undermult[i] = left_mult[i] + right_mult[i];
		}
	destCLA->rescaleRank = 2 + LCLA->rescaleRank + RCLA->rescaleRank;
	}

void CalculationManager::CalcFullCLATerminalTerminal(CondLikeArray *destCLA, const char *Ldata, const char *Rdata, MODEL_FLOAT ***Lpmat, MODEL_FLOAT ***Rpmat, const Model *mod){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.

	FLOAT_TYPE *Lpr = **Lpmat;
	FLOAT_TYPE *Rpr = **Rpmat;

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

void CalculationManager::CalcFullCLAInternalTerminal(CondLikeArray *destCLA, const CondLikeArray *LCLA, char *dat2, MODEL_FLOAT ***Lpmat, MODEL_FLOAT ***Rpmat, const Model *mod, const unsigned *ambigMap){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.

	FLOAT_TYPE *pr1 = **Lpmat;
	FLOAT_TYPE *pr2 = **Rpmat;

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


FLOAT_TYPE CalculationManager::GetScorePartialInternalRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, MODEL_FLOAT ***pmat, const Model *mod){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.

	const FLOAT_TYPE *prmat = **pmat;

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
FLOAT_TYPE CalculationManager::GetScorePartialTerminalRateHet(const CondLikeArray *partialCLA, MODEL_FLOAT ***pmat, const char *Ldata, const Model *mod){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	const FLOAT_TYPE *prmat = **pmat;
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
