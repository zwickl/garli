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

#ifndef CALCULATION_MANAGER
#define CALCULATION_MANAGER

#ifdef USE_BEAGLE
#include "beagle.h"
#endif

//extern ClaManager claMan;
//owned by each TreeNode
#include "utility.h"
#include "model.h"
class TreeNode;

struct ScoreSet{
	MODEL_FLOAT lnL;
	MODEL_FLOAT d1;
	MODEL_FLOAT d2;
	};

/*this is a generic matrix for use as a pmat or derivative matrix
multiple matrices for rates or whatever appear one after another*/
class TransMat{
	friend class TransMatSet;
	int nRates;
	int nStates;
	MODEL_FLOAT ***theMat;
public:

	TransMat(): theMat(NULL), nRates(-1), nStates(-1){}
	void Allocate(int numStates, int numRates){
		nStates = numStates;
		nRates = numRates;
		theMat = New3DArray<MODEL_FLOAT>(nRates, nStates, nStates);
		}
	~TransMat(){
		if(theMat)
			Delete3DArray<MODEL_FLOAT>(theMat);
		theMat = NULL;
		}
	MODEL_FLOAT ***GetMatrix(){
		return theMat;
		}
	};

/*package of pmat/d1mat/d2mat that will be assigned to a TransMatHolder (representing a branch)*/
class TransMatSet{
	friend class TransMatHolder;

	TransMat pmat;
	TransMat d1mat;
	TransMat d2mat;
public:
	void Allocate(int numStates, int numRates){
		pmat.Allocate(numStates, numRates);
		d1mat.Allocate(numStates, numRates);
		d2mat.Allocate(numStates, numRates);
		}
	MODEL_FLOAT ***GetPmatArray(){
		assert(pmat.theMat);
		return pmat.theMat;
		}
	MODEL_FLOAT ***GetD1Array(){
		assert(d1mat.theMat);
		return d1mat.theMat;
		}
	MODEL_FLOAT ***GetD2Array(){
		assert(d2mat.theMat);
		return d2mat.theMat;
		}
	};

/*analogous to a ClaHolder.  When a valid transmat set exists theMatSet will
keep a pointer to it.  When dirty it will be NULL*/
class TransMatHolder{
	friend class PmatManager;
	//friend class NodeClaManager;
	//friend class CalculationManager;

	//reference count
	int numAssigned;

	//these are the things that the holder must know to calculate valid transmats for itself
	Model *myMod;
	MODEL_FLOAT edgeLen;

	//will be NULL if dirty
	TransMatSet *theMatSet;

	TransMatHolder() : theMatSet(NULL), myMod(NULL), edgeLen(-1.0), numAssigned(0){};
	~TransMatHolder() {
		myMod = NULL;
		theMatSet = NULL;
		}
	void Reset(){numAssigned = 0; theMatSet = NULL; myMod = NULL;}
	
	MODEL_FLOAT ***GetPmatArray(){
		assert(theMatSet);
		assert(theMatSet->GetPmatArray()[0] && theMatSet->GetPmatArray()[0][0]);
		return theMatSet->GetPmatArray();
		}

	MODEL_FLOAT ***GetD1Array(){
		assert(theMatSet);
		return theMatSet->GetD1Array();
		}

	MODEL_FLOAT ***GetD2Array(){
		assert(theMatSet);
		return theMatSet->GetD2Array();;
		}

	const Model *GetConstModel() const{
		return myMod;
		}
	};

class PmatManager{
	int nRates;
	int nStates;
	int nMats;
	int nHolders;
	int maxUsed;
	
	//these are the actual matrices to be used in calculations, but will assigned to 
	//nodes (branches) via a PmatHolder.  There may be a limited number
	TransMatSet **allMatSets;

	//there will be enough of these such that every branch of every tree could
	//have a unique one, although many will generally be shared
	TransMatHolder *holders; 
	
	//pointers to the matrix sets
	vector<TransMatSet*> transMatSetStack;
	
	//indeces of free holders
	vector<int> holderStack;
public:
	PmatManager(int numMats, int numHolders, int numRates, int numStates) : nMats(numMats), nHolders(numHolders), nRates(numRates), nStates(numStates){
		allMatSets = new TransMatSet* [nMats];
		transMatSetStack.reserve(nMats);
		
		for(int i = nMats - 1;i >= 0;i--){
			allMatSets[i] = new TransMatSet;
			allMatSets[i]->Allocate(nStates, nRates);
			transMatSetStack.push_back(allMatSets[i]);
			}

		holders = new TransMatHolder[nHolders];
		holderStack.reserve(nHolders);
		for(int i=nHolders-1;i>=0;i--)
			holderStack.push_back(i);
		}

	~PmatManager(){
		transMatSetStack.clear();
		holderStack.clear();
		if(allMatSets != NULL){
			for(int i = 0;i < nMats;i++){
				delete allMatSets[i];
				}
			delete []allMatSets;
			}
		delete []holders;
		}

	const int GetNumStates() const { return nStates; }
	const int GetNumRates() const { return nRates; }
	bool TransMatIsAssigned(int index) const { return holders[index].theMatSet != NULL; }

	TransMatSet *GetTransMatSet(int index){
		if(holders[index].theMatSet == NULL)
			FillTransMatHolder(index);
		assert(holders[index].theMatSet != NULL);
		return holders[index].theMatSet;
		}

	void SetModel(int index, Model *mod){
		holders[index].myMod = mod;
		}

	void SetEdgelen(int index, FLOAT_TYPE e){
		holders[index].edgeLen = e;
		}

	const FLOAT_TYPE GetEdgelen(int index) const{
		return holders[index].edgeLen;
		}

	void GetEigenSolution(int index, ModelEigenSolution &sol) const{
		assert(holders[index].myMod);
		holders[index].myMod->GetEigenSolution(sol);
		}

//this includes pinv, which my scheme doesn't treat as a rate class per se
	void GetCategoryRatesForBeagle(int index, vector<FLOAT_TYPE> &r)const{
		holders[index].myMod->GetRateMultsForBeagle(r);
		}

//this includes pinv
	void GetCategoryWeightsForBeagle(int index, vector<FLOAT_TYPE> &p) const {
		holders[index].myMod->GetRateProbsForBeagle(p);
		}

	//need to include pinv, which is a separate rate as far as beagle is concerned, but not for Gar
	int NumRateCatsForBeagle(int index) const { 
		return holders[index].myMod->NumRateCatsForBeagle();
		}

	MODEL_FLOAT ***GetPmatArray(int index){
		assert(holders[index].numAssigned > 0);
		if(holders[index].theMatSet == NULL){
			FillTransMatHolder(index);
			//outman.DebugMessage("get, fill %d", index);
			}
		else{
			//outman.DebugMessage("get, %d already filled", index);
			}
		return holders[index].GetPmatArray();
		}
 
	MODEL_FLOAT ***GetD1MatArray(int index){
		assert(holders[index].numAssigned > 0);
		if(holders[index].theMatSet == NULL)
			FillTransMatHolder(index);
		return holders[index].GetD1Array();
		}

	MODEL_FLOAT ***GetD2MatArray(int index){
		assert(holders[index].numAssigned > 0);
		if(holders[index].theMatSet == NULL)
			FillTransMatHolder(index);
		return holders[index].GetD2Array();
		}

	void FillTransMatHolder(int index){
		//outman.DebugMessage("assign %d", index);
		assert(holders[index].theMatSet == NULL);
		holders[index].theMatSet = AssignFreeTransMatSet();
		}

	const TransMatHolder *GetTransMatHolder(int index) const{
		//outman.DebugMessage("get %d", index);
		return &holders[index];
		}

	TransMatHolder *GetMutableTransMatHolder(int index){
		//outman.DebugMessage("get mut %d", index);
		return &holders[index];
		}

	const Model *GetModelForTransMatHolder(int index)const{
		return holders[index].GetConstModel();
		}

	Model *GetMutableModelForTransMatHolder(int index)const{
		return holders[index].myMod;
		}

	//this will be called by a node (branch) and will return a previously unused holder index
	int AssignTransMatHolder(){
		assert(holderStack.size() > 0);
		int index=holderStack[holderStack.size()-1];
		assert(holders[index].numAssigned == 0);
		IncrementTransMatHolder(index);
		holderStack.pop_back();
		return index;
		}

	TransMatSet* AssignFreeTransMatSet(){
		//TODO
		//DEBUG need to figure out how this will work with beagle, or if recycling is even necessary with transmats
		//if(claStack.empty() == true) RecycleClas();

		assert(! transMatSetStack.empty());
		TransMatSet *mat=transMatSetStack[transMatSetStack.size()-1];
		transMatSetStack.pop_back();
		if(nMats - (int)transMatSetStack.size() > maxUsed) 
			maxUsed = nMats - (int)transMatSetStack.size();
		return mat;
		}

	void IncrementTransMatHolder(int index){
		//outman.DebugMessage("inc %d to %d", index, holders[index].numAssigned + 1);
		holders[index].numAssigned++;
		}

	void DecrementTransMatHolder(int index){
		TransMatHolder *thisHold = &holders[index];

		assert(index != -1);
		assert(thisHold->numAssigned != 0);
		//outman.DebugMessage("dec %d to %d", index, holders[index].numAssigned - 1);
		if(thisHold->numAssigned == 1){
			//if the count has fallen to zero, reclaim the holder
			//outman.DebugMessage("reclaim1 %d", index);
			holderStack.push_back(index);
			if(holders[index].theMatSet != NULL){
				//if there is a valid matrix set in the holder, reclaim it too
				//outman.DebugMessage("reclaim1 set from %d", index);
				assert(find(transMatSetStack.begin(), transMatSetStack.end(), holders[index].theMatSet) == transMatSetStack.end());
				transMatSetStack.push_back(thisHold->theMatSet);
				}
			thisHold->Reset();
			}
		else{
			thisHold->numAssigned--; 
			//TODO
			//DEBUG - what happens with this and beagle
			//this is important!
//			holders[index].tempReserved=false;
			}
		}

	int SetTransMatDirty(int index){
		//there are only two options here:
		//1. transmatSet is being made dirty, and only node node points to it 
		//	->null the holder's transmatSet pointer and return the same index
		//2. transmatSet is being made dirty, and multiple nodes point to it
		//	->remove this transmatSet from the holder (decrement) and assign a new one	
		assert(index != -1);

		TransMatHolder *thisHold = &holders[index];
		//outman.DebugMessage("dirty %d", index);
		if(thisHold->numAssigned==1){
			if(thisHold->theMatSet != NULL){
				//TODO
				//holders[index].SetReclaimLevel(0);
				//outman.DebugMessage("reclaim2 %d", index);
				
				transMatSetStack.push_back(thisHold->theMatSet);
				thisHold->theMatSet=NULL;
				}
			}
		else{
			DecrementTransMatHolder(index);
			assert(holderStack.size() > 0);
			index=holderStack[holderStack.size()-1];
			holderStack.pop_back();
			IncrementTransMatHolder(index);
			}
		return index;
		}
	};

class NodeClaManager{
	static ClaManager *claMan;
	static PmatManager *pmatMan;

public:
	//Indeces should always be valid after initial assignment
	//Dirtying is reassignment to a new ClaHolder index (initially containing no assigned CondlikeArray index)
	//or NULLing of the CondlikeArray pointer from the current ClaHolder if only one node refers to it
	int downHolderIndex;
	int ULHolderIndex;	
	int URHolderIndex;
	int transMatIndex;

	NodeClaManager() : downHolderIndex(-1), ULHolderIndex(-1), URHolderIndex(-1), transMatIndex(-1){};

	static void SetClaManager(ClaManager *cMan) {
		NodeClaManager::claMan = cMan;
		}

	static void SetPmatManager(PmatManager *pMan) {
		NodeClaManager::pmatMan = pMan;
		}

	//tips will have the cla holders all -1, but a valid transmat index
	bool IsAllocated() const {
		return !(downHolderIndex < 0 && ULHolderIndex < 0 && URHolderIndex < 0 && transMatIndex < 0);
		}

	void SetHolders(int d, int UL, int UR, int p){
		//this might be called to copy changes that were made to cla assignments NOT through the
		//NodeClaManager. So these can already be set to some other valid value
/*
		assert(downHolderIndex < 0);
		assert(ULHolderIndex < 0);
		assert(URHolderIndex < 0);
		assert(transMatIndex < 0);
*/
		downHolderIndex = d;
		ULHolderIndex = UL;	
		URHolderIndex = UR;
		transMatIndex = p;
		}

	void ObtainTransMatHolders(){
		assert(transMatIndex < 0);
		transMatIndex = pmatMan->AssignTransMatHolder();
		}

	void ObtainClaAndTransMatHolders(){
		assert(downHolderIndex < 0);
		assert(ULHolderIndex < 0);
		assert(URHolderIndex < 0);
		assert(transMatIndex < 0);
		downHolderIndex = claMan->AssignClaHolder();
		ULHolderIndex = claMan->AssignClaHolder();
		URHolderIndex = claMan->AssignClaHolder();
		transMatIndex = pmatMan->AssignTransMatHolder();
		}

	void SetDirtyUpRight(){
		URHolderIndex = claMan->SetDirty(URHolderIndex);
		}

	void SetDirtyUpLeft(){
		ULHolderIndex = claMan->SetDirty(ULHolderIndex);
		}

	void SetDirtyDown(){
		downHolderIndex = claMan->SetDirty(downHolderIndex);
		}

	void SetDirtyAll(){
		URHolderIndex = claMan->SetDirty(URHolderIndex);
		ULHolderIndex = claMan->SetDirty(ULHolderIndex);
		downHolderIndex = claMan->SetDirty(downHolderIndex);
		SetTransMatDirty();
		}

	void SetTransMatDirty(){
		if(transMatIndex >= 0)
			transMatIndex = pmatMan->SetTransMatDirty(transMatIndex);
		}

	void SetDependenciesUL(int depIndex1, int pDepIndex1, int depIndex2, int pDepIndex2){
		claMan->SetHolderDependencies(ULHolderIndex, depIndex1, pDepIndex1, depIndex2, pDepIndex2);
		}

	void SetDependenciesUR(int depIndex1, int pDepIndex1, int depIndex2, int pDepIndex2){
		claMan->SetHolderDependencies(URHolderIndex, depIndex1, pDepIndex1, depIndex2, pDepIndex2);
		}

	void SetDependenciesDown(int depIndex1, int pDepIndex1, int depIndex2, int pDepIndex2){
		claMan->SetHolderDependencies(downHolderIndex, depIndex1, pDepIndex1, depIndex2, pDepIndex2);
		}

	void SetTransMat(Model *m, FLOAT_TYPE e){
		//outman.DebugMessage("up pdeps %d", transMatIndex);

		pmatMan->SetModel(transMatIndex, m);
		pmatMan->SetEdgelen(transMatIndex, e);
		}

	inline CondLikeArray *GetClaDown(){
		bool calc = true;
		if(claMan->IsDirty(downHolderIndex)){
			if(calc==true){
				assert(0);
				//DEBUG - may have to nix this auto calc behavior when dirty with the new management system
				//Need to somehow calculate this here
				//ConditionalLikelihoodRateHet(DOWN, nd);
				}
			else claMan->FillHolder(downHolderIndex, 1);
			}

	//	if(memLevel > 1) claMan->ReserveCla(downClaIndex);
		return claMan->GetCla(downHolderIndex);
		}

	inline CondLikeArray *GetClaUpLeft(){
		bool calc = true;
		if(claMan->IsDirty(ULHolderIndex)){
			if(calc==true){
				//DEBUG
				assert(0);
				//Need to somehow calculate this here
				//ConditionalLikelihoodRateHet(UPLEFT, nd);
				}
			else claMan->FillHolder(ULHolderIndex, 2);
			}

	//	if(memLevel > 0) claMan->ReserveCla(nd->claIndexUL);

		return claMan->GetCla(ULHolderIndex);
		}

	inline CondLikeArray *GetClaUpRight(){
		bool calc = true;
		if(claMan->IsDirty(URHolderIndex)){
			if(calc==true){
				//DEBUG
				assert(0);
				//Need to somehow calculate this here
				//ConditionalLikelihoodRateHet(UPRIGHT, nd);
				}
			else claMan->FillHolder(ULHolderIndex, 2);
			}
	//	if(memLevel > 0) claMan->ReserveCla(nd->claIndexUR);
		return claMan->GetCla(ULHolderIndex);
		}

	inline TransMatSet *GetTransMatSet() const {
		return pmatMan->GetTransMatSet(transMatIndex);
		}

	void CopyHolderIndeces(const NodeClaManager *from, bool remove){
		if(remove){
			StripHolders();
			}
		else{
			int poo=2;
			}
		downHolderIndex = from->downHolderIndex;
		ULHolderIndex = from->ULHolderIndex;
		URHolderIndex = from->URHolderIndex;
		transMatIndex = from->transMatIndex;

		if(downHolderIndex >= 0){
			assert(claMan->GetNumAssigned(downHolderIndex) > 0);
			assert(claMan->GetNumAssigned(ULHolderIndex) > 0);
			assert(claMan->GetNumAssigned(URHolderIndex) > 0);
			claMan->IncrementCla(downHolderIndex);
			claMan->IncrementCla(ULHolderIndex);
			claMan->IncrementCla(URHolderIndex);
			}
		if(transMatIndex >= 0)
			pmatMan->IncrementTransMatHolder(transMatIndex);
		}

	void StripHolders(){
		if(downHolderIndex >= 0){
			claMan->DecrementCla(downHolderIndex);
			claMan->DecrementCla(ULHolderIndex);
			claMan->DecrementCla(URHolderIndex);

			downHolderIndex = -1;
			ULHolderIndex = -1;
			URHolderIndex = -1;
			}
		else{
			assert(ULHolderIndex < 0 && ULHolderIndex < 0);
			}
		if(transMatIndex >= 0){
			pmatMan->DecrementTransMatHolder(transMatIndex);
			transMatIndex = -1;
			}
		}
	};

class ClaOperation{
	friend class CalculationManager;
	friend class NodeOperation;
public:
	int destCLAIndex;
	int childCLAIndex1;
	int childCLAIndex2;
	int transMatIndex1;
	int transMatIndex2;
	int depLevel;
public:
	ClaOperation(): destCLAIndex(-1), childCLAIndex1(-1), childCLAIndex2(-1), transMatIndex1(-1), transMatIndex2(-1), depLevel(-1){};	
	ClaOperation(int d, int cla1, int cla2, int pmat1, int pmat2, int dLevel): destCLAIndex(d), childCLAIndex1(cla1), childCLAIndex2(cla2), transMatIndex1(pmat1), transMatIndex2(pmat2), depLevel(dLevel){};	
	ClaOperation(const ClaOperation &from){
		destCLAIndex = from.destCLAIndex;
		childCLAIndex1 = from.childCLAIndex1;
		childCLAIndex2 = from.childCLAIndex2;
		transMatIndex1 = from.transMatIndex1;
		transMatIndex2 = from.transMatIndex2;
		depLevel = from.depLevel;
		}
	};

class TransMatOperation{
	friend class CalculationManager;
	friend class NodeOperation;
	friend class BranchOperation;

	int destTransMatIndex;
	FLOAT_TYPE edgeLength;
	int modelIndex; // ~eigen solution
	bool calcDerivs;
public:
	TransMatOperation() : destTransMatIndex(-1), modelIndex(-1), edgeLength(-1.0), calcDerivs(false){};
	TransMatOperation(int ind, int modInd, FLOAT_TYPE len, bool d) : destTransMatIndex(ind), modelIndex(modInd), edgeLength(len), calcDerivs(d){};
	TransMatOperation(const TransMatOperation &from){
		destTransMatIndex = from.destTransMatIndex;
		edgeLength = from.edgeLength;
		modelIndex = from.modelIndex;
		calcDerivs = from.calcDerivs;
		}
	};

class NodeOperation{
	friend class CalculationManager;
	ClaOperation claOp;
	TransMatOperation transOp1;
	TransMatOperation transOp2;
public:
	NodeOperation();
	NodeOperation(const NodeOperation &from){
		claOp = from.claOp;
		transOp1 = from.transOp1;
		transOp2 = from.transOp2;
		}
	NodeOperation(const ClaOperation &c, const TransMatOperation &t1, const TransMatOperation &t2){
		claOp = c;
		transOp1 = t1;
		transOp2 = t2;
		}
	bool operator <(const NodeOperation &rhs){
		return claOp.depLevel < rhs.claOp.depLevel;
		}
	};

class ScoringOperation{
	friend class CalculationManager;
	friend class BranchOperation;

	//my functions don't write the final cla, but Beagle allows for specifying one
	int destClaIndex;
	int childClaIndex1;
	int childClaIndex2;
	int transMatIndex;
	bool derivatives;
public:
	ScoringOperation() : destClaIndex(-1), childClaIndex1(-1), childClaIndex2(-1), transMatIndex(-1), derivatives(false){};
	ScoringOperation(int dest, int child1, int child2, int pmat, bool d) : destClaIndex(dest), childClaIndex1(child1), childClaIndex2(child2), transMatIndex(pmat), derivatives(d){};
	ScoringOperation(const ScoringOperation &from){
		destClaIndex = from.destClaIndex;
		childClaIndex1 = from.childClaIndex1;
		childClaIndex2 = from.childClaIndex2;
		transMatIndex = from.transMatIndex;
		derivatives = from.derivatives;		
		}
	};

class BranchOperation{
	friend class CalculationManager;
	ScoringOperation scrOp;
	TransMatOperation transOp;
public:
	BranchOperation(const BranchOperation &from){
		scrOp = from.scrOp;
		transOp = from.transOp;
		}
	};

class BlockingOperationsSet{
public:
	list<ClaOperation> claOps;
	list<TransMatOperation> pmatOps;
	~BlockingOperationsSet(){
		claOps.clear();
		pmatOps.clear();
		}
	};

class NewBlockingOperationsSet{
public:
	list<NodeOperation> nodeOps;
	list<BranchOperation> brOps;
	~NewBlockingOperationsSet(){
		nodeOps.clear();
		brOps.clear();
		}
	};
/*
bool MyOpLessThan(const BlockingOperationsSet &lhs, const BlockingOperationsSet &rhs){
//	int d1 = lhs.claOps.size() == 0 ? -1 : lhs.claOps.begin()->depLevel;
//	int d2 = rhs.claOps.size() == 0 ? -1 : rhs.claOps.begin()->depLevel;
	if(lhs.claOps.size() == 0)
		return false;
	else if(rhs.claOps.size() == 0)
		return true;
	return 
		(lhs.claOps.begin()->depLevel < rhs.claOps.begin()->depLevel);
	}
*/
/*
This does all of the actual calculations and operations once they have been figured out at the Tree level.  It only understands about indeces
of CLA holders and pmat holders, and constructs the sets of operations to be calculated.

1. Tree funcs sweeps over the tree structure when necessary, using NodeClaManagers to dirty clas and set up the dependencies between different clas.
   Dirtying happens when blens change, topology changes or model changes (many clas get dirtied).  Topology changes don't change many dependencies
   directly, but dirtying may change the Holder indecies if multiple trees are pointing to it, and thereby change many dependencies.

2. Tree asks the CalcManager for either a lnL or lnL + derivs

3. CalcManager uses the holder dependencies to deduce what operations need to happen.

4. Operations are accumulated by the CalculationManager and carried out either by passing raw matrices and clas locally, or by passing 
   operations by indecies to Beagle.
*/

class CalculationManager{
	//the claManager may eventually be owned by calcManager, but for now it can interact indirectly
	//with the global manager, as everything did previously
	static ClaManager *claMan;
	static PmatManager *pmatMan;
	static const SequenceData *data;
	list<BlockingOperationsSet> operationSetQueue;
	list<ScoringOperation> scoreOps;
	
//	list<NewBlockingOperationsSet> newOperationSetQueue;

public:
	bool useBeagle;
	bool termOnBeagleError;
	bool rescaleBeagle;
	bool singlePrecBeagle;
	bool gpuBeagle;
	int beagleInst;
	string ofprefix;

	CalculationManager(){
		useBeagle = false; 
		gpuBeagle = false;
		beagleInst = -1; 
		termOnBeagleError = true;
		singlePrecBeagle = false;
		rescaleBeagle = false;
		}

	static void SetClaManager(ClaManager *cMan) {
		CalculationManager::claMan = cMan;
		}

	static void SetPmatManager(PmatManager *pMan) {
		CalculationManager::pmatMan = pMan;
		}

	static void SetData(SequenceData *dat){
		CalculationManager::data = dat;
		}
	
	void SetOfprefix(string &prefix){
		ofprefix = prefix;
		}

	//THE FOLLOWING DO NOT DEPEND ON BEAGLE BEING AVAILABLE

	//pass the effective root node (any non-tip), it will figure out which clas to combine
//	FLOAT_TYPE CalculateLikelihood(const TreeNode *effectiveRoot);

	//pass the node on the upper end of the branch.  may be terminal
//	ScoreSet CalculateDerivatives(const TreeNode *topOfBranch);

	ScoreSet CalculateLikelihoodAndDerivatives(const TreeNode *topOfBranch, bool calcDerivs);
	
	/*called by both CalculateLikelihood and CalculateDerivatives. interprets the focal node 
	depending on whether derivatives are requested*/
	void DetermineRequiredOperations(const TreeNode *focalNode, bool derivatives);

	/*called by DetermineRequiredOperations. using dependecies pre-set in the ClaHolders, tracks backwards recursively to 
	determine all of the pmats that must be calculated (TransMatOps) and the clas that must be combined (ClaOpts) to get 
	the initially passed holder number.  Note that at this point any notion of a tree or directionality is out of the 
	picture, only dependencies remain*/
	int AccumulateOpsOnPath(int holderInd, list<NodeOperation> &nodeOps);

	//always combine two clas and one pmat
	void PerformClaOperation(const ClaOperation *theOp);
	//calculate a whole set of cla operations at once
	void PerformClaOperationBatch(const list<ClaOperation> &theOps);

	//calculate a pmat, or a pmat, d1mat and d2mat
	void PerformTransMatOperation(const TransMatOperation *theOp);
	//calculate a whole set of transmat operations at once
	void PerformTransMatOperationBatch(const list<TransMatOperation> &theOps);
	//call back to the underlying models to calculate the pmats, then send them to beagle
	void SendTransMatsToBeagle(const list<TransMatOperation> &theOps);

	//calculate and return either the lnL alone, or the lnl, D1 and D2
	ScoreSet PerformScoringOperation(const ScoringOperation *theOp);

	//This reports back whether Beagle ended up using single precision.  It may be because it was specifically specified,
	//or just implied by other settings such as GPU.  Should be called AFTER SetBeagleDetails and InitializeBeagle, at which
	//point the details of the instance are definitely known
	bool IsSinglePrecision() {return singlePrecBeagle;}
	
	//this is just for debugging
	void OutputOperationsSummary() const;

	//these are just duplicates of the functions originally in Tree.  They are not used with beagle
	void CalcFullClaInternalInternal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA,  MODEL_FLOAT ***Lpr,  MODEL_FLOAT ***Rpr, const Model *mod);
	void CalcFullCLATerminalTerminal(CondLikeArray *destCLA, const char *Ldata, const char *Rdata,  MODEL_FLOAT ***Lpr,  MODEL_FLOAT ***Rpr, const Model *mod);
	void CalcFullCLAInternalTerminal(CondLikeArray *destCLA, const CondLikeArray *LCLA, char *data2,  MODEL_FLOAT ***pr1,  MODEL_FLOAT ***pr2, const Model *mod, const unsigned *ambigMap);
	FLOAT_TYPE GetScorePartialTerminalRateHet(const CondLikeArray *partialCLA,  MODEL_FLOAT ***prmat, const char *Ldata, const Model *mod);
	FLOAT_TYPE GetScorePartialInternalRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA,  MODEL_FLOAT ***prmat, const Model *mod);
	
	ScoreSet SumSiteValues(const FLOAT_TYPE *sitelnL, const FLOAT_TYPE *siteD1, const FLOAT_TYPE *siteD2){
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

	void Finalize(){
#ifdef USE_BEAGLE
		if(beagleInst > 0)
			beagleFinalizeInstance(beagleInst);
#endif
		}

#ifdef USE_BEAGLE
	void SetBeagleDetails(bool gpu, bool singlePrec, bool rescale, string &prefix){
		useBeagle = true;
		if(gpu){//assuming that all GPU is SP for now
			gpuBeagle = true;
			singlePrecBeagle = true;
			rescaleBeagle = true;
			}
		if(singlePrec){
			singlePrecBeagle = true;
			}
		if(rescale)
			rescaleBeagle = true;
		this->ofprefix = prefix;
		}
	
	//BEAGLE SPECIFIC FUNCS

	void InitializeBeagle(int nnod, int nClas, int nHolders, int nstates, int nchar, int nrates);

	void SendTipDataToBeagle();

	//interpret beagle error codes, output a message and bail if termOnBeagleError == true
	void CheckBeagleReturnValue(int err, const char *funcName) const;

	//simple report
	void OutputInstanceDetails(const BeagleInstanceDetails *det) const;

	/*beagle treats cla indeces 0 -> NTax - 1 as the tips, then above that the other internals.  My internal cla indexing
	starts at zero, and the negatives are just (-taxnum) NOT starting at zero.  So, convert like this: 
			T4 T3 T2 T1 int1 int2 int3
		me	-4 -3 -2 -1  0    1     2
		Be	 3  2  1  0  4    5     6	*/
	int PartialIndexForBeagle(int ind){
		return (ind < 0 ? ((-ind) - 1) : ind + data->NTax());
		}

	/*beagle doesn't know anything about pmat vs deriv mats in terms of storage.  So, I keep track of TransMatSets, which contain one of each
	these will be transformed into beagle indeces such that the pmat, d1 and d2 mats use three consecutive beagle mat indeces
	thus, for my TransMatSet 0, the beagle pmat index is 0, d1 index is 1, d2 is 2.  TransMatSet 1 is 3, 4, 5, etc. */
	int PmatIndexForBeagle(int ind){
		return (ind * 3);
		}
	int D1MatIndexForBeagle(int ind){
		return (ind * 3) + 1;
		}
	int D2MatIndexForBeagle(int ind){
		return (ind * 3) + 2;
		}

	//For scale arrays my indexing scheme and Beagle's happen to be the same
	void AccumulateRescalers(int destIndex, int childIndex1, int childIndex2);

	//for testing/debugging
	void SendClaToBeagle(int num){
		CheckBeagleReturnValue(
			beagleSetPartials(
				beagleInst,
				PartialIndexForBeagle(num),
				claMan->GetCla(num)->arr), 
			"beagleSetPartials");
		}

	void OutputBeagleSiteLikelihoods(ofstream &likes, vector<double> &siteLikesOut){
		const int *count = data->GetCounts();
		int num = 0;
		
		//for(vector<double>::iterator it = siteLikesOut.begin();it != siteLikesOut.end();it++)
		for(int c = 0;c < data->NChar();c++)
			for(int co = 0;co < count[c];co++)
				likes << c << "\t" << co << "\t" << siteLikesOut[c] << "\n";
		likes.precision(11);
		likes << SumSiteValues(&(siteLikesOut[0]), NULL, NULL).lnL << "\n";
		likes.close();
		}

	void OutputBeagleSiteValues(ofstream &out, bool derivs){
		//there is only one set of sitelikes stored by an instance, the most recent ones		
		vector<double> siteLikesOut(data->NChar());
		vector<double> siteD1Out(data->NChar());
		vector<double> siteD2Out(data->NChar());

		const int *counts = data->GetCounts();

		beagleGetSiteLogLikelihoods(beagleInst, &(siteLikesOut[0]));
		if(derivs)
			beagleGetSiteDerivatives(beagleInst, &(siteD1Out[0]), &(siteD2Out[0]));

		for(int s = 0;s < data->GapsIncludedNChar();s++){
			int packed = data->Number(s);
			out << s << "\t" << packed << "\t" << counts[packed] << "\t" << siteLikesOut[packed];
			if(derivs){
				out << "\t" << siteD1Out[packed] << "\t" << siteD2Out[packed];
				}
			out << "\n";
			}
		}

	void UpdateAllConditionals();

	void OutputBeagleTransMat(int index);
#endif 
	};
#endif
