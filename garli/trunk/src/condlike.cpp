// GARLI version 2.0 source code
// Copyright 2005-2011 Derrick J. Zwickl
// email: garli.support@gmail.com
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

#include "defs.h"
#include "condlike.h"
#include "clamanager.h"
#include "utility.h"


CondLikeArray::~CondLikeArray(){
	//with partitioning, the entire allocation is managed and deleted by the
	//condlikearrayset, so don't delete anything here
	arr = NULL;
	underflow_mult = NULL;
/*	if( arr ){
		delete []arr;
		arr=NULL;
		 }
	if(underflow_mult!=NULL) delete []underflow_mult;
*/
}

void CondLikeArray::Allocate( int ind, int nk, int ns, int nr /* = 1 */ ){
	if( arr ){
		delete []arr;
		arr=NULL;
		}
	nrates = nr;
	nsites = ns;
	nstates = nk;
	index = ind;
	arr=new FLOAT_TYPE[nk*nr*ns];
	if(arr==NULL){
		throw ErrorException("GARLI had a problem allocating memory!  Try reducing the availablememory setting.");
		}
	//DJZ 1-5-05 - the underflow mult used to be ns*nr in length,
	//but only needs to be ns
	underflow_mult=new int[ns];
	if(underflow_mult==NULL){
		throw ErrorException("GARLI had a problem allocating memory!  Try reducing the availablememory setting.");
		}
	}


CondLikeArrayHolder *ClaManager::GetMutableHolder(int index) {
	assert(index > -1 && index < numHolders);
	return &holders[index];
}

CondLikeArraySet* ClaManager::AssignFreeCla(){
	#ifdef CLA_DEBUG
	ofstream deb("cladebug.log", ios::app);
	#endif

	if(claStack.empty() == true) RecycleClas();
	
	CondLikeArraySet *arr=claStack[claStack.size()-1];
	
	assert(arr != NULL);
	claStack.pop_back();
	if(numClas - (int)claStack.size() > maxUsed) maxUsed=numClas - (int)claStack.size();
	
	return arr;
	}

void ClaManager::RecycleClas(){
	int numReclaimed=0;
	for(int i=0;i<numHolders;i++){
		if(holders[i].theSet != NULL){
			if(holders[i].GetReclaimLevel() == 2 && holders[i].tempReserved == false && holders[i].reserved == false){
				claStack.push_back(holders[i].theSet);
				holders[i].SetReclaimLevel(0);
				holders[i].theSet=NULL;
				numReclaimed++;
				}
			}
		if(memLevel < 2) 
			if(numReclaimed > 50) 
				return;
		}
	if(numReclaimed > 10) 
		return;
	for(int i=0;i<numHolders;i++){
		if(holders[i].theSet != NULL){
			if((holders[i].GetReclaimLevel() == 1 && holders[i].tempReserved == false && holders[i].reserved == false)){
				claStack.push_back(holders[i].theSet);
				holders[i].SetReclaimLevel(0);
				holders[i].theSet=NULL;
				numReclaimed++;
				}
			}
		if(numReclaimed == 20) 
			return;
		}
	if(numReclaimed==0){
		//I changed this for some reason in r1030 (April 7, 2011, which means that it was in the 2.0 release) to the ErrorException, 
		//which I should not have.  Throwing 2 will dirty the entire tree, which needs to happen in some cases when there are too many
		//CLAs in use during blen opt and certain memlevels are in effect. The error throw was causing it to bail when it shouldn't have.  
		//Reverting this June 15, 2011
		throw(2);
		//throw ErrorException("Ran out of conditional likelihood arrays. This should not really happen, but try increasing availablememory setting");
		}
	assert(numReclaimed > 0);
	}

void CondLikeArraySet::Allocate() {
	unsigned size = 0, usize = 0;
	for(vector<CondLikeArray *>::iterator cit = theSets.begin();cit != theSets.end();cit++){
		size += (*cit)->RequiredSize();
		usize += (*cit)->NChar();
		}
	try{
		rawAllocation = new FLOAT_TYPE[size];
		}
	catch(std::bad_alloc){
		throw ErrorException("Problem allocating cond. likelihood array (len = %d). Out of mem?\n\tNote: to use > 4GB of memory, you will need a 64-bit version of GARLI.", size);
		}
	try{
		rawUnder = new int[usize];
		}
	catch(std::bad_alloc){
		throw ErrorException("Problem allocating underflow multiplier array (len = %d). Out of mem?", usize);
		}
	unsigned offset = 0, uoffset = 0;
	for(vector<CondLikeArray *>::iterator cit = theSets.begin();cit != theSets.end();cit++){
		(*cit)->Assign(&rawAllocation[offset], &rawUnder[uoffset]);
		offset += (*cit)->RequiredSize();
		uoffset += (*cit)->NChar();
		}
	}

//tell the manager that another tree (node) now points to this holder.  
void ClaManager::IncrementHolder(int index) {
	assert(index > -1 && index < numHolders);
	holders[index].numAssigned++;
	if (debug_clas) {
		if (this->IsHolderReserved(index))
			outman.DebugMessage("incr(res)\t%d\t->\t%d", index, holders[index].numAssigned);
		else
			outman.DebugMessage("incr(noRes)\t%d\t->\t%d", index, holders[index].numAssigned);
	}
	CheckClaHolders();
}

//tell the manager that one fewer tree (node) now points to this holder. If that brings the count to zero, then
//reclaim the holder, as well as any cla in it
void ClaManager::DecrementHolder(int index) {
	assert(index > -1 && index < numHolders);
	CheckClaHolders();
	if (debug_clas)
		outman.DebugMessage("decr\t%d\t%d\t->\t%d", index, holders[index].numAssigned, holders[index].numAssigned - 1);

	if (holders[index].numAssigned == 1) {
		if (debug_clas) {
			assert(find(holderStack.begin(), holderStack.end(), index) == holderStack.end());
		}

		_ReclaimH(index);
		
		if (debug_clas) {
			ReportClaTotals("reclaimH", index);
		}
	}
	else {
		assert(holders[index].numAssigned != 0);
		holders[index].numAssigned--;
		//DEBUG - Why was this important?  Should this be here for beagle?
		//this is important!
		//Yes, I think that this should be done, since the holder should only be temp reseved if it
		//was in current use, and if it is being decremented then it must no longer be.  Normal
		//reserations will persist.
		RemoveTempReservation(index);
	}
}

//this will assign an unused holder index to a node that did not previously have one (no argument)
int ClaManager::AssignFreeClaHolder() {
	int newIndex = _GetFreeH();
	if (debug_clas) {
		ReportClaTotals("assignedH\0", newIndex);
	}
	return newIndex;
}

//this is essentially identical to AssignNewClaHolder, but takes an old index
//back. it is written separately for clarity 
int ClaManager::TradeInClaHolder(int oldIndex) {
	int newIndex = _GetFreeH();
	DecrementHolder(oldIndex);

	if (debug_clas) {
		ReportClaTotals("tradein\0", newIndex);
	}
	return newIndex;
}

void ClaManager::FillHolder(int index, int dir) {
	holders[index].theSet = AssignFreeCla();
	holders[index].reclaimLevel = dir;
}

int ClaManager::GetReclaimLevel(int index) {
	if (holders[index].theSet == NULL) return -1;
	return holders[index].GetReclaimLevel();
}

void ClaManager::SetReclaimLevel(int index, int lvl) {
	assert(index > -1);
	if (holders[index].theSet == NULL)
		assert(0);
	//return;
	holders[index].SetReclaimLevel(lvl);
}

void ClaManager::ReserveCla(int index, bool temp/*=true*/) {
	if (temp == true) holders[index].tempReserved = true;
	else holders[index].reserved = true;
}

void ClaManager::UnreserveCla(int index) {
	//		holders[index].tempReserved=false;
	holders[index].reserved = false;
	if (memLevel>1)
		holders[index].SetReclaimLevel(1);
}

void ClaManager::ReclaimSingleCla(int index) {
	//this simply removes the cla from a holder.  It is equivalent to just
	//dirtying it if only a single tree shares the holder
	if (holders[index].theSet == NULL) return;
	claStack.push_back(holders[index].theSet);
	holders[index].SetReclaimLevel(0);
	holders[index].theSet = NULL;
}

void ClaManager::CountClaTotals(int &clean, int &tempres, int &res, int &assigned) {
	for (int i = 0; i<numHolders; i++) {
		if (holders[i].theSet != NULL) clean++;
		if (holders[i].tempReserved == true) tempres++;
		if (holders[i].reserved == true) res++;
	}
	assigned = numHolders - (int)holderStack.size();
}

int ClaManager::GetClaNumber(int index) {
	//this is ugly, but should only be called for debugging
	if (holders[index].theSet == NULL) return -1;
	for (int i = 0; i<numClas; i++)
		if (holders[index].theSet == allClas[i]) return i;
	assert(0);
	return -1;
}

int ClaManager::CountClasInUse(int recLevel) {
	int num = 0;
	for (int i = 0; i<numHolders; i++) {
		if (holders[i].theSet != NULL)
			if (holders[i].GetReclaimLevel() == recLevel) num++;
	}
	return num;
}

CondLikeArray *ClaManager::GetCla(int index, int modnum/*=0*/) {
	assert(holders[index].theSet != NULL);
	return GetClaSet(index)->theSets[modnum];
}

CondLikeArraySet *ClaManager::GetClaSet(int index) {
	assert(holders[index].theSet != NULL);
	return holders[index].theSet;
}

const CondLikeArrayHolder *ClaManager::GetHolder(int index) {
	assert(index > -1 && index < numHolders);
	return &holders[index];
}

bool ClaManager::IsDirty(int index) {
	//dirtyness is now synonymous with a null cla pointer in the holder
	assert(index > -1);
	return (holders[index].theSet == NULL);
}


//this marks this holder as one that should not have its cla reclaimed during recycling
//TEMP RESERVATIONS are nodes that are necessary for ongoing operations - destination arrays
//and immediate dependencies
void ClaManager::TempReserveCla(int index) {
	assert(index > -1 && index < numHolders);
	holders[index].tempReserved = true;
	if (debug_clas) {
		ReportClaTotals("tempreserved", index);
	}
}

int ClaManager::SetDirty(int index) {

	//there are only two options here:
	//1. Cla is being made dirty, and only node node points to it 
	//	->null the holder's cla pointer and return the same index
	//2. Cla is being made dirty, and multiple nodes point to it
	//	->remove this node from the holder (decrement) and assign a new one	

	assert(index != -1);

	if (holders[index].numAssigned == 1) {
		if (holders[index].theSet != NULL) {
			holders[index].SetReclaimLevel(0);
			claStack.push_back(holders[index].theSet);
			holders[index].theSet = NULL;
		}
	}
	else {
		//DecrementCla(index);
		DecrementHolder(index);
		index = holderStack[holderStack.size() - 1];
		holderStack.pop_back();
		IncrementHolder(index);
	}
	return index;
}

void ClaManager::IncrementCla(int index) {
	assert(0);
	holders[index].numAssigned++;
}

void ClaManager::DecrementCla(int index) {
	//Deprecated - DecrementHolder should be used now
	assert(0);

	assert(index != -1);
	if (holders[index].numAssigned == 1) {
		holderStack.push_back(index);
		if (holders[index].theSet != NULL) {
			assert(find(claStack.begin(), claStack.end(), holders[index].theSet) == claStack.end());
			//assert(holders[index].theSet->NStates()==4);
			claStack.push_back(holders[index].theSet);
		}
		holders[index].Reset();
	}
	else {
		holders[index].numAssigned--;
		//this is important!
		holders[index].tempReserved = false;
	}
}

void ClaManager::CheckClaHolders() {
#ifndef NDEBUG
	int used = 0;
	int reclaim2 = 0;
	for (int i = 0; i<numHolders; i++) {
		if (holders[i].theSet != NULL) {
			used++;
			if (holders[i].GetReclaimLevel() == 2) reclaim2++;
		}
	}
	assert(used == numClas - claStack.size());
#endif
}

void ClaManager::MakeAllHoldersDirty() {

	for (int i = 0; i<numHolders; i++) {
		if (holders[i].theSet != NULL) {
			claStack.push_back(holders[i].theSet);
			holders[i].theSet = NULL;
		}
	}
}

//Ported from BEAGLE branch
//This returns the cla number of the cla that this holder points to
//it isn't changing anything, including indexes
int ClaManager::GetClaIndexForBeagle(int index) const {
	assert(index > -1 && index < numHolders);
	//BEAGLEMERGE
	//assert(holders[index].theArray != NULL);
	//return holders[index].theArray->Index();
	assert(holders[index].theSet != NULL);
	int cindex = holders[index].theSet->Index();
	assert(cindex >= 0);
	return cindex;
}

void ClaManager::SetHolderDependencies(int index, int depIndex1, int pDepIndex1, int depIndex2, int pDepIndex2) {
	assert(index > -1 && index < numHolders);
	assert(depIndex1 >= -(numNodes + 2) && depIndex1 < numHolders);
	assert(depIndex2 >= -(numNodes + 2) && depIndex2 < numHolders);

	holders[index].holderDep1 = depIndex1;
	holders[index].transMatDep1 = pDepIndex1;
	holders[index].holderDep2 = depIndex2;
	holders[index].transMatDep2 = pDepIndex2;
}

//this is essentially the same as ClaimClaFillIfNecessary, but it returns the assigned cla object
CondLikeArray *ClaManager::GetClaFillIfNecessary(int index, int modnum /*=0*/) {
	assert(index > -1 && index < numHolders);
	CheckClaHolders();
	//assert(holders[index].theArray != NULL);
	//DEBUG - HACK not sure if this is dangerous or not, or what should happen with the dir argument
	assert(holders[index].numAssigned > 0);
	if (holders[index].theSet == NULL) {
		FillHolder(index, 1);
		if (debug_clas) {
			ReportClaTotals("get-fill", index);
		}
	}
	else {
		if (debug_clas) {
			ReportClaTotals("get-already-filled", index);
		}
	}
	return holders[index].theSet->theSets[modnum];
}

//This is essentially the same as GetCla, but doesn't return anything and also sets the tempReserved flag
//FillHolder is only used if the holder is currently empty
void ClaManager::ClaimClaFillIfNecessary(int index, int depLevel) {
	assert(index > -1 && index < numHolders);
	assert(holders[index].numAssigned > 0);
	if (debug_clas) {
		outman.DebugMessage("claim %d dep %d", index, depLevel);
		CheckClaHolders();
	}

	if (IsHolderDirty(index)) {
		FillHolder(index, depLevel);
		if (debug_clas) {
			ReportClaTotals("claimed-fill", index);
		}
	}
	else {
		//FillHolder will also set the dep level, but we need to set it even if FillHolder doesn't need to be called
		SetDepLevel(index, depLevel);
		if (debug_clas) {
			ReportClaTotals("claimed-already-filled", index);
		}
	}
	//mark this cla as one not to be reclaimed, since it has been "claimed"
	TempReserveCla(index);
}

bool ClaManager::IsHolderDirty(int index) const {
	//dirtyness is synonymous with a null cla pointer in the holder
	assert(index > -1 && index < numHolders);
	assert(holders[index].numAssigned > 0);
	return (holders[index].theSet == NULL);
}

void ClaManager::ReportClaTotals(const char *mess, int index) const {
	if (debug_clas) {
		int clean, tempres, res, assigned;
		clean = tempres = res = assigned = 0;
		for (int i = 0; i<numHolders; i++) {
			if (holders[i].theSet != NULL)
				clean++;
			if (IsHolderTempReserved(i))
				tempres++;
			if (IsHolderReserved(i))
				res++;
		}
		assigned = numHolders - holderStack.size();
		outman.DebugMessage("%s\t%d\tstacks:\t%d\t%d\tclean\t%d\ttempres\t%d\tres\t%d\tassigned\t%d\tfreeable\t%d", mess, index, claStack.size(), holderStack.size(), clean, tempres, res, assigned, freeableHolderQueue.size());
		//the numbers can be a little out of whack because this is called mid-operation, so
		//this will only check for serious problems
		assert(abs(clean + (int)claStack.size() - numClas) <= 2);
	}
}

void ClaManager::SetDepLevel(int index, int lvl) {
	holders[index].depLevel = lvl;
	if (debug_clas) {
		char s[50];
		sprintf(s, "setdep\tlvl\t%d", lvl);
		ReportClaTotals(s, index);
	}

}

//Tell the manager that it is ok to reclaim the cla in this holder, although it might
//avoid doing so based on the dep level
//SEE ABOVE for distinction between normal and temp reservations
void ClaManager::RemoveNormalReservation(int index) {
	assert(index > -1 && index < numHolders);
	//holders[index].tempReserved=false;
	holders[index].reserved = false;
	if (memLevel>1)
		holders[index].SetReclaimLevel(1);
	if (debug_clas) {
		ReportClaTotals("unNreserved", index);
	}
}

//Tell the manager that it is ok to reclaim the cla in this holder, although it might
//avoid doing so based on the dep level
//SEE ABOVE for distinction between normal and temp reservations
void ClaManager::RemoveTempReservation(int index) {
	assert(index > -1 && index < numHolders);
	holders[index].tempReserved = false;

	if (memLevel>1)
		holders[index].SetReclaimLevel(1);
	if (debug_clas) {
		ReportClaTotals("unTreserved", index);
	}
}

//this simply removes the cla from a holder.  It is equivalent to just
//dirtying it if only a single tree shares the holder.  The holder itself
//remains valid and properly assigned, but is: 
//-unreserved (this shouldn't be getting called with reserved holders anyway)
//
void ClaManager::EmptyHolder(int index) {
	assert(index > -1 && index < numHolders);

	if (holders[index].IsDirty()) {
		if (debug_clas) {
			ReportClaTotals("alreadyEmpty\0", index);
		}
		return;
	}

	_ReclaimC(index);
	holders[index].reserved = holders[index].tempReserved = false;

	holders[index].SetReclaimLevel(0);
	if (debug_clas) {
		ReportClaTotals("emptied\0", index);
	}
}

int ClaManager::_GetFreeH() {
	CheckClaHolders();
	if (holderStack.empty() == true)
		NewRecycleClas();
	assert(holderStack.size() > 0);

	int newIndex = holderStack[holderStack.size() - 1];
	holderStack.pop_back();
	IncrementHolder(newIndex);
	return newIndex;
}

CondLikeArraySet* ClaManager::_GetFreeC() {
	//if(claStack.empty() == true) RecycleClas();
	if (claStack.empty() == true)
		NewRecycleClas();

	CondLikeArraySet *arr = claStack[claStack.size() - 1];

	assert(arr != NULL);
	claStack.pop_back();
	if (numClas - (int)claStack.size() > maxUsed) maxUsed = numClas - (int)claStack.size();
	if (debug_clas) {
		ReportClaTotals("assignedC\0", arr->Index());
	}
	return arr;
}

void ClaManager::_ReclaimH(int index) {
	if (!holders[index].IsDirty())
		_ReclaimC(index);

	holders[index].Reset();
	holderStack.push_back(index);
	if (debug_clas) {
		char s[50];
		sprintf(s, "reclaimH\t%d", index);
		ReportClaTotals(s, index);
	}
}

void ClaManager::_ReclaimC(int index) {
	claStack.push_back(holders[index].theSet);
	holders[index].theSet = NULL;
	if (debug_clas) {
		char s[50];
		sprintf(s, "reclaimC\t%d", index);
		ReportClaTotals(s, index);
	}
}

//this is called by node saying "I'm changing such that this holder is no longer valid
//FOR ME.  Please take my cla if I had one, and tell me if I can continue to use the same 
//index (because I was the only one referring to it), or assign me a new unused holder index"
int ClaManager::SetHolderDirty(int index) {
	//there are only two options here:
	//1. Cla is being made dirty, and only one node points to it 
	//	->null the holder's cla pointer and return the same index
	//2. Cla is being made dirty, and multiple nodes point to it
	//	->remove this node from the holder (decrement) and assign a new one	
	assert(index < numHolders);
	if(index == -1) {return index;}
	assert(holders[index].numAssigned > 0);
	CheckClaHolders();

	int newIndex = -1;
	if (GetNumAssigned(index) == 1) {
		//the same holder can be used in the new context
		if (holders[index].theSet != NULL) {
			//but if it had been previously calculated, grab that cla
			EmptyHolder(index);
		}
		newIndex = index;
	}
	else {
		//ok, this holder is still valid for someone else.  I'll turn this one in
		//and please give me a new one
		newIndex = TradeInClaHolder(index);
	}
	if (debug_clas) {
		char s[50];
		sprintf(s, "emptied\t%d\treassigned", index);
		//DEBUG
		//outman.DebugMessageNoCR("%s ", s);
		//ReportClaTotals(s, newIndex);
	}
	return newIndex;
}

void ClaManager::OutputClaReport() {
	//this needs work yet
	assert(0);

	/*
	ofstream cla("clareport.log");
	for (int i = 1; i<numClas; i++) {
		cla << i << "\t" << assignedClaArray[i] << "\t" << allClas[i]->nodeNum << "\t" << allClas[i]->IsDirty() << "\n";
	}
	*/
}

void ClaManager::NewRecycleClas() {
	assert(curReclaimLevel != -1);

	if (debug_clas) {
		ReportClaTotals("recycling", -1);
	}

	int numReclaimed = 0;
	int targetReclaim = min(numNodes / 2.0, 20.0);
	for (list<int>::iterator it = freeableHolderQueue.begin(); it != freeableHolderQueue.end();) {
		//if(holders[*it].depLevel > curReclaimLevel || numReclaimed > numNodes / 2)
		if (numReclaimed >= targetReclaim)
			break;
		assert(!IsHolderTempReserved(*it) && !IsHolderReserved(*it));
		EmptyHolder(*it);
		numReclaimed++;
		it++;
		freeableHolderQueue.pop_front();
	}

	if (numReclaimed >= targetReclaim) {
		ReportClaTotals("recycledFromQueueRet ", numReclaimed);
		return;
	}
	else
		ReportClaTotals("recycledFromQueueCont. ", numReclaimed);

	vector<int> level2s;
	for (int i = 0; i<numHolders; i++) {
		if (!holders[i].IsDirty()) {
			//grab any level 0 or 1 deps
			if ((holders[i].depLevel < 3) && (IsHolderTempReserved(i) == false) && (IsHolderReserved(i) == false)) {
				if (holders[i].depLevel < 2) {
					EmptyHolder(i);
					numReclaimed++;
					if (numReclaimed >= targetReclaim) {
						return;
					}
				}
				else
					level2s.push_back(i);
			}
		}
	}

	ReportClaTotals("recycled>numNodesLvl<2 ", numReclaimed);
	if (numReclaimed >= targetReclaim) {
		return;
	}

	//return the second level deps
	for (vector<int>::iterator it = level2s.begin(); it != level2s.end(); it++) {
		EmptyHolder(*it);
		numReclaimed++;
	}

	ReportClaTotals("recycled>numNodesLvl2 ", numReclaimed);
	if (numReclaimed >= targetReclaim) {
		return;
	}

	for (int i = 0; i<numHolders; i++) {
		if (!holders[i].IsDirty()) {
			//just grab anything not reserved
			if ((IsHolderTempReserved(i) == false) && (IsHolderReserved(i) == false)) {
				EmptyHolder(i);
				numReclaimed++;
			}
		}
		if (numReclaimed >= targetReclaim) {
			ReportClaTotals("recycled>numNodesAnyLevel ", numReclaimed);
			return;
		}
	}
	if (debug_clas) {
		ReportClaTotals("recycled ", numReclaimed);
		OutputClaReport();
	}
	assert(numReclaimed > 0);
}