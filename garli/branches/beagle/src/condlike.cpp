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
//
//	NOTE: Portions of this source adapted from GAML source, written by Paul O. Lewis

#include "defs.h"
#include "condlike.h"
#include "clamanager.h"
#include "utility.h"
#include "calculationmanager.h"

#ifdef CUDA_GPU
#include "cudaman.h"
extern CudaManager *cudaman;
#endif

#undef ALIGN_CLAS
#define CLA_ALIGNMENT 32

CondLikeArray::~CondLikeArray(){
	//don't want to delete shared CL from nodes.  Should only
	//be called from Population level if CONDLIKE SHARED is defined
	if( arr ){
#ifndef ALIGN_CLAS
#ifdef CUDA_GPU
	if(cudaman->GetPinnedMemoryEnabled())
		FreePinnedMemory(arr);
	else
#endif
		delete []arr;


#else
		DeleteAlignedArray(arr);
#endif
		arr=NULL;
		 }
	if(underflow_mult!=NULL) 
		delete []underflow_mult;
	underflow_mult = NULL;
}

void CondLikeArray::Allocate( int ind, int nk, int ns, int nr /* = 1 */ ){
	if( arr ){
#ifndef ALIGN_CLAS
#ifdef CUDA_GPU
	if(cudaman->GetPinnedMemoryEnabled())
		FreePinnedMemory(arr);
	else
#endif
		delete []arr;

#else
		DeleteAlignedArray(arr);
#endif
		arr=NULL;
		}
	nrates = nr;
	nsites = ns;
	nstates = nk;
	index = ind;
#ifndef ALIGN_CLAS
#ifdef CUDA_GPU
if(cudaman->GetPinnedMemoryEnabled())
	AllocatePinnedMemory((void**)&arr, sizeof(FLOAT_TYPE)*nk*nr*ns);
else
#endif
	arr=new FLOAT_TYPE[nk*nr*ns];


#else
	arr = NewAlignedArray<FLOAT_TYPE>(nk*nr*ns, CLA_ALIGNMENT);
#endif
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

	const CondLikeArrayHolder *ClaManager::GetHolder(int index){
		assert(index > -1 && index < numHolders);
		return &holders[index];
		}

	CondLikeArrayHolder *ClaManager::GetMutableHolder(int index){
		assert(index > -1 && index < numHolders);
		return &holders[index];
		}

		//this will assign an unused holder index to a node that did not previously have one (no argument)
	int ClaManager::AssignFreeClaHolder(){
		int newIndex = _GetFreeH();
		if(debug_clas){
			ReportClaTotals("assignedH\0", newIndex);
			}
		return newIndex;
		}

	//this is essentially identical to AssignNewClaHolder, but takes an old index
	//back. it is written separately for clarity 
	int ClaManager::TradeInClaHolder(int oldIndex){
		int newIndex = _GetFreeH();
		DecrementHolder(oldIndex);

		if(debug_clas){
			ReportClaTotals("tradein\0", newIndex);
			}
		return newIndex;
		}

	void ClaManager::FillHolder(int index, int depLevel){
		assert(index > -1 && index < numHolders);
		CheckClaHolders();
		assert(holders[index].numAssigned > 0);

		holders[index].theArray = _GetFreeC();
		holders[index].SetDepLevel(depLevel);
		if(debug_clas){
			char s[100];
			sprintf(s, "filled\tdep\t%d\tind", depLevel);
			ReportClaTotals(s, index);
			}
		}

	//this simply removes the cla from a holder.  It is equivalent to just
	//dirtying it if only a single tree shares the holder.  The holder itself
	//remains valid and properly assigned, but is: 
	//-unreserved (this shouldn't be getting called with reserved holders anyway)
	//
	void ClaManager::EmptyHolder(int index){
		assert(index > -1 && index < numHolders);
	
		if(holders[index].IsDirty()){
			if(debug_clas){
				ReportClaTotals("alreadyEmpty\0", index);
				}
			return;
			}

		_ReclaimC(index);
		holders[index].reserved = holders[index].tempReserved = false;

		holders[index].SetReclaimLevel(0);
		if(debug_clas){
			ReportClaTotals("emptied\0", index);
			}
		}

	//this is called by node saying "I'm changing such that this holder is no longer valid
	//FOR ME.  Please take my cla if I had one, and tell me if I can continue to use the same 
	//index (because I was the only one referring to it), or assign me a new unused holder index"
	int ClaManager::SetHolderDirty(int index){
		//there are only two options here:
		//1. Cla is being made dirty, and only one node points to it 
		//	->null the holder's cla pointer and return the same index
		//2. Cla is being made dirty, and multiple nodes point to it
		//	->remove this node from the holder (decrement) and assign a new one	
		assert(index > -1 && index < numHolders);
		assert(holders[index].numAssigned > 0);
		CheckClaHolders();
		
		int newIndex = -1;
		if(NumAssigned(index) == 1){
			//the same holder can be used in the new context
			if(holders[index].theArray != NULL){
				//but if it had been previously calculated, grab that cla
				EmptyHolder(index);
				}
			newIndex = index;
			}
		else{
			//ok, this holder is still valid for someone else.  I'll turn this one in
			//and please give me a new one
			newIndex = TradeInClaHolder(index);
			}
		if(debug_clas){
			char s[50];
			sprintf(s, "emptied\t%d\treassigned", index);
			//DEBUG
			//outman.DebugMessageNoCR("%s ", s);
			//ReportClaTotals(s, newIndex);
			}
		return newIndex;
		}

	//this marks this holder as one that should not have its cla reclaimed during recycling
	//NORMAL RESERVATIONS are nodes that are reserved at a higher level, for example because
	//they belong to the bestIndiv
	void ClaManager::ReserveCla(int index){
		assert(index > -1 && index < numHolders);
		holders[index].reserved = true;
		if(debug_clas){
			ReportClaTotals("reserved", index);
			}
		}

	//this marks this holder as one that should not have its cla reclaimed during recycling
	//TEMP RESERVATIONS are nodes that are necessary for ongoing operations - destination arrays
	//and immediate dependencies
	void ClaManager::TempReserveCla(int index){
		assert(index > -1 && index < numHolders);
		holders[index].tempReserved = true;
		if(debug_clas){
			ReportClaTotals("tempreserved", index);
			}
		}

	//Tell the manager that it is ok to reclaim the cla in this holder, although it might
	//avoid doing so based on the dep level
	//SEE ABOVE for distinction between normal and temp reservations
	void ClaManager::RemoveNormalReservation(int index){
		assert(index > -1 && index < numHolders);
		//holders[index].tempReserved=false;
		holders[index].reserved=false;
		if(memLevel>1)
			holders[index].SetReclaimLevel(1);
		if(debug_clas){
			ReportClaTotals("unNreserved", index);
			}
		}

	//Tell the manager that it is ok to reclaim the cla in this holder, although it might
	//avoid doing so based on the dep level
	//SEE ABOVE for distinction between normal and temp reservations
	void ClaManager::RemoveTempReservation(int index){
		assert(index > -1 && index < numHolders);
		holders[index].tempReserved=false;
		
		if(memLevel>1)
			holders[index].SetReclaimLevel(1);
		if(debug_clas){
			ReportClaTotals("unTreserved", index);
			}
		}

	void ClaManager::MarkOperationAsFreeable(ClaOperation &op){
		AddToFreeableQueue(op.childClaIndex1);
		AddToFreeableQueue(op.childClaIndex2);
		}

	//this is essentially the same as ClaimClaFillIfNecessary, but it returns the assigned cla object
	CondLikeArray *ClaManager::GetClaFillIfNecessary(int index){
		assert(index > -1 && index < numHolders);
		CheckClaHolders();
		//assert(holders[index].theArray != NULL);
		//DEBUG - HACK not sure if this is dangerous or not, or what should happen with the dir argument
		assert(holders[index].numAssigned > 0);
		if(holders[index].theArray == NULL){
			FillHolder(index, 1);
			if(debug_clas){
				ReportClaTotals("get-fill", index);
				}
			}
		else{
			if(debug_clas){
				ReportClaTotals("get-already-filled", index);
				}
			}
		return holders[index].theArray;
		}

	//This is essentially the same as GetCla, but doesn't return anything and also sets the tempReserved flag
	//FillHolder is only used if the holder is currently empty
	void ClaManager::ClaimClaFillIfNecessary(int index, int depLevel){
		assert(index > -1 && index < numHolders);
		assert(holders[index].numAssigned > 0);
		if(debug_clas){
			outman.DebugMessage("claim %d dep %d", index, depLevel);
			CheckClaHolders();
			}
		
		if(IsHolderDirty(index)){
			FillHolder(index, depLevel);
			if(debug_clas){
				ReportClaTotals("claimed-fill", index);
				}
			}
		else{
			//FillHolder will also set the dep level, but we need to set it even if FillHolder doesn't need to be called
			SetDepLevel(index, depLevel);
			if(debug_clas){
				ReportClaTotals("claimed-already-filled", index);
				}
			}
		//mark this cla as one not to be reclaimed, since it has been "claimed"
		TempReserveCla(index);
		}

	//tell the manager that another tree (node) now points to this holder.  
	void ClaManager::IncrementHolder(int index){
		assert(index > -1 && index < numHolders);
		holders[index].numAssigned++;
		if(debug_clas){
			if(this->IsHolderReserved(index))
				outman.DebugMessage("incr(res)\t%d\t->\t%d", index, holders[index].numAssigned);
			else
				outman.DebugMessage("incr(noRes)\t%d\t->\t%d", index, holders[index].numAssigned);
			}
		CheckClaHolders();
		}

	//tell the manager that one fewer tree (node) now points to this holder. If that brings the count to zero, then
	//reclaim the holder, as well as any cla in it
	void ClaManager::DecrementHolder(int index){
		assert(index > -1 && index < numHolders);
		CheckClaHolders();
		if(debug_clas)
			outman.DebugMessage("decr\t%d\t%d\t->\t%d", index, holders[index].numAssigned, holders[index].numAssigned-1);

		if(holders[index].numAssigned==1){
			if(debug_clas){
				assert(find(holderStack.begin(), holderStack.end(), index) == holderStack.end());
				}
			holderStack.push_back(index);
			if(debug_clas){
				ReportClaTotals("reclaimH", index);
				}
			EmptyHolder(index);
			holders[index].Reset();
			}
		else{
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

	void ClaManager::RecycleClas(){
		int numReclaimed=0;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL){
				if(holders[i].DepLevel() == 2 && holders[i].tempReserved == false && holders[i].reserved == false){
					EmptyHolder(i);
					numReclaimed++;
					}
				}
			if(memLevel<2) if(numReclaimed>50) return;
			}
		if(numReclaimed>10) return;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL){
				if((holders[i].DepLevel() == 1 && holders[i].tempReserved == false && holders[i].reserved == false)){
					EmptyHolder(i);
					numReclaimed++;
					}
				}
			if(numReclaimed==20) return;
			}
		if(numReclaimed==0){
			throw(2);
			}
		assert(numReclaimed > 0);
		}

	void ClaManager::NewRecycleClas(){
		assert(curReclaimLevel != -1);

		if(debug_clas){
			ReportClaTotals("recycling", -1);
			}

		int numReclaimed=0;
		int targetReclaim = min(numNodes / 2.0, 20.0);
		for(list<int>::iterator it = freeableHolderQueue.begin(); it != freeableHolderQueue.end();){
			//if(holders[*it].depLevel > curReclaimLevel || numReclaimed > numNodes / 2)
			if(numReclaimed >= targetReclaim)
				break;
			assert(!IsHolderTempReserved(*it) && !IsHolderReserved(*it));
			EmptyHolder(*it);
			numReclaimed++;
			it++;
			freeableHolderQueue.pop_front();
			}

		if(numReclaimed >= targetReclaim){
			ReportClaTotals("recycledFromQueueRet ", numReclaimed);
			return;
			}
		else
			ReportClaTotals("recycledFromQueueCont. ", numReclaimed);

		vector<int> level2s;
		for(int i=0;i<numHolders;i++){
			if(!holders[i].IsDirty()){
				//grab any level 0 or 1 deps
				if((holders[i].depLevel < 3) && (IsHolderTempReserved(i) == false) && (IsHolderReserved(i) == false)){
					if(holders[i].depLevel < 2){
						EmptyHolder(i);
						numReclaimed++;
						if(numReclaimed >= targetReclaim){
							return;
							}
						}
					else 
						level2s.push_back(i);
					}
				}
			}

		ReportClaTotals("recycled>numNodesLvl<2 ", numReclaimed);
		if(numReclaimed >= targetReclaim){
			return;
			}

		//return the second level deps
		for(vector<int>::iterator it = level2s.begin(); it != level2s.end();it++){
			EmptyHolder(*it);
			numReclaimed++;
			}

		ReportClaTotals("recycled>numNodesLvl2 ", numReclaimed);
		if(numReclaimed >= targetReclaim){
			return;
			}

		for(int i=0;i<numHolders;i++){
			if(!holders[i].IsDirty()){
				//just grab anything not reserved
				if((IsHolderTempReserved(i) == false) && (IsHolderReserved(i) == false)){
					EmptyHolder(i);
					numReclaimed++;
					}
				}
			if(numReclaimed >= targetReclaim){
				ReportClaTotals("recycled>numNodesAnyLevel ", numReclaimed);
				return;
				}
			}
		if(debug_clas){
			ReportClaTotals("recycled ", numReclaimed);
			OutputClaReport();
			}
		assert(numReclaimed > 0);
		}


//DEBUGGING FUNCS

	int ClaManager::GetClaNumber(int index) const{
		//this is ugly, but should only be called for debugging
		assert(index > -1 && index < numHolders);
		if(holders[index].theArray == NULL) return -1;
		for(int i=0;i<numClas;i++)
			if(holders[index].theArray == allClas[i]) return i;
		assert(0);
		return -1;
		}

	int ClaManager::CountClasInUse(int recLevel) const{
		int num=0;
		for(int i=0;i<numHolders;i++){
			if(! IsHolderDirty(i))
				if(holders[i].DepLevel() == recLevel) 
					num++;
			}
		return num;
		}

	void ClaManager::ReportClaTotals(const char *mess, int index) const{
		if(debug_clas){
			int clean, tempres, res, assigned;
			clean = tempres = res = assigned = 0;
			for(int i=0;i<numHolders;i++){
				if(holders[i].theArray != NULL)
					clean++;
				if(IsHolderTempReserved(i)) 
					tempres++;
				if(IsHolderReserved(i))
					res++;
				}		
			assigned = numHolders - holderStack.size();
			outman.DebugMessage("%s\t%d\tstacks:\t%d\t%d\tclean\t%d\ttempres\t%d\tres\t%d\tassigned\t%d\tfreeable\t%d", mess, index, claStack.size(), holderStack.size(), clean, tempres, res, assigned, freeableHolderQueue.size());
			//the numbers can be a little out of whack because this is called mid-operation, so
			//this will only check for serious problems
			assert(abs(clean + (int) claStack.size() - numClas) <= 2);
			}
		}

	void ClaManager::CheckClaHolders() const{
		//DEBUG
		//if(debug_clas){
		if(0){
			int used=0;
			int reclaim2=0;
			for(int i=0;i<numHolders;i++){
				if(holders[i].theArray != NULL){
					used++;
					if(holders[i].DepLevel() == 2) reclaim2++;
					}
				}
			assert(used == numClas - claStack.size());
			for(int i=0;i<numHolders;i++){
				if(find(holderStack.begin(),holderStack.end(), i) == holderStack.end()){
					assert(holders[i].numAssigned > 0);
					}
				else{
					assert(holders[i].numAssigned == 0);
					assert(holders[i].theArray == NULL);
					}
				}
			}
		}
	
	void ClaManager::OutputClaReport() const{
		ofstream cla("clareport.log");
		cla << "hIndex\tclaIndex\tnumAss\treclaimLvl\tdepLvl\tRes\tTres\n";
		for(int i=0;i< numHolders;i++){
			cla << i << "\t" << (holders[i].theArray == NULL ? -1 : holders[i].theArray->Index()) << "\t" << holders[i].numAssigned << "\t";
			cla << holders[i].reclaimLevel << "\t" << holders[i].depLevel << "\t";
			cla << (holders[i].reserved ? "1" : "0") << "\t" << (holders[i].tempReserved ? "1" : "0")<< "\n";
			}
		cla.close();
		}
