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

	void ClaManager::FillHolder(int index, int dir){
		assert(index > -1 && index < numHolders);
		CheckClaHolders();
		assert(holders[index].numAssigned > 0);

		holders[index].theArray = _GetFreeC();
		holders[index].SetReclaimLevel(dir);
		if(debug_clas){
			char s[100];
			sprintf(s, "filled\tdir\t%d\tind", dir);
			ReportClaTotals(s, index);
			}
		}

	//this simply removes the cla from a holder.  It is equivalent to just
	//dirtying it if only a single tree shares the holder.  Nothing happens 
	//to the holder itself, which remains valid and properly assigned
	void ClaManager::EmptyHolder(int index){
		assert(index > -1 && index < numHolders);

		holders[index].reserved = holders[index].tempReserved = false;
		if(holders[index].theArray==NULL){
			if(debug_clas){
				ReportClaTotals("alreadyEmpty\0", index);
				}
			return;
			}
		claStack.push_back(holders[index].theArray);
		holders[index].SetReclaimLevel(0);
		holders[index].theArray=NULL;
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
		CheckClaHolders();
		assert(index > -1 && index < numHolders);
		assert(holders[index].numAssigned > 0);
		
		int newIndex = -1;
		if(holders[index].numAssigned == 1){
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
			ReportClaTotals(s, newIndex);
			}
		return newIndex;
		}

	//this marks this holder as one that should not have its cla reclaimed during recycling
	//SEE BELOW for distinction between normal and temp reservations
	void ClaManager::ReserveCla(int index, bool temp/*=true*/){
		assert(index > -1 && index < numHolders);
		if(temp==true) holders[index].tempReserved=true;
		else holders[index].reserved=true;
		if(debug_clas){
			if(temp)
				ReportClaTotals("tempreserved", index);
			else
				ReportClaTotals("reserved", index);
			}
		}

	//Tell the manager that it is ok to reclaim the cla in this holder, although it might
	//avoid doing so based on the dep level
	//TEMP RESERVATIONS are nodes that are necessary for ongoing operations - destination arrays
	//and immediate dependencies
	void ClaManager::RemoveTempReservation(int index){
		assert(index > -1 && index < numHolders);
		holders[index].tempReserved=false;
		if(memLevel>1)
			holders[index].SetReclaimLevel(1);
		if(debug_clas){
			ReportClaTotals("unTreserved", index);
			}
		}

	//Tell the manager that it is ok to reclaim the cla in this holder, although it might
	//avoid doing so based on the dep level
	//NORMAL RESERVATIONS are nodes that are reserved at a higher level, for example because
	//they belong to the bestIndiv
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

	//This is essentially the same as GetCla, but doesn't return anything and also sets the reserved flag
	//FillHolder is only used if the holder is currently empty
	void ClaManager::ClaimClaFillIfNecessary(int index, int depLevel){
		assert(index > -1 && index < numHolders);
		if(debug_clas)
			outman.DebugMessage("claim %d dep %d", index, depLevel);
		CheckClaHolders();
		assert(holders[index].numAssigned > 0);
		//FillHolder will also set the dep level, but we need to set it even if FillHolder doesn't need to be called
		holders[index].depLevel = depLevel;
		if(holders[index].theArray == NULL){
			FillHolder(index, depLevel);
			if(debug_clas){
				ReportClaTotals("claimed-fill", index);
				}
			}
		else{
			if(debug_clas){
				ReportClaTotals("claimed-already-filled", index);
				}
			}
		//mark this cla as one not to be reclaimed
		ReserveCla(index, true);
		}

	//tell the manager that another tree (node) now points to this holder.  
	void ClaManager::IncrementHolder(int index){
		assert(index > -1 && index < numHolders);
		holders[index].numAssigned++;
		if(debug_clas)
			outman.DebugMessage("incr\t%d\t->\t%d", index, holders[index].numAssigned);
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
			//holders[index].tempReserved=false;
			}
		}

	void ClaManager::RecycleClas(){
		int numReclaimed=0;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL){
				if(holders[i].GetReclaimLevel() == 2 && holders[i].tempReserved == false && holders[i].reserved == false){
					EmptyHolder(i);
					numReclaimed++;
					}
				}
			if(memLevel<2) if(numReclaimed>50) return;
			}
		if(numReclaimed>10) return;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL){
				if((holders[i].GetReclaimLevel() == 1 && holders[i].tempReserved == false && holders[i].reserved == false)){
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
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL){
				if(holders[i].depLevel < curReclaimLevel && holders[i].tempReserved == false && holders[i].reserved == false){
					EmptyHolder(i);
					numReclaimed++;
					}
				}
			if(numReclaimed > numNodes){
				ReportClaTotals("recycled>numNodes ", numReclaimed);
				return;
				}
			}
		if(debug_clas){
			ReportClaTotals("recycled ", numReclaimed);
			}
		if(debug_clas)
			OutputClaReport();
		assert(numReclaimed > 0);
		}