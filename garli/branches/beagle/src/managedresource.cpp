

#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

#include "defs.h"
#include "model.h"
#include "errorexception.h"
#include "outputman.h"
#include "managedresource.h"

//this will assign an unused holder index to a node that did not previously have one (no argument)
int ResourceManager::AssignFreeHolder(){
	int newIndex = _GetFreeHIndex();
	if(debug_resources){
		ReportHolderTotals("assignedH\0", newIndex);
		}
	return newIndex;
	}

//this is essentially identical to AssignNewClaHolder, but takes an old index
//back. it is written separately for clarity 
int ResourceManager::TradeInHolder(int oldIndex){
	int newIndex = _GetFreeHIndex();
	DecrementHolder(oldIndex);

	if(debug_resources){
		ReportHolderTotals("tradein\0", newIndex);
		}
	return newIndex;
	}

void ResourceManager::FillHolder(int index, int dir){
	assert(index > -1 && index < numHolders);
	CheckHolders();
	assert(holders[index].numAssigned > 0);

	holders[index].theResource = _GetFreeR();
	//holders[index].SetReclaimLevel(dir);
	if(debug_resources){
		char s[100];
		sprintf(s, "filled\tdir\t%d\tind", dir);
		ReportHolderTotals(s, index);
		}
	}

//tell the manager that another tree (node) now points to this holder.  
void ResourceManager::IncrementHolder(int index){
	assert(index > -1 && index < numHolders);
	CheckHolders();
	if(debug_resources)
		outman.DebugMessage("incr\t%d\t%d\t->\t%d", index, holders[index].numAssigned, holders[index].numAssigned+1);
	holders[index].numAssigned++;
	}

//tell the manager that one fewer tree (node) now points to this holder. If that brings the count to zero, then
//reclaim the holder, as well as any cla in it
void ResourceManager::DecrementHolder(int index){
	assert(index > -1 && index < numHolders);
	CheckHolders();
	if(debug_resources)
		outman.DebugMessage("decr\t%d\t%d\t->\t%d", index, holders[index].numAssigned, holders[index].numAssigned-1);

	if(holders[index].numAssigned==1){
		assert(find(holderIndexStack.begin(), holderIndexStack.end(), index) == holderIndexStack.end());
		holderIndexStack.push_back(index);
		if(debug_resources){
			ReportHolderTotals("reclaimH", index);
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

//this simply removes the cla from a holder.  It is equivalent to just
//dirtying it if only a single tree shares the holder.  Nothing happens 
//to the holder itself, which remains valid and properly assigned
void ResourceManager::EmptyHolder(int index){
	assert(index > -1 && index < numHolders);

	holders[index].reserved = holders[index].tempReserved = false;
	if(holders[index].theResource==NULL){
		if(debug_resources){
			ReportHolderTotals("alreadyEmpty\0", index);
			}
		return;
		}
	resourceStack.push_back(holders[index].theResource);
	//TODO
	//holders[index].SetReclaimLevel(0);
	holders[index].theResource=NULL;
	if(debug_resources){
		ReportHolderTotals("emptied\0", index);
		}
	}

//this marks this holder as one that should not have its cla reclaimed during recycling
//SEE BELOW for distinction between normal and temp reservations
void ResourceManager::ReserveResource(int index, bool temp/*=true*/){
	assert(index > -1 && index < numHolders);
	if(temp==true) 
		holders[index].tempReserved=true;
	else 
		holders[index].reserved=true;
	if(debug_resources){
		if(temp)
			ReportHolderTotals("tempreserved", index);
		else
			ReportHolderTotals("reserved", index);
		}
	}

//this is called by node saying "I'm changing such that this holder is no longer valid
//FOR ME.  Please take my cla if I had one, and tell me if I can continue to use the same 
//index (because I was the only one referring to it), or assign me a new unused holder index"
int ResourceManager::SetHolderDirty(int index){
	//there are only two options here:
	//1. Cla is being made dirty, and only one node points to it 
	//	->null the holder's cla pointer and return the same index
	//2. Cla is being made dirty, and multiple nodes point to it
	//	->remove this node from the holder (decrement) and assign a new one	
	CheckHolders();
	assert(index > -1 && index < numHolders);
	assert(holders[index].numAssigned > 0);
	
	int newIndex = -1;
	if(holders[index].numAssigned == 1){
		//the same holder can be used in the new context
		if(holders[index].theResource != NULL){
			//but if it had been previously calculated, grab that cla
			EmptyHolder(index);
			}
		newIndex = index;
		}
	else{
		//ok, this holder is still valid for someone else.  I'll turn this one in
		//and please give me a new one
		newIndex = TradeInHolder(index);
		}
	if(debug_resources){
		char s[50];
		sprintf(s, "emptied\t%d\treassigned", index);
		ReportHolderTotals(s, newIndex);
		}
	return newIndex;
	}

//Tell the manager that it is ok to reclaim the cla in this holder, although it might
//avoid doing so based on the dep level
//TEMP RESERVATIONS are nodes that are necessary for ongoing operations - destination arrays
//and immediate dependencies
void ResourceManager::RemoveTempReservation(int index){
	assert(index > -1 && index < numHolders);
	holders[index].tempReserved=false;
	//TODO
//	if(memLevel>1)
//		holders[index].SetReclaimLevel(1);
	if(debug_resources){
		ReportHolderTotals("unTreserved", index);
		}
	}

//Tell the manager that it is ok to reclaim the cla in this holder, although it might
//avoid doing so based on the dep level
//NORMAL RESERVATIONS are nodes that are reserved at a higher level, for example because
//they belong to the bestIndiv
void ResourceManager::RemoveNormalReservation(int index){
	assert(index > -1 && index < numHolders);
	//holders[index].tempReserved=false;
	holders[index].reserved=false;
	//TODO
//	if(memLevel>1)
//		holders[index].SetReclaimLevel(1);
	if(debug_resources){
		ReportHolderTotals("unNreserved", index);
		}
	}

//this is essentially the same as ClaimClaFillIfNecessary, but it returns the assigned cla object
Resource *ResourceManager::GetHolderFillIfNecessary(int index){
	assert(index > -1 && index < numHolders);
	CheckHolders();
	//assert(holders[index].theArray != NULL);
	//DEBUG - HACK not sure if this is dangerous or not, or what should happen with the dir argument
	assert(holders[index].numAssigned > 0);
	if(holders[index].theResource == NULL){
		FillHolder(index, 1);
		if(debug_resources){
			ReportHolderTotals("get-fill", index);
			}
		}
	else{
		if(debug_resources){
			ReportHolderTotals("get-already-filled", index);
			}
		}
	return holders[index].theResource;
	}

//This is essentially the same as GetCla, but doesn't return anything and also sets the reserved flag
//FillHolder is only used if the holder is currently empty
void ResourceManager::ClaimHolderFillIfNecessary(int index, int depLevel){
	assert(index > -1 && index < numHolders);
	outman.DebugMessage("claim %d dep %d", index, depLevel);
	CheckHolders();
	assert(holders[index].numAssigned > 0);
	//FillHolder will also set the dep level, but we need to set it even if FillHolder doesn't need to be called
	holders[index].depLevel = depLevel;
	if(holders[index].theResource == NULL){
		FillHolder(index, depLevel);
		if(debug_resources){
			ReportHolderTotals("claimed-fill", index);
			}
		}
	else{
		if(debug_resources){
			ReportHolderTotals("claimed-already-filled", index);
			}
		}
	//mark this cla as one not to be reclaimed
	ReserveResource(index, true);
	}

void ResourceManager::RecycleResources(){
	OutputResourceReport();
	assert(curReclaimLevel != -1);

	if(debug_resources){
		ReportHolderTotals("recycling", -1);
		}
	int numReclaimed=0;
	for(int i=0;i<numHolders;i++){
		if(holders[i].theResource != NULL){
			if(holders[i].depLevel < curReclaimLevel && holders[i].tempReserved == false && holders[i].reserved == false){
				EmptyHolder(i);
				numReclaimed++;
				}
			}
		if(memLevel<2){
			if(numReclaimed>50){
				if(debug_resources){
					ReportHolderTotals("recycled1", numReclaimed);
					}
				return;
				}
			}
		}
	if(debug_resources){
		ReportHolderTotals("recycled2", numReclaimed);
		}
	OutputResourceReport();
	assert(numReclaimed > 0);
	}

//for debugging
void ResourceManager::ReportHolderTotals(const char *mess, int index) const{
	if(debug_resources){
		int clean, tempres, res, assigned;
		clean = tempres = res = assigned = 0;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theResource != NULL)
				clean++;
			if(holders[i].tempReserved == true) 
				tempres++;
			if(holders[i].reserved == true) 
				res++;
			}		
		assigned = numHolders - holderIndexStack.size();
		outman.DebugMessage("%s\t%d\tstacks:\t%d\t%d\tclean\t%d\ttempres\t%d\tres\t%d\tassigned\t%d", mess, index, resourceStack.size(), holderIndexStack.size(), clean, tempres, res, assigned);
		}
	}

void ResourceManager::CheckHolders() const{
	if(debug_resources){
		int used=0;
		int reclaim2=0;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theResource != NULL){
				used++;
				if(holders[i].GetReclaimLevel() == 2) reclaim2++;
				}
			}
		assert(used == numResources - resourceStack.size());
		for(int i=0;i<numHolders;i++){
			if(find(holderIndexStack.begin(),holderIndexStack.end(), i) == holderIndexStack.end()){
				assert(holders[i].numAssigned > 0);
				}
			else{
				assert(holders[i].numAssigned == 0);
				assert(holders[i].theResource == NULL);
				}
			}
		}
	}

void ResourceManager::OutputResourceReport() const{
	ofstream cla("resourceReport.log");
	cla << "hIndex\tclaIndex\tnumAss\treclaimLvl\tdepLvl\tRes\tTres\n";
	for(int i=0;i< numHolders;i++){
		cla << i << "\t" << (holders[i].IsDirty() ? -1 : holders[i].theResource->Index()) << "\t" << holders[i].numAssigned << "\t";
		cla << holders[i].reclaimLevel << "\t" << holders[i].depLevel << "\t";
		cla << (holders[i].reserved ? "1" : "0") << "\t" << (holders[i].tempReserved ? "1" : "0")<< "\n";
		}
	cla.close();
	}

void CondLikeArray2::Allocate( int ind, int nk, int ns, int nr ){
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

CondLikeArray2::~CondLikeArray2(){
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

void CondLikeArrayHolder2::SetHolderDependencies(int depIndex1, int pDepIndex1, int depIndex2, int pDepIndex2){
	holderDep1 = depIndex1;
	transMatDep1 = pDepIndex1;
	holderDep2 = depIndex2;
	transMatDep2 = pDepIndex2;
	}

void ClaManager2::SetHolderDependencies(int index, int depIndex1, int pDepIndex1, int depIndex2, int pDepIndex2){
	assert(index > -1 && index < numHolders);
	assert(depIndex1 >= -(numNodes + 2) && depIndex1 < numHolders);
	assert(depIndex2 >= -(numNodes + 2) && depIndex2 < numHolders);

	CondLikeArrayHolder2 *hold = (CondLikeArrayHolder2*) &holders[index];
	hold->SetHolderDependencies(depIndex1, pDepIndex1, depIndex2, pDepIndex2);
/*
	holders[index].holderDep1 = depIndex1;
	holders[index].transMatDep1 = pDepIndex1;
	holders[index].holderDep2 = depIndex2;
	holders[index].transMatDep2 = pDepIndex2;
*/	}


