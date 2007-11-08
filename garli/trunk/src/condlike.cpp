// GARLI version 0.96b4 source code
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

#include "defs.h"
#include "condlike.h"
#include "clamanager.h"
#include "utility.h"

#undef ALIGN_CLAS
#define CLA_ALIGNMENT 32

CondLikeArray::~CondLikeArray(){
	//don't want to delete shared CL from nodes.  Should only
	//be called from Population level if CONDLIKE SHARED is defined
	if( arr ){
#ifndef ALIGN_CLAS
		delete []arr;
#else
		DeleteAlignedArray(arr);
#endif
		arr=NULL;
		 }
	if(underflow_mult!=NULL) delete []underflow_mult;
}

void CondLikeArray::Allocate( int nk, int ns, int nr /* = 1 */ ){
	if( arr ){
#ifndef ALIGN_CLAS
		delete []arr;
#else
		DeleteAlignedArray(arr);
#endif
		arr=NULL;
		}
	nrates = nr;
	nsites = ns;
	nstates = nk;
#ifndef ALIGN_CLAS
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


	CondLikeArray* ClaManager::AssignFreeCla(){
		#ifdef CLA_DEBUG
		ofstream deb("cladebug.log", ios::app);
		#endif

		if(claStack.empty() == true) RecycleClas();
		
		CondLikeArray *arr=claStack[claStack.size()-1];
		
		assert(arr != NULL);
		claStack.pop_back();
		if(numClas - (int)claStack.size() > maxUsed) maxUsed=numClas - (int)claStack.size();
		
		return arr;
		}

	void ClaManager::RecycleClas(){
		int numReclaimed=0;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL){
				if(holders[i].GetReclaimLevel() == 2 && holders[i].tempReserved == false && holders[i].reserved == false){
					assert(holders[i].theArray->NStates()==4);
					claStack.push_back(holders[i].theArray);
					holders[i].SetReclaimLevel(0);
					holders[i].theArray=NULL;
					numReclaimed++;
					}
				}
			if(memLevel<2) if(numReclaimed>50) return;
			}
		if(numReclaimed>10) return;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL){
				if((holders[i].GetReclaimLevel() == 1 && holders[i].tempReserved == false && holders[i].reserved == false)){
					assert(holders[i].theArray->NStates()==4);
					claStack.push_back(holders[i].theArray);
					holders[i].SetReclaimLevel(0);
					holders[i].theArray=NULL;
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