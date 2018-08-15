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

void CondLikeArray::Allocate( int nk, int ns, int nr /* = 1 */ ){
	if( arr ){
		delete []arr;
		arr=NULL;
		}
	nrates = nr;
	nsites = ns;
	nstates = nk;
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
