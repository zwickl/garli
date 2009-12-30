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

#ifndef CLA_MANAGER
#define CLA_MANAGER

#include <vector>
#include <algorithm>
#include <cassert>
#include "condlike.h"
#include "model.h"
using namespace std;

extern int memLevel;
extern ModelSpecification modSpec;

#ifdef UNIX
	#include <sys/mman.h>
#endif

#undef CLA_DEBUG

class ClaManager{
	int numNodes;//the number of nodes in each tree
	int numRates;
	int numClas;
	int numHolders;
	int maxUsed;
	CondLikeArray **allClas; //these are the actual arrays to be used in calculations, but will assigned to 
							 //nodes via a CondLikeArrayHolder.  There may be a limited number
							 
	CondLikeArrayHolder *holders; //there will be enough of these such that every node and direction could
								  //have a unique one, although many will generally be shared
								  
	vector<CondLikeArray *> claStack;
	vector<int> holderStack;

	//DEBUG
	bool debug_clas;
	
	public:	
	ClaManager(int nnod, int nClas, int nHolders, int nchar, int nrates) : numNodes(nnod), numClas(nClas), numHolders(nHolders), numRates(nrates){
		maxUsed=0;
		allClas=new CondLikeArray*[numClas];
		claStack.reserve(numClas);
		for(int i=numClas-1;i>=0;i--){
			allClas[i]=new CondLikeArray;
			allClas[i]->Allocate(modSpec.nstates, nchar, numRates);
			claStack.push_back(allClas[i]);
			}
		holders = new CondLikeArrayHolder[numHolders];
		holderStack.reserve(numHolders);
		for(int i=numHolders-1;i>=0;i--)
			holderStack.push_back(i);
#ifdef CLA_DEBUG
		debug_clas = true;
#else
		debug_clas = false;
#endif

		}

	~ClaManager(){
		if(allClas!=NULL){
			for(int i=0;i<numClas;i++){
				delete allClas[i];
				}
			delete []allClas;
			}
		delete []holders;
		}
	
	int NumClas() {return numClas;}
	int MaxUsedClas() {return maxUsed;}
	int NumFreeClas() {return (int) claStack.size();}
	int NumFreeHolders() {return (int) holderStack.size();}
	int NumHolders() {return numHolders;}


	int AssignClaHolder();
	CondLikeArray* AssignFreeCla();
	void FillHolder(int index, int dir); //sorry Mark

	int GetReclaimLevel(int index);
	int GetNumAssigned(int index) {return holders[index].numAssigned;}
	void ReserveCla(int index, bool temp=true);
	void ClearTempReservation(int index) {holders[index].tempReserved=false;}
	void UnreserveCla(int index);
	bool IsClaReserved(int index) {return holders[index].reserved;}
	bool IsClaTempReserved(int index) {return holders[index].tempReserved;};
	void ReclaimSingleCla(int index);
	void CountClaTotals(int &clean, int &tempres, int &res, int &assigned);
	void RecycleClas();
	int GetClaNumber(int index);
	int CountClasInUse(int recLevel);
	CondLikeArray *GetCla(int index);	
	const CondLikeArrayHolder *GetHolder(int index);	
	CondLikeArrayHolder *GetMutableHolder(int index);
	bool IsDirty(int index);
	int SetDirty(int index);
	void IncrementCla(int index);
	void DecrementCla(int index);
	void CheckClaHolders();
	void MakeAllHoldersDirty();
	void GetHolderUsageCountTotals(vector<int> &);
	void SetHolderDependencies(int index, int depIndex1, int pDepIndex1, int depIndex2, int pDepIndex2);
	};
	
	inline int ClaManager::AssignClaHolder(){
		CheckClaHolders();
		assert(holderStack.size() > 0);
		int index=holderStack[holderStack.size()-1];
		assert(holders[index].numAssigned == 0);
		IncrementCla(index);
		holderStack.pop_back();
		if(debug_clas)
			outman.DebugMessage("assign\t%d\tsizes:\t%d\t%d", index, claStack.size(), holderStack.size());
		return index;
		}
	
	inline void ClaManager::FillHolder(int index, int dir){
		CheckClaHolders();
		if(debug_clas)
			outman.DebugMessage("fill holder\t%d\tsizes:\t%d\t%d", index, claStack.size(), holderStack.size());
		assert(holders[index].numAssigned > 0);
		holders[index].theArray = AssignFreeCla();
		holders[index].reclaimLevel=dir;
		}

	inline int ClaManager::GetReclaimLevel(int index){
		if(holders[index].theArray == NULL) return -1;
		return holders[index].GetReclaimLevel();
		}

	inline void ClaManager::ReserveCla(int index, bool temp/*=true*/){
		if(temp==true) holders[index].tempReserved=true;
		else holders[index].reserved=true;
		}

	inline void ClaManager::UnreserveCla(int index){
//		holders[index].tempReserved=false;
		holders[index].reserved=false;
		if(memLevel>1)
			holders[index].SetReclaimLevel(1);
		}

	inline void ClaManager::ReclaimSingleCla(int index){
		//this simply removes the cla from a holder.  It is equivalent to just
		//dirtying it if only a single tree shares the holder
		if(holders[index].theArray==NULL) return;
		claStack.push_back(holders[index].theArray);
		holders[index].SetReclaimLevel(0);
		holders[index].theArray=NULL;				
		}

	inline void ClaManager::CountClaTotals(int &clean, int &tempres, int &res, int &assigned){
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL) clean++;
			if(holders[i].tempReserved ==true) tempres++;
			if(holders[i].reserved==true) res++;
			}		
		assigned = this->numHolders - holderStack.size();
		}
	
	inline int ClaManager::GetClaNumber(int index){
		//this is ugly, but should only be called for debugging
		if(holders[index].theArray == NULL) return -1;
		for(int i=0;i<numClas;i++)
			if(holders[index].theArray == allClas[i]) return i;
		assert(0);
		return -1;
		}

	inline int ClaManager::CountClasInUse(int recLevel){
		int num=0;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL)
				if(holders[i].GetReclaimLevel() == recLevel) num++;
			}
		return num;
		}

	inline CondLikeArray *ClaManager::GetCla(int index){
		CheckClaHolders();
		//assert(holders[index].theArray != NULL);
		//DEBUG - HACK not sure if this is dangerous or not, or what should happen with the dir argument
		assert(holders[index].numAssigned > 0);
		if(holders[index].theArray == NULL){
			FillHolder(index, 1);
//			outman.DebugMessage("fill and get cla from holder\t%d", index);
			}
		else{
			if(debug_clas)
				outman.DebugMessage("get existing cla from holder\t%d\tsizes:\t%d\t%d", index, claStack.size(), holderStack.size());
			}
		return holders[index].theArray;
		}
	
	inline const CondLikeArrayHolder *ClaManager::GetHolder(int index){
		return &holders[index];
		}

	inline CondLikeArrayHolder *ClaManager::GetMutableHolder(int index){
		return &holders[index];
		}

	inline void ClaManager::SetHolderDependencies(int index, int depIndex1, int pDepIndex1, int depIndex2, int pDepIndex2){
		holders[index].holderDep1 = depIndex1;
		holders[index].transMatDep1 = pDepIndex1;
		holders[index].holderDep2 = depIndex2;
		holders[index].transMatDep2 = pDepIndex2;
		}

	inline bool ClaManager::IsDirty(int index){
		//dirtyness is now synonymous with a null cla pointer in the holder
		assert(index > -1);
		assert(holders[index].numAssigned > 0);
		return (holders[index].theArray == NULL);	
		}

	inline int ClaManager::SetDirty(int index){
		if(debug_clas)
			outman.DebugMessage("dirty\t%d\tassigned\t%d", index, holders[index].numAssigned);
		CheckClaHolders();
		//there are only two options here:
		//1. Cla is being made dirty, and only node node points to it 
		//	->null the holder's cla pointer and return the same index
		//2. Cla is being made dirty, and multiple nodes point to it
		//	->remove this node from the holder (decrement) and assign a new one	
	
		assert(index > -1);
		assert(holders[index].numAssigned > 0);

		if(holders[index].numAssigned==1){
			if(holders[index].theArray != NULL){
				holders[index].SetReclaimLevel(0);
				claStack.push_back(holders[index].theArray);
				if(debug_clas)
					outman.DebugMessage("reclaim cla from holder\t%d\tsizes:\t%d\t%d", index, claStack.size(), holderStack.size());
				holders[index].theArray=NULL;
				}
			}
		else{
			DecrementCla(index);
			index=holderStack[holderStack.size()-1];
			assert(holders[index].numAssigned == 0);
			assert(holders[index].theArray == NULL);
			IncrementCla(index);
			holderStack.pop_back();
			}
		if(debug_clas)
			outman.DebugMessage("reassigned\t%d\tsizes:\t%d\t%d", index, claStack.size(), holderStack.size());
		assert(holders[index].numAssigned > 0);
		return index;
		}

	inline void ClaManager::IncrementCla(int index){
		//DEBUG
		CheckClaHolders();
		if(debug_clas)
			outman.DebugMessage("incr\t%d\t%d\t->\t%d", index, holders[index].numAssigned, holders[index].numAssigned+1);
		assert(index >= 0);
		holders[index].numAssigned++;
		}

	inline void ClaManager::DecrementCla(int index){

		CheckClaHolders();
		if(debug_clas)
			outman.DebugMessage("decr\t%d\t%d\t->\t%d", index, holders[index].numAssigned, holders[index].numAssigned-1);
		assert(index > -1);
		if(holders[index].numAssigned==1){
			assert(find(holderStack.begin(), holderStack.end(), index) == holderStack.end());
			holderStack.push_back(index);
			if(debug_clas)
				outman.DebugMessage("reclaim holder\t%d\tsizes:\t%d\t%d", index, claStack.size(), holderStack.size());
			if(holders[index].theArray != NULL){
				assert(find(claStack.begin(), claStack.end(), holders[index].theArray) == claStack.end());
				//assert(holders[index].theArray->NStates()==4);
				claStack.push_back(holders[index].theArray);
				}
			holders[index].Reset();
			}
		else{
			assert(holders[index].numAssigned != 0);
			holders[index].numAssigned--; 
			//this is important!
			holders[index].tempReserved=false;
			}
		}

	inline void ClaManager::CheckClaHolders(){
		if(debug_clas){
			int used=0;
			int reclaim2=0;
			for(int i=0;i<numHolders;i++){
				if(holders[i].theArray != NULL){
					used++;
					if(holders[i].GetReclaimLevel() == 2) reclaim2++;
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
	
	inline void ClaManager::GetHolderUsageCountTotals(vector<int> &counts){
		for(int h = 0;h < numHolders;h++){
			counts.push_back(holders[h].numAssigned);
			}
		}

	inline void ClaManager::MakeAllHoldersDirty(){
		
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL){
				claStack.push_back(holders[i].theArray);
				holders[i].theArray=NULL;
				}
			}
		}

/*
	void MarkReclaimable(int index, int val, bool observeCounts=true){
		assert(0);
		#ifdef CLA_DEBUG
		ofstream deb("cladebug.log", ios::app);
		deb << index << "\tmarked reclaimed\n";
		#endif
//		if(holders[index].theArray != NULL && (observeCounts==false || holders[index].nodes.size()==1)){
//		if(holders[index].theArray != NULL && (observeCounts==false || holders[index].numAssigned==1)){
//		if(holders[index].theArray != NULL && holders[index].reserved==false){
//			holders[index].SetReclaimLevel(val);

//		if(holders[index].theArray != NULL){
//			if(holders[index].reserved==false) holders[index].SetReclaimLevel(val);
			holders[index].tempReserved = false;
//			}
		}
*/	

/*		
	void CheckClaManager(int checktot){
		int tot=0;
		//verify that the total number of clas assigned for any particular node is less than the number of trees
		for(int i=1;i<size;i++){
			tot+=assignedClaArray[i];
			}
		}
*/
/*	
	void CheckAssignedNumber(int chk, int node, int index){
		assert(chk==assignedClaArray[node*numCopies+index]);
		}
*.		
//	int NumCopies(){ return numCopies;}
	int NumNodes(){return numNodes;}
*/ /*	
	void OutputAssignedClaArray(){
		static int count = 0;
		//ofstream out("claindeces.txt", ios::app);
		ofstream out("claindeces.txt");
		out << "\ncount " << count++ << endl;
		for(int i=0;i<numNodes;i++){
			for(int j=0;j<numCopies;j++){
				out << i << "\t" << j << "\t" << assignedClaArray[i*numCopies+j] << endl;
				}
			}
		out << endl;
		}
	void ClearAllClas(){
		for(int i=0;i<numNodes*numCopies;i++){
			allClas[i]->SetDirty(true);
			assignedClaArray[i]=0;
			}
		}
*/	
/*	void OutputClaReport(){
		ofstream cla("clareport.log");
		for(int i=1;i<numClas;i++){
			cla << i << "\t" << assignedClaArray[i] << "\t" << allClas[i]->nodeNum << "\t" << allClas[i]->IsDirty() << "\n";
			}
		}
*/
/*
	int NumAssigned(int index){
		return assignedClaArray[index];
		}	
*/
/*
	void OutputClaInfo(ostream &str, int index, int nd){
		str << nd << "\t" << index << "\t" << allClas[index]->nodeNum << "\t" << allClas[index]->IsDirty(nd) << "\t" << assignedClaArray[index] << "\t" << allClas[index]->GetReclaimLevel() << "\n";
		}
*/

#endif