// GARLI version 0.96b8 source code
// Copyright 2005-2008 Derrick J. Zwickl
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

#ifdef UNIX
	#include <sys/mman.h>
#endif

#undef CLA_DEBUG

class ClaSpecifier{
	public:
	int claIndex; //this is just the number of the corresponding cla - there is a 1 to 1 correspondence
	int modelIndex;
	int dataIndex;
	ClaSpecifier(int c, int m, int d):claIndex(c), modelIndex(m), dataIndex(d){};
	};

extern vector<ClaSpecifier> claSpecs;

class ClaManager{
	int numNodes;//the number of nodes in each tree
	int numRates;
	int numClas;
	int numHolders;
	int maxUsed;
	CondLikeArraySet **allClas; //these are the actual sets of arrays to be used in calculations, but will assigned to 
							 //nodes via a CondLikeArrayHolder.  There may be a limited number						 

	CondLikeArrayHolder *holders; //there will be enough of these such that every node and direction could
								  //have a unique one, although many will generally be shared
								  
	vector<CondLikeArraySet *> claStack;
	vector<int> holderStack;
	
	public:	
	//PARTITION	
	//ClaManager(int nnod, int nClas, int nHolders, int nchar, int nrates) : numNodes(nnod), numClas(nClas), numHolders(nHolders), numRates(nrates){
/*	ClaManager(int nnod, int numClas, int nHolders, const ModelPartition *mods, const DataPartition *data) : numNodes(nnod), numHolders(nHolders){
		maxUsed=0;
		allClas=new CondLikeArraySet*[numClas];
		claStack.reserve(numClas);
		for(int i=numClas-1;i>=0;i--){
			allClas[i]=new CondLikeArraySet;
			//for(vector<Model *>::iterator modit = mods->models.begin();modit != mods->models.end();modit++){
			for(int m = 0;m < mods->NumModels();m++){
				const Model *thisMod = mods->GetModel(m);
				CondLikeArray *thisCLA = new CondLikeArray(data->GetSubset(thisMod->dataIndex)->NChar(), thisMod->NStates(), thisMod->NRateCats(), thisMod->dataIndex);
				allClas[i]->AddCLA(thisCLA);
				allClas[i]->Allocate();
				}
			claStack.push_back(allClas[i]);
			}
		holders = new CondLikeArrayHolder[numHolders];
		holderStack.reserve(numHolders);
		for(int i=numHolders-1;i>=0;i--)
			holderStack.push_back(i);
		}
*/
	ClaManager(int nnod, int nClas, int nHolders, const ModelPartition *mods, const DataPartition *data) : numNodes(nnod), numClas(nClas), numHolders(nHolders){
		maxUsed=0;
		allClas=new CondLikeArraySet*[numClas];
		claStack.reserve(numClas);
		for(int i=numClas-1;i>=0;i--){
			allClas[i]=new CondLikeArraySet;
			//for(vector<Model *>::iterator modit = mods->models.begin();modit != mods->models.end();modit++){
			//for(int m = 0;m < mods->NumModels();m++){
			for(vector<ClaSpecifier>::iterator specs = claSpecs.begin();specs != claSpecs.end();specs++){
				const Model *thisMod = mods->GetModel((*specs).modelIndex);
				CondLikeArray *thisCLA = new CondLikeArray(data->GetSubset((*specs).dataIndex)->NChar(), thisMod->NStates(), thisMod->NRateCats());
				allClas[i]->AddCLA(thisCLA);
				}
			allClas[i]->Allocate();
			claStack.push_back(allClas[i]);
			}
		holders = new CondLikeArrayHolder[numHolders];
		holderStack.reserve(numHolders);
		for(int i=numHolders-1;i>=0;i--)
			holderStack.push_back(i);
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


	int AssignClaHolder();
	CondLikeArraySet* AssignFreeCla();
	void FillHolder(int index, int dir); //sorry Mark

	int GetReclaimLevel(int index);
	int GetNumAssigned(int index) {return holders[index].numAssigned;}
	void ReserveCla(int index, bool temp=true);
	void ClearTempReservation(int index) {holders[index].tempReserved=false;}
	void UnreserveCla(int index);
	bool IsClaReserved(int index) {return holders[index].reserved;}
	bool IsClaTempReserved(int index) {return holders[index].tempReserved;};
	void ReclaimSingleCla(int index);
	void CountClaTotals(int &clean, int &tempres, int &res);
	void RecycleClas();
	int GetClaNumber(int index);
	int CountClasInUse(int recLevel);
	CondLikeArraySet *GetCla(int index);	
	const CondLikeArrayHolder *GetHolder(int index);	
	bool IsDirty(int index);
	int SetDirty(int index);
	void IncrementCla(int index);
	void DecrementCla(int index);
	void CheckClaHolders();
	void MakeAllHoldersDirty();
	};
	
	inline int ClaManager::AssignClaHolder(){
		assert(holderStack.size() > 0);
		int index=holderStack[holderStack.size()-1];
		IncrementCla(index);
		holderStack.pop_back();
		return index;
		}
	
	inline void ClaManager::FillHolder(int index, int dir){
		holders[index].theSet = AssignFreeCla();
		holders[index].reclaimLevel=dir;
		}

	inline int ClaManager::GetReclaimLevel(int index){
		if(holders[index].theSet == NULL) return -1;
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
		if(holders[index].theSet==NULL) return;
		claStack.push_back(holders[index].theSet);
		holders[index].SetReclaimLevel(0);
		holders[index].theSet=NULL;				
		}

	inline void ClaManager::CountClaTotals(int &clean, int &tempres, int &res){
		for(int i=0;i<numHolders;i++){
			if(holders[i].theSet != NULL) clean++;
			if(holders[i].tempReserved ==true) tempres++;
			if(holders[i].reserved==true) res++;
			}		
		}
	
	inline int ClaManager::GetClaNumber(int index){
		//this is ugly, but should only be called for debugging
		if(holders[index].theSet == NULL) return -1;
		for(int i=0;i<numClas;i++)
			if(holders[index].theSet == allClas[i]) return i;
		assert(0);
		return -1;
		}

	inline int ClaManager::CountClasInUse(int recLevel){
		int num=0;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theSet != NULL)
				if(holders[i].GetReclaimLevel() == recLevel) num++;
			}
		return num;
		}

	inline CondLikeArraySet *ClaManager::GetCla(int index){
		assert(holders[index].theSet != NULL);
		return holders[index].theSet;
		}
	
	inline const CondLikeArrayHolder *ClaManager::GetHolder(int index){
		return &holders[index];
		}
	
	inline bool ClaManager::IsDirty(int index){
		//dirtyness is now synonymous with a null cla pointer in the holder
		assert(index > -1);
		return (holders[index].theSet == NULL);	
		}

	inline int ClaManager::SetDirty(int index){
		//there are only two options here:
		//1. Cla is being made dirty, and only node node points to it 
		//	->null the holder's cla pointer and return the same index
		//2. Cla is being made dirty, and multiple nodes point to it
		//	->remove this node from the holder (decrement) and assign a new one	
	
		assert(index != -1);

		if(holders[index].numAssigned==1){
			if(holders[index].theSet != NULL){
				holders[index].SetReclaimLevel(0);
				claStack.push_back(holders[index].theSet);
				holders[index].theSet=NULL;
				}
			}
		else{
			DecrementCla(index);
			index=holderStack[holderStack.size()-1];
			holderStack.pop_back();
			IncrementCla(index);
			}
		return index;
		}

	inline void ClaManager::IncrementCla(int index){
		holders[index].numAssigned++;
		}

	inline void ClaManager::DecrementCla(int index){
		assert(index != -1);
		if(holders[index].numAssigned==1){
			holderStack.push_back(index);
			if(holders[index].theSet != NULL){
				assert(find(claStack.begin(), claStack.end(), holders[index].theSet) == claStack.end());
				//assert(holders[index].theSet->NStates()==4);
				claStack.push_back(holders[index].theSet);
				}
			holders[index].Reset();
			}
		else{
			holders[index].numAssigned--; 
			//this is important!
			holders[index].tempReserved=false;
			}
		}

	inline void ClaManager::CheckClaHolders(){
		int used=0;
		int reclaim2=0;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theSet != NULL){
				used++;
				if(holders[i].GetReclaimLevel() == 2) reclaim2++;
				}
			}
		assert(used == numClas - claStack.size());
		}
	
	inline void ClaManager::MakeAllHoldersDirty(){
		
		for(int i=0;i<numHolders;i++){
			if(holders[i].theSet != NULL){
				claStack.push_back(holders[i].theSet);
				holders[i].theSet=NULL;
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
*/		
//	int NumCopies(){ return numCopies;}
//	int NumNodes(){return numNodes;}
/*	
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
