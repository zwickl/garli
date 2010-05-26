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

class ClaManager{
	int numNodes;//the number of nodes in each tree
	int numRates;
	int numClas;
	int numHolders;
	int maxUsed;
	int curReclaimLevel;
	CondLikeArray **allClas; //these are the actual arrays to be used in calculations, but will assigned to 
							 //nodes via a CondLikeArrayHolder.  There may be a limited number
							 
	CondLikeArrayHolder *holders; //there will be enough of these such that every node and direction could
								  //have a unique one, although many will generally be shared

	vector<CondLikeArray *> claStack;
	//this is for use in beagle mode, when the objects being dealt with are indeces rather than pointers to CondLikeArrays
	vector<int> claIndexStack;

	vector<int> holderStack;

	public:	
	bool debug_clas;
		
	ClaManager(int nnod, int nClas, int nHolders, int nchar, int nrates) : curReclaimLevel(-1), numNodes(nnod), numClas(nClas), numHolders(nHolders), numRates(nrates){
		maxUsed=0;
		allClas=new CondLikeArray*[numClas];
		claStack.reserve(numClas);
		for(int i=numClas-1;i>=0;i--){
			allClas[i]=new CondLikeArray;
			allClas[i]->Allocate(i, modSpec.nstates, nchar, numRates);
			claStack.push_back(allClas[i]);
			claIndexStack.push_back(i);
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


	int AssignFreeClaHolder();
	int TradeInClaHolder(int oldIndex);
	void FillHolder(int index, int dir); //sorry Mark
	void EmptyHolder(int index);

	int GetReclaimLevel(int index);
	int GetNumAssigned(int index) {return holders[index].numAssigned;}
	void ReserveCla(int index, bool temp=true);
	void ClearTempReservation(int index) {
		if(debug_clas)
			outman.DebugMessage("remv temp res %d", index);
		holders[index].tempReserved=false;
		}

	void RemoveNormalReservation(int index);
	void RemoveTempReservation(int index);
	bool IsClaReserved(int index) {return holders[index].reserved;}
	bool IsClaTempReserved(int index) {return holders[index].tempReserved;};
	void ReportClaTotals(const char *mess, int index) const;
	void RecycleClas();
	void NewRecycleClas();
	int GetClaNumber(int index);
	int CountClasInUse(int recLevel);
	CondLikeArray *GetClaFillIfNecessary(int index);
	void ClaimClaFillIfNecessary(int index, int depLevel);
	int GetClaIndexForBeagle(int index) const;
	const CondLikeArrayHolder *GetHolder(int index);	
	CondLikeArrayHolder *GetMutableHolder(int index);
	bool IsHolderDirty(int index);
	int SetHolderDirty(int index);
	void IncrementHolder(int index);
	void DecrementHolder(int index);
	void CheckClaHolders();
	void MakeAllHoldersDirty();
	void GetHolderUsageCountTotals(vector<int> &);
	void SetHolderDependencies(int index, int depIndex1, int pDepIndex1, int depIndex2, int pDepIndex2);
	void SetCurReclaimLevel(int lvl){curReclaimLevel = lvl;}
	void OutputClaReport();
	void SetDepLevel(int index, int lev);

	private:
	int _GetFreeH();
	CondLikeArray* _GetFreeC();
	};
	
	inline int ClaManager::_GetFreeH(){
		CheckClaHolders();
		assert(holderStack.size() > 0);

		int newIndex=holderStack[holderStack.size()-1];		
		holderStack.pop_back();
		IncrementHolder(newIndex);
		return newIndex;
		}

	inline CondLikeArray* ClaManager::_GetFreeC(){
		//if(claStack.empty() == true) RecycleClas();
		if(claStack.empty() == true)	
			NewRecycleClas();

		CondLikeArray *arr=claStack[claStack.size()-1];

		assert(arr != NULL);
		claStack.pop_back();
		if(numClas - (int)claStack.size() > maxUsed) maxUsed=numClas - (int)claStack.size();
		if(debug_clas){
			ReportClaTotals("assignedC\0", arr->Index());
			}
		return arr;
		}

	inline void ClaManager::SetDepLevel(int index, int lvl){
		holders[index].depLevel = lvl;
		if(debug_clas){
			char s[50];
			sprintf(s, "setdep\tlvl\t%d", lvl);
			ReportClaTotals(s, index);
			}
		}

	inline int ClaManager::GetReclaimLevel(int index){
		assert(index > -1 && index < numHolders);
		if(holders[index].theArray == NULL) return -1;
		return holders[index].GetReclaimLevel();
		}

	//inline void ClaManager::CountClaTotals(int &clean, int &tempres, int &res, int &assigned){
	inline void ClaManager::ReportClaTotals(const char *mess, int index) const{
		if(debug_clas){
			int clean, tempres, res, assigned;
			clean = tempres = res = assigned = 0;
			for(int i=0;i<numHolders;i++){
				if(holders[i].theArray != NULL) clean++;
				if(holders[i].tempReserved ==true) tempres++;
				if(holders[i].reserved==true) res++;
				}		
			assigned = numHolders - holderStack.size();
			outman.DebugMessage("%s\t%d\tstacks:\t%d\t%d\tclean\t%d\ttempres\t%d\tres\t%d\tassigned\t%d", mess, index, claStack.size(), holderStack.size(), clean, tempres, res, assigned);
			}
		}
	
	inline int ClaManager::GetClaNumber(int index){
		//this is ugly, but should only be called for debugging
		assert(index > -1 && index < numHolders);
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

	//This returns the cla number of the cla that this holder points to
	//it isn't changing anything, including indexes
	inline int ClaManager::GetClaIndexForBeagle(int index) const{
		assert(index > -1 && index < numHolders);
		assert(holders[index].theArray != NULL);
		return holders[index].theArray->Index();
		}

	inline const CondLikeArrayHolder *ClaManager::GetHolder(int index){
		assert(index > -1 && index < numHolders);
		return &holders[index];
		}

	inline CondLikeArrayHolder *ClaManager::GetMutableHolder(int index){
		assert(index > -1 && index < numHolders);
		return &holders[index];
		}

	inline void ClaManager::SetHolderDependencies(int index, int depIndex1, int pDepIndex1, int depIndex2, int pDepIndex2){
		assert(index > -1 && index < numHolders);
		assert(depIndex1 >= -(numNodes + 2) && depIndex1 < numHolders);
		assert(depIndex2 >= -(numNodes + 2) && depIndex2 < numHolders);

		holders[index].holderDep1 = depIndex1;
		holders[index].transMatDep1 = pDepIndex1;
		holders[index].holderDep2 = depIndex2;
		holders[index].transMatDep2 = pDepIndex2;
		}

	inline bool ClaManager::IsHolderDirty(int index){
		//dirtyness is now synonymous with a null cla pointer in the holder
		assert(index > -1 && index < numHolders);
		assert(holders[index].numAssigned > 0);
		return (holders[index].theArray == NULL);	
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
	
	inline void ClaManager::OutputClaReport(){
		ofstream cla("clareport.log");
		cla << "hIndex\tclaIndex\tnumAss\treclaimLvl\tdepLvl\tRes\tTres\n";
		for(int i=0;i< numHolders;i++){
			cla << i << "\t" << (holders[i].theArray == NULL ? -1 : holders[i].theArray->Index()) << "\t" << holders[i].numAssigned << "\t";
			cla << holders[i].reclaimLevel << "\t" << holders[i].depLevel << "\t";
			cla << (holders[i].reserved ? "1" : "0") << "\t" << (holders[i].tempReserved ? "1" : "0")<< "\n";
			}
		cla.close();
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