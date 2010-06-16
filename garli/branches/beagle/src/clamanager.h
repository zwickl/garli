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

class ClaOperation;

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

	vector<CondLikeArray *> claStack; //these are the object being managed

	vector<int> holderStack; //these are holders that are unassigned

	list<int> freeableHolderQueue; //these are holders that are assigned (or at least were when they were put in here),
								   //have a valid cla in them and could be recycled.
								   //they wiil be added as a queue, with lowest level dependencies appearing first, etc.
								   //the curreclaim level indicates how far up the queue can be recycled.  If that isn't
								   //enough then brute force iteration over holders will be necessary

public:	
	bool debug_clas;
	ClaManager(){
		//shouldn't be calling the default constructor
		assert(0);
		}
	ClaManager(int nnod, int nClas, int nHolders, int nchar, int nrates) : curReclaimLevel(-1), numNodes(nnod), numClas(nClas), numHolders(nHolders), numRates(nrates){
		maxUsed=0;
		allClas=new CondLikeArray*[numClas];
		claStack.reserve(numClas);
		for(int i=numClas-1;i>=0;i--){
			allClas[i]=new CondLikeArray;
			allClas[i]->Allocate(i, modSpec.nstates, nchar, numRates);
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

private:
	//"primative" functions - the only ones that should be directly dealing with the stacks themselves
	int _GetFreeH();
	CondLikeArray* _GetFreeC();
	void _ReclaimH(int index);
	void _ReclaimC(int index);

public:
	//querying functions
	int NumClas() const {return numClas;}
	int MaxUsedClas() const {return maxUsed;}
	int NumFreeClas() const {return (int) claStack.size();}
	int NumFreeHolders() const {return (int) holderStack.size();}
	int NumHolders() const {return numHolders;}
	int NumAssigned(int index) const {return holders[index].numAssigned;}
	int ReclaimLevel(int index) const;
	int DepLevel(int index) const {return holders[index].depLevel;}
	bool IsHolderReserved(int index) const {return holders[index].reserved;}
	bool IsHolderTempReserved(int index) const {return holders[index].tempReserved;};
	bool IsHolderDirty(int index) const;
	
	//debugging funcs
	int CountClasInUse(int recLevel) const;
	void CheckClaHolders() const;
	void ReportClaTotals(const char *mess, int index) const;
	int GetClaNumber(int index) const;	void OutputClaReport() const;

	int AssignFreeClaHolder();
	int TradeInClaHolder(int oldIndex);
	void FillHolder(int index, int dir);
	void EmptyHolder(int index);

	void ReserveCla(int index);
	void TempReserveCla(int index);

	void RemoveNormalReservation(int index);
	void RemoveTempReservation(int index);

	void MarkOperationAsFreeable(ClaOperation &op);

	void RecycleClas();
	void NewRecycleClas();
	CondLikeArray *GetClaFillIfNecessary(int index);
	void ClaimClaFillIfNecessary(int index, int depLevel);
	int GetClaIndexForBeagle(int index) const;
	const CondLikeArrayHolder *GetHolder(int index);	
	CondLikeArrayHolder *GetMutableHolder(int index);
	
	int SetHolderDirty(int index);
	void IncrementHolder(int index);
	void DecrementHolder(int index);
	void SetHolderDependencies(int index, int depIndex1, int pDepIndex1, int depIndex2, int pDepIndex2);
	void SetCurReclaimableLevel(int lvl){curReclaimLevel = lvl;}
	void SetDepLevel(int index, int lev);
	void ClearFreeableQueue() {freeableHolderQueue.clear();}
	void AddToFreeableQueue(int index){
		if(index > -1){
			RemoveTempReservation(index);
			if(! IsHolderReserved(index))
				freeableHolderQueue.push_back(index);
			}
		}
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

	inline void ClaManager::_ReclaimH(int index){
		if( ! holders[index].IsDirty() )
			_ReclaimC(index);

		holders[index].Reset();
		holderStack.push_back(index);
		}
	
	inline void ClaManager::_ReclaimC(int index){
		claStack.push_back(holders[index].theArray);
		holders[index].theArray=NULL;
		}

	inline void ClaManager::SetDepLevel(int index, int lvl){
		holders[index].depLevel = lvl;
		if(debug_clas){
			char s[50];
			sprintf(s, "setdep\tlvl\t%d", lvl);
			ReportClaTotals(s, index);
			}
		}
/*
	inline int ClaManager::GetReclaimLevel(int index) const {
		assert(index > -1 && index < numHolders);
		if(holders[index].theArray == NULL) 
			return -1;
		return 
			holders[index].GetReclaimLevel();
		}
*/
	//This returns the cla number of the cla that this holder points to
	//it isn't changing anything, including indexes
	inline int ClaManager::GetClaIndexForBeagle(int index) const{
		assert(index > -1 && index < numHolders);
		assert(holders[index].theArray != NULL);
		return holders[index].theArray->Index();
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

	inline bool ClaManager::IsHolderDirty(int index) const {
		//dirtyness is synonymous with a null cla pointer in the holder
		assert(index > -1 && index < numHolders);
		assert(holders[index].numAssigned > 0);
		return (holders[index].theArray == NULL);	
		}

#endif