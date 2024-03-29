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

// --------
// Example showing how the array is laid out:
//	------------------------------------------
//	4 bases, 3 sites, 5 rate categories
//
//	          11111111112222222222333333333344444444445555555555
//	012345678901234567890123456789012345678901234567890123456789
//	------------------------------------------------------------
//	012301230123012301230123012301230123012301230123012301230123 <== base b
//	000011112222333344440000111122223333444400001111222233334444 <== rate r
//	000000000000000000001111111111111111111122222222222222222222 <== site s
//	indexing:             ^ this is s*4*5 + r*4 + b = 22
//
// If there is only one rate category...
//	4 bases, 3 sites
//
//	          11
//	012345678901
//	------------
//	012301230123 <== base b
//	000000000000 <== rate r
//	000011112222 <== site s
//	indexing: ^ this is s*4*1 + b = 10
//

#ifndef __CONDLIKE_H
#define __CONDLIKE_H

#include <vector>
#include <stddef.h>
using namespace std;

#include "defs.h"

//******************************************************************************
//  CondLikeArray
//
class CondLikeArray{
	//this is a CLA for a single model
	friend class CondLikeArrayIterator;

	unsigned nsites, nrates, nstates;
	int index;
	public:
		FLOAT_TYPE* arr;
		int* underflow_mult;
		unsigned rescaleRank;
		CondLikeArray(int nsit, int nsta, int nrat)
			: nsites(nsit), nrates(nrat), nstates(nsta), arr(NULL), underflow_mult(NULL), rescaleRank(1){}
		CondLikeArray()
			: nsites(0), nrates(0), nstates(0), arr(0), underflow_mult(0), rescaleRank(1){}
		~CondLikeArray();
		int NStates() const {
			return nstates;
			}
		int NChar() const {return nsites;}
		int NRateCats() const {return nrates;}
		int RequiredSize() const {return nsites * nstates * nrates;}
		void Assign(FLOAT_TYPE *alloc, int * under) {arr = alloc; underflow_mult = under;}
		int Index() { return index;}

		void Allocate( int nk, int ns, int nr = 1 );
		void Allocate(int ind, int nk, int ns, int nr = 1);
	};

class CondLikeArraySet{
	//this is a set of CLAs, one for each model
public:
		vector<CondLikeArray *> theSets;
		FLOAT_TYPE *rawAllocation;
		int *rawUnder;
		int index;

		CondLikeArraySet() : rawAllocation(NULL), rawUnder(NULL), index(-1){};
		CondLikeArraySet(int ind) : rawAllocation(NULL), rawUnder(NULL), index(ind) {};
		~CondLikeArraySet() {
			for(int i = 0;i < theSets.size();i++)
				delete theSets[i];
			theSets.clear();
			delete []rawAllocation;
			delete []rawUnder;
			}

		void Allocate();
		void AddCLA(CondLikeArray *cla ){
			theSets.push_back(cla);
			}
		CondLikeArray *GetCLA(int index){
			return theSets[index];
			}
		int Index() const{
			//assert(index >= 0);
			return index;
			}
	};

class CondLikeArrayHolder{
	friend class ClaManager;
	public:
	short numAssigned;
	short reclaimLevel;
	bool tempReserved;
	bool reserved;
	//CondLikeArray *theArray;
	CondLikeArraySet *theSet;
	int holderDep1;
	int holderDep2;
	int transMatDep1;
	int transMatDep2;
	int depLevel;
	
	//BEAGLEMERGE2018
	//These are mostly directly from the orig beagle branch, but replacing theArray with theSet//

	//CondLikeArrayHolder() : theSet(NULL), numAssigned(0), reclaimLevel(0), reserved(false) , tempReserved(false){}
	//~CondLikeArrayHolder() {theSet = NULL;}
	CondLikeArrayHolder() : depLevel(-1), theSet(NULL), numAssigned(0), reclaimLevel(0), reserved(false), tempReserved(false) {
		holderDep1 = holderDep2 = transMatDep1 = transMatDep2 = -1;
	}
	//the arrays are actually allocated and deleted by the claManager, not the holders
	~CondLikeArrayHolder() {
		theSet = NULL;
	};
	int HolderDep1() const { return holderDep1; }
	int HolderDep2() const { return holderDep2; }
	int TransMatDep1() const { return transMatDep1; }
	int TransMatDep2() const { return transMatDep2; }
	bool IsDirty() const { return theSet == NULL; }
	int DepLevel() const { return depLevel; }
	////

	int GetReclaimLevel() {return reclaimLevel;}
	void SetReclaimLevel(int lvl) {reclaimLevel = lvl;}
	void SetDepLevel(int lvl) { depLevel = lvl; }
	
	//BEAGLEMERGE2018
	//void Reset(){reclaimLevel=0;numAssigned=0,tempReserved=false;reserved=false;theSet=NULL;}
	void Reset() {
		depLevel = -1;
		theSet = NULL;
		numAssigned = 0;
		reclaimLevel = 0;
		tempReserved = false;
		reserved = false;
		holderDep1 = holderDep2 = transMatDep1 = transMatDep2 = -1;
	}
	////
	};
#endif

