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
using namespace std;

#include "defs.h"
class Model;

//******************************************************************************
//  CondLikeArray
//
class CondLikeArray
{
	friend class CondLikeArrayIterator;

	unsigned nsites, nrates, nstates;
	int index;
	public:
		FLOAT_TYPE* arr;
		int* underflow_mult;
		unsigned rescaleRank;
		CondLikeArray()
			: nsites(0), nrates(0), nstates(0), arr(0), underflow_mult(0), rescaleRank(1), index(-1){}
		~CondLikeArray();
		int NStates() {return nstates;}
		int NSites() {return nsites;}
		int Index() {return index;}

		void Allocate( int ind, int nk, int ns, int nr = 1 );
	};

class CondLikeArrayHolder{
	friend class ClaManager;

	short numAssigned;
	short reclaimLevel;
	bool tempReserved;
	bool reserved;
	CondLikeArray *theArray;
	int holderDep1;
	int holderDep2;
	int transMatDep1;	
	int transMatDep2;
	int depLevel;

public:

	CondLikeArrayHolder() : depLevel(-1), theArray(NULL), numAssigned(0), reclaimLevel(0), reserved(false) , tempReserved(false){
		holderDep1 = holderDep2 = transMatDep1 = transMatDep2 = -1;
		}
	//the arrays are actually allocated and deleted by the claManager, not the holders
	~CondLikeArrayHolder() {
		theArray = NULL;
		};
	int HolderDep1() const {return holderDep1;}
	int HolderDep2() const {return holderDep2;}
	int TransMatDep1() const {return transMatDep1;}
	int TransMatDep2() const {return transMatDep2;}
	bool IsDirty() const {return theArray == NULL;}
	int DepLevel() const {return depLevel;}

private:
	void SetReclaimLevel(int lvl) {reclaimLevel = lvl;}
	void SetDepLevel(int lvl) {depLevel = lvl;}

	void Reset(){
		depLevel = -1;
		theArray=NULL;
		numAssigned = 0;
		reclaimLevel = 0;
		tempReserved=false;
		reserved=false;
		holderDep1 = holderDep2 = transMatDep1 = transMatDep2 = -1;
		}
	};
#endif
