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

//******************************************************************************
//  CondLikeArray
//
class CondLikeArray
{
	friend class CondLikeArrayIterator;

	unsigned nsites, nrates, nstates;
	public:
		FLOAT_TYPE* arr;
		int* underflow_mult;
		unsigned rescaleRank;
		CondLikeArray()
			: nsites(0), nrates(0), nstates(0), arr(0), underflow_mult(0), rescaleRank(1){}
		~CondLikeArray();
		int NStates() {return nstates;}
		int NSites() {return nsites;}

		void Allocate( int nk, int ns, int nr = 1 );
	};

class CondLikeArrayHolder{
	public:
	short numAssigned;
	short reclaimLevel;
	bool tempReserved;
	bool reserved;
	CondLikeArray *theArray;
	CondLikeArrayHolder() : theArray(NULL), numAssigned(0), reclaimLevel(0), reserved(false) , tempReserved(false){}
	~CondLikeArrayHolder() {};
	int GetReclaimLevel() {return reclaimLevel;}
	void SetReclaimLevel(int lvl) {reclaimLevel = lvl;}
	void Reset(){reclaimLevel=0;numAssigned=0,tempReserved=false;reserved=false;theArray=NULL;}
	};
#endif

