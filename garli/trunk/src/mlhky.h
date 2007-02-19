// GARLI version 0.951 source code
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

#ifndef __MLHKY_H
#define __MLHKY_H

#include <vector>
using namespace std;

#include "memchk.h"
#include "datamatr.h"
#include "defs.h"

class HKYData : public DNAData
{

	public:
		//DZ 9-2-04
		vector<char*> ambigStrings;
#ifdef OPEN_MP
		vector<unsigned*> ambigToCharMap;
#endif
		FLOAT_TYPE *empStateFreqs;
		
		HKYData() : DNAData() {}
		HKYData( int ntax, int nchar ) : DNAData( ntax, nchar ) {empStateFreqs=NULL;}
		~HKYData() {
			for(vector<char*>::iterator delit=ambigStrings.begin();delit!=ambigStrings.end();delit++)
				delete [](*delit);
			}
		void CalcEmpiricalFreqs();
		void GetEmpiricalFreqs(FLOAT_TYPE *f) const{
			for(int i=0;i<4;i++) f[i]=empStateFreqs[i];
			}
	
		void MakeAmbigStrings();
		char *GetAmbigString(int i) const{
			return ambigStrings[i];
			}
#ifdef OPEN_MP
		unsigned *GetAmbigToCharMap(int i) const{
			return ambigToCharMap[i];
			}
#endif
		// overrides of base class's virtual fuctions
		FLOAT_TYPE Freq( unsigned char d, int = 0);
};


#endif

