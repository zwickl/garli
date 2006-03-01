// GARLI version 0.93 source code
// Copyright  2005 by Derrick J. Zwickl
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//  Derrick Zwickl
//	Integrative Biology, UT
//	1 University Station, C0930
//	Austin, TX  78712
//  email: zwickl@mail.utexas.edu
//
//	Note: In 2006  moving to NESCENT (The National
//	Evolutionary Synthesis Center) for a postdoc

//	NOTE: Portions of this source adapted from GAML source, written by Paul O. Lewis

#ifndef __MLHKY_H
#define __MLHKY_H

#include <vector>
using namespace std;

#include "datamatr.h"

class HKYData : public DNAData
{

	public:
		//DZ 9-2-04
		vector<char*> ambigStrings;
		
		HKYData() : DNAData() {}
		HKYData( int ntax, int nchar ) : DNAData( ntax, nchar ) {}
		~HKYData() {
			for(vector<char*>::iterator delit=ambigStrings.begin();delit!=ambigStrings.end();delit++)
				delete [](*delit);
			}
		void CalcEmpiricalFreqs( double* p);

		void MakeAmbigStrings();
		char *GetAmbigString(int i){
			return ambigStrings[i];
			}
		// overrides of base class's virtual fuctions
		double Freq( unsigned char d, int = 0);
};


#endif

