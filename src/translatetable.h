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


#ifndef _TRANSTABLE
#define _TRANSTABLE

#include <cassert>
#include <string>
#include <iostream>

using namespace std;

class  HKYData;

class TranslateTable
{
	int nTax;
	char** taxonName;

	int Check( int i )
		{ return (i >= 0 && i<nTax ? 1 : 0 ); }
	void Alloc();
	void Destroy();
	void SetName( int i, char* s );

	public:
		TranslateTable( int n ) : nTax(n) { Alloc(); }
		TranslateTable( HKYData* d );
		~TranslateTable() { Destroy(); }

		void SetTaxonName( int i, char* s );
		char* GetTaxonName( int i )
			{ assert( Check(i-1) ); return taxonName[i-1]; }
		int GetNameLength( int i )
			{ assert( Check(i-1) ); return strlen( taxonName[i-1] ); }
		int Find( const char* s );

		friend ostream& operator<<( ostream& out, TranslateTable& tt );
};
#endif

