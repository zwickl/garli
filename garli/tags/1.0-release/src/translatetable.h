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

#ifndef _TRANSTABLE
#define _TRANSTABLE

#include <cassert>
#include <string>
#include <cstring>
#include <iostream>

using namespace std;

class  SequenceData;

class TranslateTable
{
	int nTax;
	char** taxonName;

	int Check( int i )
		{ return (i >= 0 && i<nTax ? 1 : 0 ); }
	void Alloc();
	void Destroy();
	void SetName( int i, const char* s );

	public:
		TranslateTable( int n ) : nTax(n) { Alloc(); }
		TranslateTable( SequenceData* d );
		~TranslateTable() { Destroy(); }

		void SetTaxonName( int i, const char* s );
		char* GetTaxonName( int i )
			{ assert( Check(i-1) ); return taxonName[i-1]; }
		int GetNameLength( int i )
			{ assert( Check(i-1) ); return (int)strlen( taxonName[i-1] ); }
		int Find( const char* s );

		friend ostream& operator<<( ostream& out, TranslateTable& tt );
};
#endif

