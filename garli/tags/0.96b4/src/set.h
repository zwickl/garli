// set.h
// Copyright © 1998 by Paul O. Lewis
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//   Paul O. Lewis, Assistant Professor
//   167 Castetter Hall
//   Department of Biology
//   The University of New Mexico
//   Albuquerque, NM 87131-1091
//   Phone: (505) 277-6681
//   Fax: (505) 277-0304
//   email: lewisp@unm.edu
//   http://biology.unm.edu/~lewisp/pol.html
//
// Note: moving January 1, 1999, to the Department of Ecology and
// Evolutionary Biology, University of Connecticut
//
// Associated source code file: "set.cpp"
//

#ifndef __SET_H
#define __SET_H

#include <assert.h>

// Note: the data type unsigned char and the macro MISSING_DATA
// should be defined here the same way they are defined
// in the file "datamatr.h"

#define MISSING_DATA	(0xf)		

class Set
{
        int sz;
        int next;
	int* arr;

	void Realloc( int newsz );

	public:
		Set() : sz(0), arr(0), next(0) {}
		Set( int startsz );
		~Set();

		int Size() const { return next; }
		int Empty() const { return (sz==0); }

		Set& operator+=( const int i );
		Set& operator-=( const int i );
		int operator[]( const int i ) const;
};

class DNASet
{
	int set;

	public:
		enum {
			BASE_A = 0x01,
			BASE_C = 0x02,
			BASE_G = 0x04,
			BASE_T = 0x08,
			BASE_MISSING = 0x0f
		};

		DNASet() : set(0) {}
		DNASet( int s ) : set(s) {}

		int Empty() const { return (set==0); }
		void Flush() { set = 0; }

		DNASet& operator=( const DNASet& d )
			{ set = d.set; return *this; }
		DNASet& operator+=( const unsigned char i );
		DNASet& operator-=( const unsigned char i );

		// set intersection is mapped to the bitwise AND operator
		DNASet& operator|=( DNASet& b );
		friend DNASet operator|( DNASet& a, DNASet& b );

		// set union is mapped to the bitwise OR operator
		DNASet& operator&=( DNASet& b );
		friend DNASet operator&( DNASet& a, DNASet& b );
};

inline Set& Set::operator+=( const int i )
{
	if( next == sz ) Realloc( sz + 5 );
	arr[next++] = i;
	return *this;
}

inline int Set::operator[]( const int i ) const
{
	assert( i >= 0 );
	assert( i < next );
	return arr[i];
}

inline DNASet& DNASet::operator+=( const unsigned char i )
{
	// Note: this function is designed to take values of type unsigned char
	// see file datamatr.h before changing the relationship between
	// bases and unsigned char values, specifically the function DNAData::DatumToChar
	//
	// BUGBUG: unsigned char should be defined in its own header file along with
	// the necessary conversion functions such as DatumToChar and CharToDatum
	if( i == 0 )
		set |= BASE_A;
	else if( i == 1 )
		set |= BASE_C;
	else if( i == 2 )
		set |= BASE_G;
	else if( i == 3 )
		set |= BASE_T;
	else if( i == MISSING_DATA )
		set |= BASE_MISSING;

	return *this;
}

inline DNASet& DNASet::operator-=( const unsigned char i )
{
	// Note: this function is designed to take values of type unsigned char
	// see file datamatr.h before changing the relationship between
	// bases and unsigned char values, specifically the function DNAData::DatumToChar
	if( i == 0 )
		set &= ~BASE_A;
	else if( i == 1 )
		set &= ~BASE_C;
	else if( i == 2 )
		set &= ~BASE_G;
	else if( i == 3 )
		set &= ~BASE_T;
	else if( i == MISSING_DATA )
		set &= ~BASE_MISSING;	// this does nothing, but is aesthetically pleasing!

	return *this;
}

#endif
