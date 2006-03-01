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

#include <cstring>

using namespace std;

#include "set.h"

DNASet& DNASet::operator|=( DNASet& b )
{
	set |= b.set;
	return *this;
}

DNASet operator|( DNASet& a, DNASet& b )
{
	return DNASet( a.set | b.set );
}

DNASet& DNASet::operator&=( DNASet& b )
{
	set &= b.set;
	return *this;
}

DNASet operator&( DNASet& a, DNASet& b )
{
	return DNASet( a.set & b.set );
}

Set::Set( int startsz ) : sz(startsz), next(0)
{
	arr=new int[startsz];
}

Set::~Set()
{
	delete []arr;
}

Set& Set::operator-=( const int i )
{
        for( int j = 0; j < next; j++ ) {
                if( arr[j] != i ) continue;

                // arr[j] equals i
                arr[j] = arr[--next];
                arr[next] = 0;
                break;
        }
        return *this;
}

void Set::Realloc( int newsz )
{
        if( newsz <= sz ) return;
        int* newarr;
		newarr=new int[newsz];
        memset( newarr, 0x00, newsz*sizeof(int) );
        for( int i = 0; i < sz; i++ )
                newarr[i] = arr[i];
        delete []arr;
		sz = newsz;
        arr = newarr;
}


