// GARLI version 0.94 source code
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

#include "translatetable.h"
#include "datamatr.h"
#include "mlhky.h"
#include "defs.h"

void TranslateTable::SetTaxonName( int i, char* s )
{
        assert( Check(i-1) );
        assert( nTax );
        SetName(i-1, s);
}

TranslateTable::TranslateTable( HKYData* d ) : nTax(0)
{
        assert( d );
        nTax = d->NTax();
	Alloc();
        for( int i = 0; i < nTax; i++ )
                SetTaxonName( i+1, d->TaxonLabel(i) );
}

void TranslateTable::Alloc()
{
        assert( nTax );
        MEM_NEW_ARRAY(taxonName,char*,nTax);
        for( int i = 0; i < nTax; i++ )
                taxonName[i] = 0;
}

void TranslateTable::Destroy()
{
        for( int i = 0; i < nTax; i++ ) {
		    int nmlen = (int)strlen( taxonName[i] );
            assert(nmlen > 0);
			MEM_DELETE_ARRAY(taxonName[i]); // taxonName[i] has length nmlen+1
			}
        MEM_DELETE_ARRAY(taxonName); // taxonName has length nTax
        taxonName = 0;
        nTax = 0;
}

void TranslateTable::SetName( int i, char* s )
{
        assert(s);
        int nmlen;
        if( taxonName[i] ) {
                nmlen = (int)strlen( taxonName[i] );
                MEM_DELETE_ARRAY(taxonName[i]); // taxonName[i] has length nmlen+1
        }
        nmlen = (int)strlen(s);
	MEM_NEW_ARRAY(taxonName[i],char,nmlen+1);
        assert( taxonName[i] );
        strcpy( taxonName[i], s );
}

int TranslateTable::Find( const char* s )
{
        assert(s);
        int taxonNumber = 0;
        for( int i = 0; i < nTax; i++ ) {
                if( strcmp( taxonName[i], s ) == 0 ) {
                        taxonNumber = i+1;
                        break;
                }
	}

	return taxonNumber;
}

ostream& operator<<( ostream& out, TranslateTable& tt )
{
	out << "translate" << endl;
	for( int i = 0; i < tt.nTax-1; i++ ) {
		out << "  " << (i+1) << ' ' << tt.taxonName[i] << ',' << endl;
	}
	out << "  " << tt.nTax << ' ' << tt.taxonName[tt.nTax-1] << endl;
	out << "  ;" << endl;
	return out;
}
