// GARLI version 0.96b8 source code
// Copyright 2005-2008 Derrick J. Zwickl
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

#include "defs.h"
#include "translatetable.h"
#include "datamatr.h"
#include "sequencedata.h"


void TranslateTable::SetTaxonName( int i, const char* s ){
	assert( Check(i-1) );
	assert( nTax );
	SetName(i-1, s);
	}

TranslateTable::TranslateTable( SequenceData* d ) : nTax(0){
	assert( d );
	nTax = d->NTax();
	Alloc();
	for( int i = 0; i < nTax; i++ )
		SetTaxonName( i+1, d->TaxonLabel(i) );
	}

void TranslateTable::Alloc(){
	assert( nTax );
	MEM_NEW_ARRAY(taxonName,char*,nTax);
	for( int i = 0; i < nTax; i++ )
	taxonName[i] = 0;
	}

void TranslateTable::Destroy(){
	for( int i = 0; i < nTax; i++ ) {
		int nmlen = (int)strlen( taxonName[i] );
		assert(nmlen > 0);
		MEM_DELETE_ARRAY(taxonName[i]); // taxonName[i] has length nmlen+1
		}
	MEM_DELETE_ARRAY(taxonName); // taxonName has length nTax
	taxonName = 0;
	nTax = 0;
	}

void TranslateTable::SetName( int i, const char* s ){
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

int TranslateTable::Find( const char* s ){
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

ostream& operator<<( ostream& out, TranslateTable& tt ){
	out << "translate" << endl;
	for( int i = 0; i < tt.nTax-1; i++ ) {
		out << "  " << (i+1) << ' ' << tt.taxonName[i] << ',' << endl;
		}
	out << "  " << tt.nTax << ' ' << tt.taxonName[tt.nTax-1] << endl;
	out << "  ;" << endl;
	return out;
	}
