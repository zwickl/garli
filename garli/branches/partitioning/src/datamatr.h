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
//
//	NOTE: Portions of this source adapted from GAML source, written by Paul O. Lewis

#ifndef __DATAMATR_H
#define __DATAMATR_H

#include <iostream>
#include <cassert>
#include <stdio.h>
#include <math.h>
using namespace std;

#include "ncl.h"
#include "errorexception.h"

class GarliReader;

typedef FLOAT_TYPE** DblPtrPtr;
#define MAX_STATES (8*sizeof(unsigned char))
#define FIRST_STATE	(0x01)
#define LAST_STATE	(0x80)
#define MISSING_DATA	(0xf)		// all bits set to 1

#if defined( CPLUSPLUS_EXCEPTIONS )
#  define THROW_BADSTATE(a) throw XBadState(a)
#else
#  define THROW_BADSTATE(a) BadState(a)
#endif

// Note: the class below has pure virtual member functions
class DataMatrix{
protected:
	int		nTax;
	int		nChar;     //after compression
	int		totalNChar; //columns in matrix, with all missing columns removed
	int 	gapsIncludedNChar; //the actual number of columns in the data matrix read in
								//only used when outputting something relative to input alignment
	int		dense;		//whether the data has been sorted and identical patterns combined
	
	unsigned char**         matrix;
	int*		count;
	int*		origCounts;
	int*		number; //maping of chars to columns
						//Indeces are original char numbers,
						//contents are the packed column representing that char
						//both start at 0, so offset upon output
						//This used to represent something else (I think)
	char**          taxonLabel;
	char**          taxonColor;
	int		nMissing;
	int		nConstant;
	int		nInformative;
	int		nAutapomorphic;
	int 	lastConstant;//DJZ
	int 	*constStates;//the state (or states) that a constant site contains
	FLOAT_TYPE*	stateDistr;
	int		stateDistrComputed;
	int		currentBootstrapSeed;
	
	protected:
		int*	numStates;
		int		dmFlags;
		int     maxNumStates;

	protected:
		char	info[80];
		virtual void	SwapCharacters( int i, int j );
		virtual int	ComparePatterns( const int i, const int j ) const;
		void	BSort( int byCounts = 0 );
		void	DebugSaveQSortState( int top, int bottom, int ii, int jj, int xx, const char* title );
		void	QSort( int top, int bottom );
		void	ReplaceTaxonLabel( int i, const char* s );
		void	ReplaceTaxonColor( int i, const char* s );

	public:
		enum {STANDARD, DNA, RNA, PROTEIN };
		enum {	// for dmFlags variable
			allvariable = 0x0001,	// all characters in data matrix variable
			vnstates    = 0x0002,	// characters vary in maximum number of states
			ambigstates = 0x0004	// ambiguous states found in data matrix
		};
		enum {
			PT_MISSING		= 0x0000,
			PT_CONSTANT		= 0x0001,
			PT_INFORMATIVE		= 0x0002,
			PT_AUTAPOMORPHIC	= 0x0004,
			PT_VARIABLE		= 0x0008
		};

	public:
		DataMatrix() : dmFlags(0), dense(0), nTax(0), nChar(0), matrix(0), count(0)
			, number(0), taxonLabel(0), numStates(0), stateDistr(0)
			, nMissing(0), nConstant(0), nInformative(0), nAutapomorphic(0), stateDistrComputed(0),
			lastConstant(-1), constStates(0), origCounts(0), currentBootstrapSeed(0)
			{ memset( info, 0x00, 80 ); }
		DataMatrix( int ntax, int nchar )
			: nTax(ntax), nChar(nchar), dmFlags(0), dense(0), matrix(0), count(0)
			, number(0), taxonLabel(0), numStates(0), stateDistr(0)
			, nMissing(0), nConstant(0), nInformative(0), nAutapomorphic(0), stateDistrComputed(0),
			lastConstant(-1), constStates(0), origCounts(0), currentBootstrapSeed(0)
			{ memset( info, 0x00, 80 ); NewMatrix(ntax, nchar); }
		virtual ~DataMatrix();

		// pure virtual functions - must override in derived class
		virtual unsigned char	CharToDatum( char ch )					= 0;
		virtual unsigned char CharToBitwiseRepresentation( char ch ) 	= 0;
		virtual char	DatumToChar( unsigned char d )	                = 0;
		virtual unsigned char	FirstState()		                    = 0;
		virtual unsigned char	LastState()		                        = 0;
//		virtual FLOAT_TYPE	Freq( unsigned char, int = 0)	                        = 0;

		// virtual functions - can override in derived class
		virtual FLOAT_TYPE	TransitionProb( int /*i*/, int /*j*/
			, int /*site*/, FLOAT_TYPE /*brlen*/) { return 0.0; }
		virtual int NumStates(int j) const
			{ return ( numStates && (j < nChar) ? numStates[j] : 0 ); }

		FLOAT_TYPE prNumStates( int n ) const;

		// functions for quizzing dmFlags
		int InvarCharsExpected() { return !(dmFlags & allvariable); }
		int VariableNumStates() { return (dmFlags & vnstates); }
		int AmbiguousStates() { return (dmFlags & ambigstates); }

		// functions for getting the data in and out
		int GetToken( istream& in, char* tokenbuf, int maxlen, bool acceptComments=true );
		int GetToken( FILE *in, char* tokenbuf, int maxlen);
		int ReadPhylip( const char* filename);
		int ReadFasta( const char* filename);
		int Save( const char* filename, char* newfname = 0, char* nxsfname = 0 );

		char*	DataType() { return info; }
		int     unsigned charToInt( unsigned char d ) { return (int)d; }

		int NTax() const { return nTax; }
		void SetNTax(int ntax) { nTax = ntax; }

		virtual int NChar() const { return nChar; }
		int TotalNChar() const { return totalNChar; }
		int GapsIncludedNChar() const { return gapsIncludedNChar; }
		void SetNChar(int nchar) { nChar = nchar; }

		void Flush() { NewMatrix( 0, 0 ); }
		int Dense() const { return dense; }
		
		int Number(int j) const
			//this appears to have been a bug.  Should be gapsIncludedNChar
			//{ return ( number && (j < totalNChar) ? number[j] : 0 ); }
			{ return ( number && (j < gapsIncludedNChar) ? number[j] : 0 ); }

		virtual int Count(int j) const
			{ return ( count && (j < nChar) ? count[j] : 0 ); }
		virtual const int *GetCounts() const {return count;}
		const int *GetConstStates() const {return constStates;}
		void SetCount(int j, int c)
			{ if( count && (j < nChar) ) count[j] = c; }

		void SetNumStates(int j, int c)
			{ if( numStates && (j < nChar) ) numStates[j] = c; }

		const char* TaxonLabel(int i) const{
			return ( taxonLabel && (i < nTax) ? taxonLabel[i] : 0 );
			}
		void SetTaxonLabel(int i, const char* s);
		
		int TaxonNameToNumber(const NxsString &name) const{
			for(int i=0;i<nTax;i++){
				if(strcmp(name.c_str(), TaxonLabel(i)) == 0) return i+1;//indeces run 0->ntax-1, taxon numbers 1->ntax
				}
			return -1;
			}
		void CopyNamesFromOtherMatrix(const DataMatrix *dat){
			assert(taxonLabel);
			for(int t=0;t<nTax;t++) SetTaxonLabel(t, dat->TaxonLabel(t));
			}

		void BeginNexusTreesBlock(ofstream &treeout) const;
		virtual void CreateMatrixFromNCL(NxsCharactersBlock *) = 0;
	
		virtual unsigned char Matrix( int i, int j ) const {
			assert( matrix );
			assert( i >= 0 );
			assert( i < nTax );
			assert( j >= 0 );
			assert( j < nChar );
			return (unsigned char)matrix[i][j];
			}

		unsigned char *GetRow( int i) const {
			assert( matrix );
			assert( i >= 0 );
			assert( i < nTax );
			return matrix[i];
		}
		virtual void SetMatrix( int i, int j, unsigned char c )
			{ if( matrix && (i < nTax) && (j < nChar) ) matrix[i][j] = c; }

		int MatrixExists() const { return ( matrix && nTax>0 && nChar>0 ? 1 : 0 ); }
		int NMissing() const { return nMissing; }
		int NConstant() const { return nConstant; }
		int LastConstant() const {return lastConstant;}
		int NInformative() const { return nInformative; }
		int NAutapomorphic() const { return nAutapomorphic; }

		DataMatrix& operator =(const DataMatrix&);

		void Sort( int byCounts = 0 ){
			byCounts;
			QSort( 0, NChar()-1 );
			}
		virtual int PatternType( int , int*, unsigned char *) const;	// returns PT_XXXX constant indicating type of pattern
		void Summarize();       // fills in nConstant, nInformative, and nAutapomorphic data members
		virtual void Collapse();
		virtual void Pack();
		void NewMatrix(int nt, int nc);	// flushes old matrix, creates new one
		int PositionOf( char* s ) const; // returns pos (0..nTax-1) of taxon named s
		void AllocPr( DblPtrPtr& pr );
		void DeletePr( DblPtrPtr& pr );
		void DumpCounts( const char* s );
		void WriteCollapsedData();  //DZ
		void SaveNexus(const char* filename, int iosFlags /* = 0 */); //DZ
		void DetermineConstantSites();
		int Serialize(char**, int*);  // cjb
		int Deserialize(const char*, const int);  // cjb
		bool operator==(const DataMatrix& rhs) const; // cjb - to test serialization
		void ExplicitDestructor();  // cjb - totally clear the DataMatrix and revert it to its original state as if it was just constructed
		void CheckForIdenticalTaxonNames();

	public:	// exception classes
#if defined( CPLUSPLUS_EXCEPTIONS )
		class XBadState {
			public:
				XBadState( char c ) : ch(c) {}
				char ch;
		};
#else
      void Abort( char* msg )
         { cerr << endl << "Error:" << msg << endl << "Program terminated." << endl; 
         #ifdef POWERMAC_VERSION
         	assert(0);
         	cout<<"quit"<<endl;
         	char c;
         	cin>>c;
         	throw 1;
         #else
         	exit(1); 
         #endif
         }
      void BadState( char ch ) { char s[80]; sprintf(s, "Bad character state (%c)", ch); Abort(s); }
      
      void ReserveOriginalCounts(){
      		if(origCounts==NULL) origCounts=new int[nChar];
      		for(int i=0;i<nChar;i++){
      			origCounts[i]=count[i];
      			}
      		}
      void RestoreOriginalCounts(){
			if(origCounts==NULL) return;
      		for(int i=0;i<nChar;i++){
      			count[i]=origCounts[i];
      			}
      		}
      void Reweight(FLOAT_TYPE prob);
      long BootstrapReweight(int seed, FLOAT_TYPE resampleProportion);
      
#endif
};

//
// BinaryData implements a simple two-state poisson process model
// described by Felsenstein, J. 1981. A likelihood approach to character
// weighting and what it tells us about parsimony and compatibility.
// Biol. J. Linnean Soc. 16: 183-196.
//
class BinaryData : public DataMatrix
{
	public:
		BinaryData() : DataMatrix()
                        { maxNumStates=2; dmFlags = allvariable; strcpy( info, "binary" ); }
		BinaryData( int ntax, int nchar ) : DataMatrix( ntax, nchar )
                        { maxNumStates=2; strcpy( info, "binary" ); }
		~BinaryData() {}

		// overrides of base class's virtual fuctions
		char		DatumToChar( unsigned char d );
		virtual unsigned char 	CharToDatum( char ch );
		virtual unsigned char	FirstState() { return FIRST_STATE; }
		virtual unsigned char	LastState() { return (FIRST_STATE<<1); }
//		virtual FLOAT_TYPE	Freq( unsigned char, int = 0) { return 0.5; }
#if defined( RATIO_FOR_EACH_BRANCH )
		virtual void	CalcPrMatrix( FLOAT_TYPE** pr, FLOAT_TYPE brlen, FLOAT_TYPE kappa );
#else
		virtual void	CalcPrMatrix( FLOAT_TYPE** pr, FLOAT_TYPE brlen );
#endif
		virtual int	NumStates(int) const { return 2; }
};

inline unsigned char BinaryData::CharToDatum( char ch )
{
	unsigned char datum;

	if( ch == '0' || ch == '-' )
		datum = 0x01;
	else if( ch == '1' || ch == '+' )
		datum = 0x02;
	else if( ch == '?' )
		datum = MISSING_DATA;
	else
      THROW_BADSTATE(ch);

	return datum;
}

inline char BinaryData::DatumToChar( unsigned char d )
{
	char ch = 'X';	// ambiguous

	if( d == MISSING_DATA )
		ch = '?';
	else if( d == 0x01 )
		ch = '0';
	else if( d == 0x02 )
		ch = '1';

	return ch;
}

#if defined( RATIO_FOR_EACH_BRANCH )
inline void BinaryData::CalcPrMatrix( FLOAT_TYPE** pr, FLOAT_TYPE brlen, FLOAT_TYPE /*kappa*/ )
#else
inline void BinaryData::CalcPrMatrix( FLOAT_TYPE** pr, FLOAT_TYPE brlen )
#endif
{
	FLOAT_TYPE prNoChg = (FLOAT_TYPE) (0.5 + 0.5 * exp((FLOAT_TYPE) -2.0 * brlen ));
	pr[0][0] = pr[1][1] = prNoChg;
	FLOAT_TYPE prChg = (FLOAT_TYPE) (0.5 - 0.5 * exp((FLOAT_TYPE) -2.0 * brlen ));
	pr[0][1] = pr[1][0] = prChg;
}

//
// MorphData is identical to BinaryData except that it allows each
// character to have a different maximum number of states.  As a
// result of this generality, it is much slower than BinaryData.
//
class MorphData : public DataMatrix
{
	public:
		MorphData() : DataMatrix()
                        { maxNumStates=2; dmFlags = ( vnstates | allvariable ); strcpy( info, "morph" ); }
		MorphData( int ntax, int nchar )
			: DataMatrix( ntax, nchar )
                        { maxNumStates=2; strcpy( info, "morph" ); }
		~MorphData() {}

		// overrides of base class's virtual fuctions
		char		DatumToChar( unsigned char d );
		virtual unsigned char 	CharToDatum( char ch );
		virtual unsigned char	FirstState() { return 0x01; }
		virtual unsigned char	LastState() { return 0x08; }
//		virtual FLOAT_TYPE	Freq( unsigned char, int i = 0) { return (1.0 / (FLOAT_TYPE)numStates[i]); }
#if defined( RATIO_FOR_EACH_BRANCH )
		virtual void	CalcPrMatrix( FLOAT_TYPE** /*pr*/, FLOAT_TYPE /*brlen*/, FLOAT_TYPE /*kappa*/ ) {}
#else
		virtual void	CalcPrMatrix( FLOAT_TYPE** /*pr*/, FLOAT_TYPE /*brlen*/ ) {}
#endif

		// transition probability from state i to state j for character k
		// given branch length brlen
		// i and j range from 0, 1, 2, 3 (i.e.,they are not bits)
		// this function called by a node only if its pr matrix does not exist
		virtual FLOAT_TYPE	TransitionProb( int i, int j, int k, FLOAT_TYPE brlen );
};

// prob. i --> j when there are k total states possible
inline FLOAT_TYPE MorphData::TransitionProb( int i, int j, int k, FLOAT_TYPE brlen )
{
	FLOAT_TYPE prob = 0.0;
	assert( k > 1 );
   FLOAT_TYPE betat = (FLOAT_TYPE) (k * brlen / ( k - 1.0 ));
   prob = (FLOAT_TYPE)(( 1.0 - exp( -betat ) ) / k);
   if( i == j )
      prob = (FLOAT_TYPE)(1.0 - ( (FLOAT_TYPE)k - 1.0) * prob);
	return prob;
}

inline unsigned char MorphData::CharToDatum( char ch )
{
	unsigned char datum;

	if( ch == '0' )
		datum = 0x01;
	else if( ch == '1' )
		datum = 0x02;
	else if( ch == '2' )
		datum = 0x04;
	else if( ch == '3' )
		datum = 0x08;
	else if( ch == '4' )
		datum = 0x10;
	else if( ch == '5' )
		datum = 0x20;
	else if( ch == '6' )
		datum = 0x40;
	else if( ch == '7' )
		datum = 0x80;
	else if( ch == '?' )
		datum= MISSING_DATA;
	else
      THROW_BADSTATE(ch);

	return datum;
}

inline char MorphData::DatumToChar( unsigned char d )
{
	char ch = 'X';	// ambiguous

	if( d == MISSING_DATA )
		ch = '?';
	else if( d == 0x01 )
		ch = '0';
	else if( d == 0x02 )
		ch = '1';
	else if( d == 0x04 )
		ch = '2';
	else if( d == 0x08 )
		ch = '3';
	else if( d == 0x10 )
		ch = '4';
	else if( d == 0x20 )
		ch = '5';
	else if( d == 0x40 )
		ch = '6';
	else if( d == 0x80 )
		ch = '7';

	return ch;
}

#endif

