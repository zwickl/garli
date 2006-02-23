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

#ifndef __DATAMATR_H
#define __DATAMATR_H

#include <iostream>
#include <cassert>
#include <stdio.h>
#include <math.h>
using namespace std;

#include "stricl.h"

typedef double** DblPtrPtr;
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
class DataMatrix
{
	int		nTax;
	int		nChar;     //after compression
	int		totalNChar; //columns in matrix
	int		dense;
	
	unsigned char**         matrix;
	int*		count;
	int*		origCounts;
	int*		number;
	char**          taxonLabel;
	char**          taxonColor;	//POL-2/19/98
	int		nConstant;
	int 	lastConstant;//DJZ
	int 	*constBases;//DJZ
	int		nInformative;
	int		nAutapomorphic;
	double*         stateDistr;
	int		stateDistrComputed;
	
	protected:
		int*	numStates;
		int	dmFlags;
		int     maxNumStates;

	protected:
		char	info[80];
		void	SwapCharacters( int i, int j );
		int	ComparePatterns( const int i, const int j ) const;
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
			PT_CONSTANT		= 0x0001,
			PT_INFORMATIVE		= 0x0002,
			PT_AUTAPOMORPHIC	= 0x0004,
			PT_VARIABLE		= 0x0008
		};

	public:
		DataMatrix() : dmFlags(0), dense(0), nTax(0), nChar(0), matrix(0), count(0)
			, number(0), taxonLabel(0), numStates(0), stateDistr(0)
			, nConstant(0), nInformative(0), nAutapomorphic(0), stateDistrComputed(0),
			constBases(0), origCounts(0)
			{ memset( info, 0x00, 80 ); }
		DataMatrix( int ntax, int nchar )
			: nTax(ntax), nChar(nchar), dmFlags(0), dense(0), matrix(0), count(0)
			, number(0), taxonLabel(0), numStates(0), stateDistr(0)
			, nConstant(0), nInformative(0), nAutapomorphic(0), stateDistrComputed(0),
			constBases(0), origCounts(0)
			{ memset( info, 0x00, 80 ); NewMatrix(ntax, nchar); }
		~DataMatrix();

		// pure virtual functions - must override in derived class
		virtual unsigned char	CharToDatum( char ch )	                        = 0;
		virtual unsigned char CharToBitwiseRepresentation( char ch ) 			=0;
		virtual char	DatumToChar( unsigned char d )	                        = 0;
		virtual unsigned char	FirstState()		                        = 0;
		virtual unsigned char	LastState()		                        = 0;
//		virtual double	Freq( unsigned char, int = 0)	                        = 0;
/*
#if defined( RATIO_FOR_EACH_BRANCH )
		virtual void	CalcPrMatrix( double** pr, double brlen, double kappa )       = 0;
#else
		virtual void	CalcPrMatrix( double** pr, double brlen )       = 0;
#endif
*/
		// virtual functions - can override in derived class
		virtual double	TransitionProb( int /*i*/, int /*j*/
			, int /*site*/, double /*brlen*/) { return 0.0; }
		virtual int NumStates(int j) const
			{ return ( numStates && (j < nChar) ? numStates[j] : 0 ); }

		double prNumStates( int n ) const;

		// functions for quizzing dmFlags
		int InvarCharsExpected() { return !(dmFlags & allvariable); }
		int VariableNumStates() { return (dmFlags & vnstates); }
		int AmbiguousStates() { return (dmFlags & ambigstates); }

		// functions for getting the data in and out
		int GetToken( istream& in, char* tokenbuf, int maxlen );
		int Read( const char* filename, char* left_margin = 0 );
		int Save( const char* filename, char* newfname = 0, char* nxsfname = 0 );

		char*	DataType() { return info; }
		int     unsigned charToInt( unsigned char d ) { return (int)d; }

		int NTax() const { return nTax; }
		void SetNTax(int ntax) { nTax = ntax; }

		int NChar() const { return nChar; }
		void SetNChar(int nchar) { nChar = nchar; }

		void Flush() { NewMatrix( 0, 0 ); }
		int Dense() const { return dense; }
		
		int Number(int j) const
			{ return ( number && (j < nChar) ? number[j] : 0 ); }

		int Count(int j) const
			{ return ( count && (j < nChar) ? count[j] : 0 ); }
		const int *GetCounts() const {return count;}
		const int *GetConstBases() const {return constBases;}
		void SetCount(int j, int c)
			{ if( count && (j < nChar) ) count[j] = c; }

		void SetNumStates(int j, int c)
			{ if( numStates && (j < nChar) ) numStates[j] = c; }

		char* TaxonLabel(int i) const{
			return ( taxonLabel && (i < nTax) ? taxonLabel[i] : 0 );
			}
		void SetTaxonLabel(int i, const char* s);
		
		int TaxonNameToNumber(const NxsString &name) const{
			for(int i=0;i<nTax;i++){
				if(Strcmp(name, TaxonLabel(i)) == 0) return i+1;//indeces run 0->ntax-1, taxon numbers 1->ntax
				}
			return -1;
			}
		
		void BeginNexusTreesBlock(ofstream &treeout);

		// added 3/20/1998 in order to implement the Muse and Gaut codon model
//		virtual int NCodons() { return 0; }
//		virtual unsigned char CodonState( int i, int j ) { return (unsigned char)0; }
		
		unsigned char Matrix( int i, int j ) const {
//			assert( matrix );
//			assert( i >= 0 );
//			assert( i < nTax );
//			assert( j >= 0 );
//			assert( j < nChar );
			return (unsigned char)matrix[i][j];
		}
				//MTH
		unsigned char *GetRow( int i) const {
			assert( matrix );
			assert( i >= 0 );
			assert( i < nTax );
			return matrix[i];
		}
		void SetMatrix( int i, int j, unsigned char c )
			{ if( matrix && (i < nTax) && (j < nChar) ) matrix[i][j] = c; }

		int MatrixExists() const { return ( matrix && nTax>0 && nChar>0 ? 1 : 0 ); }
		int NConstant() { return nConstant; }
		int LastConstant() const {return lastConstant;}
		int NInformative() { return nInformative; }
		int NAutapomorphic() { return nAutapomorphic; }

		DataMatrix& operator =(const DataMatrix&);

		void Sort( int byCounts = 0 ){
			byCounts;
			QSort( 0, nChar-1 );
			}
		int PatternType( int , int*, unsigned char *) const;	// returns PT_XXXX constant indicating type of pattern
		void Summarize();       // fills in nConstant, nInformative, and nAutapomorphic data members
		void Collapse();
		void Pack();
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
      		for(int i=0;i<nChar;i++){
      			count[i]=origCounts[i];
      			}
      		}
      void Reweight(double prob);
      void BootstrapReweight();
      
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
//		virtual double	Freq( unsigned char, int = 0) { return 0.5; }
#if defined( RATIO_FOR_EACH_BRANCH )
		virtual void	CalcPrMatrix( double** pr, double brlen, double kappa );
#else
		virtual void	CalcPrMatrix( double** pr, double brlen );
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
inline void BinaryData::CalcPrMatrix( double** pr, double brlen, double /*kappa*/ )
#else
inline void BinaryData::CalcPrMatrix( double** pr, double brlen )
#endif
{
	double prNoChg = 0.5 + 0.5 * exp( -2.0 * brlen );
	pr[0][0] = pr[1][1] = prNoChg;
	double prChg = 0.5 - 0.5 * exp( -2.0 * brlen );
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
//		virtual double	Freq( unsigned char, int i = 0) { return (1.0 / (double)numStates[i]); }
#if defined( RATIO_FOR_EACH_BRANCH )
		virtual void	CalcPrMatrix( double** /*pr*/, double /*brlen*/, double /*kappa*/ ) {}
#else
		virtual void	CalcPrMatrix( double** /*pr*/, double /*brlen*/ ) {}
#endif

		// transition probability from state i to state j for character k
		// given branch length brlen
		// i and j range from 0, 1, 2, 3 (i.e.,they are not bits)
		// this function called by a node only if its pr matrix does not exist
		virtual double	TransitionProb( int i, int j, int k, double brlen );
};

// prob. i --> j when there are k total states possible
inline double MorphData::TransitionProb( int i, int j, int k, double brlen )
{
	double prob = 0.0;
	assert( k > 1 );
   double betat = (double)k * brlen / ( (double)k - 1.0 );
   prob = ( 1.0 - exp( -betat ) ) / (double)k;
   if( i == j )
      prob = 1.0 - ( (double)k - 1.0) * prob;
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

// ********************* new material above *****************************

class DNAData : public DataMatrix
{
	public:
		DNAData() : DataMatrix()
			{ maxNumStates=4; strcpy( info, "DNA" ); }
		DNAData( int ntax, int nchar ) : DataMatrix( ntax, nchar )
			{ maxNumStates=4; strcpy( info, "DNA" ); }
		~DNAData() {}

		// overrides of base class's virtual fuctions
		virtual unsigned char 	CharToDatum( char ch );
		virtual unsigned char CharToBitwiseRepresentation( char ch );
		virtual char	DatumToChar( unsigned char d );
		virtual unsigned char	FirstState() { return 0; }
		virtual unsigned char	LastState() { return 3; }
		virtual int	NumStates(int) const { return 4; }

		
};

inline unsigned char DNAData::CharToDatum( char ch )
{
	unsigned char datum;

	if( ch == 'A' || ch == 'a' )
		datum = 0;
	else if( ch == 'C' || ch == 'c' )
		datum = 1;
	else if( ch == 'G' || ch == 'g' )
		datum  = 2;
	else if( ch == 'T' || ch == 't' )
		datum = 3;
	else if( ch == '?' || ch == '-' )
		datum = MISSING_DATA;
	else if( strchr( "rRyYmMkKsSwWhHbBvVdDnN", ch ) ) {
		datum = MISSING_DATA;
		dmFlags |= ambigstates;
	}
	else
		THROW_BADSTATE(ch);

	return datum;
}

inline unsigned char DNAData::CharToBitwiseRepresentation( char ch )
{
	unsigned char datum=0;
	switch(ch){
		case 'A' : datum=1; break;
    	case 'C' : datum=2; break;
       	case 'G' : datum=4; break;
       	case 'T' : datum=8; break;
       	case 'U' : datum=8; break;
       	case 'M' : datum=3; break;
       	case 'R' : datum=5; break;
       	case 'S' : datum=6; break;
       	case 'V' : datum=7; break;
       	case 'W' : datum=9; break;
       	case 'Y' : datum=10; break;
       	case 'H' : datum=11; break;
      	case 'K' : datum=12; break;
       	case 'D' : datum=13; break;
       	case 'B' : datum=14; break;
       	case 'N' : datum=15; break;

		case 'a' : datum=1; break;
    	case 'c' : datum=2; break;
       	case 'g' : datum=4; break;
       	case 't' : datum=8; break;
	case 'u' : datum=8; break;
       	case 'm' : datum=3; break;
       	case 'r' : datum=5; break;
       	case 's' : datum=6; break;
       	case 'v' : datum=7; break;
       	case 'w' : datum=9; break;
       	case 'y' : datum=10; break;
       	case 'h' : datum=11; break;
      	case 'k' : datum=12; break;
       	case 'd' : datum=13; break;
       	case 'b' : datum=14; break;
       	case 'n' : datum=15; break;
       	case '-' : datum=15; break;
  		case '?' : datum=15; break;
       	default  : assert(0);
		}
	return datum;
}


inline char DNAData::DatumToChar( unsigned char d )
{
	char ch;
	switch(d){
		case 1 : ch='A'; break;
    	case 2 : ch='C'; break;
       	case 4 : ch='G'; break;
       	case 8 : ch='T'; break;
 //      	case 8 : ch='U'; break;
       	case 3 : ch='M'; break;
       	case 5 : ch='R'; break;
       	case 6 : ch='S'; break;
       	case 7 : ch='V'; break;
       	case 9 : ch='W'; break;
       	case 10 : ch='Y'; break;
       	case 11 : ch='H'; break;
      	case 12 : ch='K'; break;
       	case 13 : ch='D'; break;
       	case 14 : ch='B'; break;
       	case 15 : ch='?'; break;
		default  : assert(0);
		}
/*		
	char ch = 'X';

	if( d == MISSING_DATA )
		ch = '?';
	else if( d == 0 )
		ch = 'A';
	else if( d == 1 )
		ch = 'C';
	else if( d == 2 )
		ch = 'G';
	else if( d == 3 )
		ch = 'T';
*/
	return ch;
}

#endif

