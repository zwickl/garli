// GARLI version 2.0 source code
// Copyright 2005-2011 Derrick J. Zwickl
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

#include <string>
#include <cstring>
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

class SitePattern{
public:
	int count;
	int origCount;
	int numStates;
	int constStates;
	static int numTax;
	static int maxNumStates;
	vector<unsigned char> stateVec;
	vector<int> siteNumbers;
	enum patternType{
		MISSING = 1,
		CONSTANT = 2,
		UNINFORM_VARIABLE = 3,
		INFORMATIVE = 4
		}type;

	SitePattern(){Reset();}
	SitePattern(const SitePattern &rhs){
		Reset();
		siteNumbers = rhs.siteNumbers;
		stateVec = rhs.stateVec;
		count = rhs.count;
		origCount = rhs.origCount;
		numStates = rhs.numStates;
		constStates = rhs.constStates;
		}
	~SitePattern(){
		stateVec.clear();
		siteNumbers.clear();
		}
	void Reset(){
		count = origCount = numStates = constStates = -1;
		stateVec.clear();
		siteNumbers.clear();
//		if(numTax > 0)
//			stateVec.reserve(numTax);
		}
	static void SetStatics(int nt, int ns){
		numTax = nt;
		maxNumStates = ns;
		}

	//bool PatternLessThan(const SitePattern &lhs, const SitePattern &rhs) const;
	bool operator==(const SitePattern &rhs) const;
	bool operator<(const SitePattern &rhs) const;
	void AddChar(const unsigned char c){
		stateVec.push_back(c);
		}
	void SetCount(int c) {
		count = origCount = c;
		}
	int CalcPatternTypeAndNumStates(vector<unsigned int> &stateCounts);
	int MinScore(set<unsigned char> patt, int bound, unsigned char bits=15, int prevSc=0) const;
	};

class PatternManager{
	friend class DataMatrix;

	int numTax;
	int maxNumStates;
	int numUniquePats;				//NChar

	int totNumChars;				//gapsIncludedNChar
	int numNonMissingChars;			//totalNChar
	int numMissingChars;			//nMissing
	int numConstantChars;			//nConstant
	int numInformativeChars;		//nInformative
	int numUninformVariableChars;	//nVarUninform
	
	int lastConstant;				//lastConstant
	bool compressed;				//dense
	list<SitePattern> patterns;
	list<SitePattern> uniquePatterns;
	vector<int> constStates;

	~PatternManager(){
		patterns.clear();
		uniquePatterns.clear();
		constStates.clear();
		}
	virtual void NewCollapse();
	virtual void NewPack();
	virtual void NewSort();
	virtual void NewDetermineConstantSites();

public:
	void Initialize(int nt, int max){
		Reset();
		numTax = nt;
		maxNumStates = max;
		SitePattern::maxNumStates = max;
		SitePattern::numTax = nt;
		}
	void Reset(){
		numTax = maxNumStates = totNumChars = numNonMissingChars = numUniquePats = numMissingChars = numConstantChars = numInformativeChars = lastConstant = numUninformVariableChars = 0;
		compressed = false;
		patterns.clear();
		uniquePatterns.clear();
		constStates.clear();
		}
	void AddPattern(const SitePattern &add){
		patterns.push_back(add);
		}
	//these are named along the lines of the old DataMatrix members
	int NChar() const {
		if(uniquePatterns.empty())
			return -1;
		else
			return uniquePatterns.size();
		}
	void ProcessPatterns();
	void CalcPatternTypesAndNumStates();
	//funcs for getting info back out of the patman into the datamatrix object
	void FillNumberVector(vector<int> &nums) const;
	void FillTaxaXCharMatrix(unsigned char **mat) const;
	void FillNumStatesVector(vector<int> &ns) const;
	void FillCountVector(vector<int> &counts) const;
	void FillConstStatesVector(vector<int> &cs) const;
	void FillIntegerValues(int &nMissing, int &nConstant, int &nVarUninform, int &nInformative, int &lastConstant, int &gapsIncludedNChar, int &totNChar, int &NChar) const;
	};

// Note: the class below has pure virtual member functions
class DataMatrix{
protected:
	int		nTax;
	int		nTaxAllocated;//allocate more than nTax to allow for the addition of dummy taxa
							//this will only be used during allocating and deallocation
							//if a dummy taxon is created then nTax will be incremented
	int		nChar;     //after compression
	int		totalNChar; //columns in matrix, with all missing columns removed
	int 	gapsIncludedNChar; //the actual number of columns in the data matrix read in
								//only used when outputting something relative to input alignment
	int		dense;		//whether the data has been sorted and identical patterns combined
	int		nonZeroCharCount;	//this is the number of character patterns that have non-zero
								//counts after bootstrap resampling.  Zero count characters can
								//be avoided in the conditional likelihood calcs, but how this
								//is done varies depending on the context
	
	unsigned char**         matrix;
	PatternManager patman;

	int*		count;
	int*		origCounts;
	//maping of chars to columns. indeces are original char numbers, values are the packed column representing that char
	//both start at 0, so offset upon output
	int*		number; 
	/*in the partitioned context number maps the columns of the original partition subset to the columns of the
	compressed matrix.  So, number[j] is the column of the packed matrix that represents column j of the 
	partition subset.  So, this may have no relationship to the original data matrix before the subsets were
	even made.  origDataNumber then maps the columns of the uncompressed subset to the original full datamatrix.
	Thus, number[j] is the compressed column that represents uncompressed subset column j (many-to-one mapping)
	origDataNumber[j] is the column of the orignal matrix that corresponds to uncompressed subset column j (one-to-one mapping)
	example (zero offset): partition by codon position, so sub1 = {0, 3, 6, ...}, sub2 = {1, 4, 7, ...} and sub3 = {2, 5, 8, ...}
	each subset is its own datamatrix object, with its own number and origDataNumber arrays.
	so, sub1->number[0] is the column of the compressed sub1 matrix that represents the first column of sub1
	    (same for sub2 and sub3)
	    sub1->number[1] is the column of the compressed sub2 matrix that represents the second column of sub2
	    (same for sub2 and sub3)
		sub1->origDataNumber[0] = 0
		sub2->origDataNumber[0] = 1
		sub1->origDataNumber[1] = 3
		sub2->origDataNumber[1] = 4
		etc.
	the values in number must the shuffled around as the matrix is compressed
	the values in origDataNumber are set when SetMatrix is called, and don't change thereafter
	*/
	int*		origDataNumber;

	//These are new correlates to the old dynamicaly allocated arrays.  They will be filled from
	//the pattern manager.
	vector<int> newNumber;
	vector<int> newNumStates;
	vector<int> newCount;
	vector<int> newOrigCounts;
	vector<int> newConstStates;
	vector<string> newTaxonLabel;

	char**          taxonLabel;
	int		nMissing;
	int		nConstant;
	int		nInformative;
	int		nVarUninform;
	int 	lastConstant;//DJZ
	int 	*constStates;//the state (or states) that a constant site contains
	unsigned char fullyAmbigChar;
	unsigned numConditioningPatterns;

	protected:
		int*	numStates;
		int     maxNumStates;
		bool	useDefaultWeightsets;
		string	wtsetName;
		bool	usePatternManager;

	protected:
		char	info[80];
		virtual void	SwapCharacters( int i, int j );
		virtual int	ComparePatterns( const int i, const int j ) const;
		void	BSort( int byCounts = 0 );
		void	DebugSaveQSortState( int top, int bottom, int ii, int jj, int xx, const char* title );
		void	QSort( int top, int bottom );
		void	ReplaceTaxonLabel( int i, const char* s );

	public:
		enum {
			PT_MISSING		= 0x0000,
			PT_CONSTANT		= 0x0001,
			PT_INFORMATIVE		= 0x0002,
			PT_VARIABLE		= 0x0004
		};

	public:
		DataMatrix() : dense(0), nTax(0), nChar(0), matrix(0), count(0),
			number(0), taxonLabel(0), numStates(0),
			nMissing(0), nConstant(0), nInformative(0), nVarUninform(0),
			lastConstant(-1), constStates(0), origCounts(0),
			fullyAmbigChar(15), useDefaultWeightsets(true), usePatternManager(false),
			nTaxAllocated(0), origDataNumber(0), numConditioningPatterns(0)
			{ memset( info, 0x00, 80 ); }
		DataMatrix( int ntax, int nchar )
			: nTax(ntax), nChar(nchar), dense(0), matrix(0), count(0),
			number(0), taxonLabel(0), numStates(0),
			nMissing(0), nConstant(0), nInformative(0), nVarUninform(0),
			lastConstant(-1), constStates(0), origCounts(0),
			fullyAmbigChar(15), useDefaultWeightsets(true), usePatternManager(false),
			nTaxAllocated(0), origDataNumber(0), numConditioningPatterns(0)
			{ memset( info, 0x00, 80 ); NewMatrix(ntax, nchar); }
		virtual ~DataMatrix();

		// pure virtual functions - must override in derived class
		virtual unsigned char	CharToDatum( char ch )	const			= 0;
		virtual unsigned char CharToBitwiseRepresentation( char ch ) const= 0;
		virtual char	DatumToChar( unsigned char d )	const           = 0;
		virtual unsigned char	FirstState()		    const           = 0;
		virtual unsigned char	LastState()		        const           = 0;
		virtual void CalcEmpiricalFreqs() = 0;
//		virtual FLOAT_TYPE	Freq( unsigned char, int = 0)	                        = 0;

		// virtual functions - can override in derived class
		virtual FLOAT_TYPE	TransitionProb( int /*i*/, int /*j*/
			, int /*site*/, FLOAT_TYPE /*brlen*/) { return 0.0; }
		virtual int NumStates(int j) const
			{ return ( numStates && (j < nChar) ? numStates[j] : 0 ); }

		void SetUsePatternManager(bool tf) {usePatternManager = tf;}
		bool GetUsePatternManager() const {return usePatternManager;}
		void ProcessPatterns();
		void OutputDataSummary() const;

		void GetDataFromPatternManager();
		// functions for getting the data in and out
		int GetToken( istream& in, char* tokenbuf, int maxlen, bool acceptComments=true );
		int GetToken( FILE *in, char* tokenbuf, int maxlen);
		int ReadPhylip( const char* filename);
		int ReadFasta( const char* filename);
		int Save( const char* filename, char* newfname = 0, char* nxsfname = 0 );

		char*	DataType() { return info; }
		int     unsigned charToInt( unsigned char d ) const { return (int)d; }

		int NTax() const { return nTax; }
		void SetNTax(int ntax) { nTax = ntax; }

		virtual int NChar() const { return nChar; }
		int TotalNChar() const { return totalNChar; }
		int GapsIncludedNChar() const { return gapsIncludedNChar; }
		void SetNChar(int nchar) { nChar = nchar; }
		unsigned NumConditioningPatterns() const{return numConditioningPatterns;}

		int BootstrappedNChar() {return nonZeroCharCount;} 
		void Flush() { NewMatrix( 0, 0 ); }
		int Dense() const { return dense; }
		
		//argument here is column number from uncompressed subset
		//return val is compressed pattern representing that column
		int Number(int j) const{
			if(newNumber.size() > 0)
				return newNumber[j];
			return ( number && (j < gapsIncludedNChar) ? number[j] : 0 );
			}

		//argument here is column number from uncompressed subset
		//return val is column from original full matrix before partitioning
		int OrigDataNumber(int j) const
			{ return ( origDataNumber && (j < gapsIncludedNChar) ? origDataNumber[j] : 0 ); }

		virtual int Count(int j) const{ 
			if(newCount.size() > 0)
				return newCount[j];
			return ( count && (j < nChar) ? count[j] : 0 ); 
			}
		virtual int CountByOrigIndex(int j) const{ 
			if(newCount.size() > 0)
				if(newNumber.size() > 0){
					assert(newCount.size() > j);
					return newCount[newNumber[j]];
					}
			return ( count && (j < nChar) ? count[number[j]] : 0 ); 
			}
		virtual const int *GetCounts() const {
			if(newCount.size() > 0)
				return &(newCount[0]);
			return count;
			}
		const int *GetConstStates() const {
			if(newConstStates.size() > 0)
				return &(newConstStates[0]);
			return constStates;
			}
		void SetCount(int j, int c){
			if(newCount.size() > 0){
				assert(newCount.size() > j);
				newCount[j] = c;
				}
			else
				if( count && (j < nChar) ) 
					count[j] = c; 
			}
		void SetNumStates(int j, int c){ 
			if( numStates && (j < nChar) ) numStates[j] = c;
			}
		const char* TaxonLabel(int i) const{
			return ( taxonLabel && (i < nTax) ? taxonLabel[i] : 0 );
			}
		void SetTaxonLabel(int i, const char* s);
		
		int TaxonNameToNumber(const NxsString &name) const;

		void CopyNamesFromOtherMatrix(const DataMatrix *dat){
			assert(taxonLabel);
			for(int t=0;t<nTax;t++) SetTaxonLabel(t, dat->TaxonLabel(t));
			}

		void BeginNexusTreesBlock(ofstream &treeout) const;
		void BeginNexusTreesBlock(string &trans) const;

		virtual void CreateMatrixFromNCL(const NxsCharactersBlock *, NxsUnsignedSet &charset) = 0;
		
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
		virtual void SetMatrix( int i, int j, unsigned char c){
			if(matrix && (i < nTax) && (j < nChar))
				matrix[i][j] = c;
			}
		void SetOriginalDataNumber(const int subsetMatColumn, const int origMatColumn){
			origDataNumber[subsetMatColumn] = origMatColumn;			
			}

		int MatrixExists() const { return ( matrix && nTax>0 && nChar>0 ? 1 : 0 ); }
		int NMissing() const { return nMissing; }
		int NConstant() const { return nConstant; }
		int LastConstant() const {return lastConstant;}
		int NInformative() const { return nInformative; }
		int NVarUninform() const { return nVarUninform; }

		DataMatrix& operator =(const DataMatrix&);

		void Sort( int byCounts = 0 ){
			byCounts;
			QSort( 0, NChar()-1 );
			}
		virtual int PatternType( int , unsigned int *) const;	// returns PT_XXXX constant indicating type of pattern
		void Summarize();       // fills in nConstant, nInformative, and nAutapomorphic data members
		virtual void Collapse();
		void EliminateAdjacentIdenticalColumns();
		virtual void Pack();
		void NewMatrix(int nt, int nc);	// flushes old matrix, creates new one
		void ResizeCharacterNumberDependentVariables(int nCh);
		int PositionOf( char* s ) const; // returns pos (0..nTax-1) of taxon named s
		void DumpCounts( const char* s );
		void WriteCollapsedData();  //DZ
		void SaveNexus(const char* filename, int iosFlags /* = 0 */); //DZ
		virtual void DetermineConstantSites();
		void ExplicitDestructor();  // cjb - totally clear the DataMatrix and revert it to its original state as if it was just constructed
		void CheckForIdenticalTaxonNames();
		bool DidUseDefaultWeightsets() const {return (wtsetName.length() > 0);}
		string WeightsetName() const { return wtsetName;}

		//for determining parsimony informative chars
		int MinScore(set<unsigned char> patt, int bound, unsigned char bits=15, int sc=0) const;
		void GetStringOfOrigDataColumns(string &str) const;

	public:	
      void ReserveOriginalCounts(){
      		if(origCounts==NULL) origCounts=new int[nChar];
      		for(int i=0;i<nChar;i++){
				if(newCount.size() > 0){
					assert(newCount.size() > i);
					origCounts[i] = newCount[i];
					}
				else
					origCounts[i] = count[i];
      			}
      		}
      void RestoreOriginalCounts(){
			if(origCounts==NULL) return;
      		for(int i=0;i<nChar;i++){
				if(newCount.size() > 0){
					assert(newCount.size() > i);
					newCount[i] = origCounts[i];
					}
				else
      				count[i] = origCounts[i];
      			}
      		}
      void Reweight(FLOAT_TYPE prob);
      virtual int BootstrapReweight(int seedToUse, FLOAT_TYPE resampleProportion);
	  void CountMissingCharsByColumn(vector<int> &vec);
	  void MakeWeightSetString(NxsCharactersBlock &charblock, string &wtstring, string name);
      void MakeWeightSetString(std::string &wtstring, string name);
};

#endif

