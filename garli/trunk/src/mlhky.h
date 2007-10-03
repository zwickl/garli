// GARLI version 0.96b4 source code
// Copyright  2005-2006 by Derrick J. Zwickl
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//  Derrick Zwickl
//	National Evolutionary Synthesis Center
//	2024 W. Main Street, Suite A200
//	Durham, NC 27705
//  email: zwickl@nescent.org
//
//	NOTE: Portions of this source adapted from GAML source, written by Paul O. Lewis

#ifndef __MLHKY_H
#define __MLHKY_H

#include <vector>
using namespace std;

#include "datamatr.h"
#include "defs.h"

class HKYData : public DNAData
{

	public:
		//DZ 9-2-04
		vector<char*> ambigStrings;
#ifdef OPEN_MP
		vector<unsigned*> ambigToCharMap;
#endif
		FLOAT_TYPE *empStateFreqs;
		
		HKYData() : DNAData() {empStateFreqs=NULL;}
		HKYData( int ntax, int nchar ) : DNAData( ntax, nchar ) {empStateFreqs=NULL;}
		~HKYData() {
			for(vector<char*>::iterator delit=ambigStrings.begin();delit!=ambigStrings.end();delit++)
				delete [](*delit);
			if(empStateFreqs != NULL) delete []empStateFreqs;
			}
		virtual void CalcEmpiricalFreqs();
		void CalcEmpiricalAAFreqs();
		virtual void GetEmpiricalFreqs(FLOAT_TYPE *f) const{
			for(int i=0;i<maxNumStates;i++) f[i]=empStateFreqs[i];
			}
	
		void MakeAmbigStrings();
		char *GetAmbigString(int i) const{
			return ambigStrings[i];
			}
#ifdef OPEN_MP
		unsigned *GetAmbigToCharMap(int i) const{
			return ambigToCharMap[i];
			}
#endif
		// overrides of base class's virtual fuctions
		FLOAT_TYPE Freq( unsigned char d, int = 0);
};

class GeneticCode{
	//mapping from codon number (ordered AAA, AAC, AAG, AAT, ACA, etc) to 
	//amino acid number (0-19). Stop codons are 20.
	int codonTable[64];
	int map64to61[64];

	public:
	GeneticCode(){};
	void SetStandardCode(){
		codonTable[ 0 ]= 8;
		codonTable[ 1 ]= 11;
		codonTable[ 2 ]= 8;
		codonTable[ 3 ]= 11;
		codonTable[ 4 ]= 16;
		codonTable[ 5 ]= 16;
		codonTable[ 6 ]= 16;
		codonTable[ 7 ]= 16;
		codonTable[ 8 ]= 14;
		codonTable[ 9 ]= 15;
		codonTable[ 10 ]= 14;
		codonTable[ 11 ]= 15;
		codonTable[ 12 ]= 7;
		codonTable[ 13 ]= 7;
		codonTable[ 14 ]= 10;
		codonTable[ 15 ]= 7;
		codonTable[ 16 ]= 13;
		codonTable[ 17 ]= 6;
		codonTable[ 18 ]= 13;
		codonTable[ 19 ]= 6;
		codonTable[ 20 ]= 12;
		codonTable[ 21 ]= 12;
		codonTable[ 22 ]= 12;
		codonTable[ 23 ]= 12;
		codonTable[ 24 ]= 14;
		codonTable[ 25 ]= 14;
		codonTable[ 26 ]= 14;
		codonTable[ 27 ]= 14;
		codonTable[ 28 ]= 9;
		codonTable[ 29 ]= 9;
		codonTable[ 30 ]= 9;
		codonTable[ 31 ]= 9;
		codonTable[ 32 ]= 3;
		codonTable[ 33 ]= 2;
		codonTable[ 34 ]= 3;
		codonTable[ 35 ]= 2;
		codonTable[ 36 ]= 0;
		codonTable[ 37 ]= 0;
		codonTable[ 38 ]= 0;
		codonTable[ 39 ]= 0;
		codonTable[ 40 ]= 5;
		codonTable[ 41 ]= 5;
		codonTable[ 42 ]= 5;
		codonTable[ 43 ]= 5;
		codonTable[ 44 ]= 17;
		codonTable[ 45 ]= 17;
		codonTable[ 46 ]= 17;
		codonTable[ 47 ]= 17;
		codonTable[ 48 ]= 20;
		codonTable[ 49 ]= 19;
		codonTable[ 50 ]= 20;
		codonTable[ 51 ]= 19;
		codonTable[ 52 ]= 15;
		codonTable[ 53 ]= 15;
		codonTable[ 54 ]= 15;
		codonTable[ 55 ]= 15;
		codonTable[ 56 ]= 20;
		codonTable[ 57 ]= 1;
		codonTable[ 58 ]= 18;
		codonTable[ 59 ]= 1;
		codonTable[ 60 ]= 9;
		codonTable[ 61 ]= 4;
		codonTable[ 62 ]= 9;
		codonTable[ 63 ]= 4;
		
		map64to61[0]=0;
		map64to61[1]=1;
		map64to61[2]=2;
		map64to61[3]=3;
		map64to61[4]=4;
		map64to61[5]=5;
		map64to61[6]=6;
		map64to61[7]=7;
		map64to61[8]=8;
		map64to61[9]=9;
		map64to61[10]=10;
		map64to61[11]=11;
		map64to61[12]=12;
		map64to61[13]=13;
		map64to61[14]=14;
		map64to61[15]=15;
		map64to61[16]=16;
		map64to61[17]=17;
		map64to61[18]=18;
		map64to61[19]=19;
		map64to61[20]=20;
		map64to61[21]=21;
		map64to61[22]=22;
		map64to61[23]=23;
		map64to61[24]=24;
		map64to61[25]=25;
		map64to61[26]=26;
		map64to61[27]=27;
		map64to61[28]=28;
		map64to61[29]=29;
		map64to61[30]=30;
		map64to61[31]=31;
		map64to61[32]=32;
		map64to61[33]=33;
		map64to61[34]=34;
		map64to61[35]=35;
		map64to61[36]=36;
		map64to61[37]=37;
		map64to61[38]=38;
		map64to61[39]=39;
		map64to61[40]=40;
		map64to61[41]=41;
		map64to61[42]=42;
		map64to61[43]=43;
		map64to61[44]=44;
		map64to61[45]=45;
		map64to61[46]=46;
		map64to61[47]=47;
		map64to61[48]=-1;
		map64to61[49]=48;
		map64to61[50]=-1;
		map64to61[51]=49;
		map64to61[52]=50;
		map64to61[53]=51;
		map64to61[54]=52;
		map64to61[55]=53;
		map64to61[56]=-1;
		map64to61[57]=54;
		map64to61[58]=55;
		map64to61[59]=56;
		map64to61[60]=57;
		map64to61[61]=58;
		map64to61[62]=59;
		map64to61[63]=60;
		}
		
	void SetVertMitoCode(){
		SetStandardCode();
		codonTable[56] = 18;
		codonTable[8] = 20;
		codonTable[47] = 20;

		map64to61[0]=0;
		map64to61[1]=1;
		map64to61[2]=2;
		map64to61[3]=3;
		map64to61[4]=4;
		map64to61[5]=5;
		map64to61[6]=6;
		map64to61[7]=7;
		map64to61[8]=-1;
		map64to61[9]=8;
		map64to61[10]=9;
		map64to61[11]=10;
		map64to61[12]=11;
		map64to61[13]=12;
		map64to61[14]=13;
		map64to61[15]=14;
		map64to61[16]=15;
		map64to61[17]=16;
		map64to61[18]=17;
		map64to61[19]=18;
		map64to61[20]=19;
		map64to61[21]=20;
		map64to61[22]=21;
		map64to61[23]=22;
		map64to61[24]=23;
		map64to61[25]=24;
		map64to61[26]=25;
		map64to61[27]=26;
		map64to61[28]=27;
		map64to61[29]=28;
		map64to61[30]=29;
		map64to61[31]=30;
		map64to61[32]=31;
		map64to61[33]=32;
		map64to61[34]=33;
		map64to61[35]=34;
		map64to61[36]=35;
		map64to61[37]=36;
		map64to61[38]=37;
		map64to61[39]=38;
		map64to61[40]=39;
		map64to61[41]=40;
		map64to61[42]=41;
		map64to61[43]=42;
		map64to61[44]=43;
		map64to61[45]=44;
		map64to61[46]=45;
		map64to61[47]=-1;
		map64to61[48]=46;
		map64to61[49]=47;
		map64to61[50]=-1;
		map64to61[51]=48;
		map64to61[52]=49;
		map64to61[53]=50;
		map64to61[54]=51;
		map64to61[55]=52;
		map64to61[56]=53;
		map64to61[57]=54;
		map64to61[58]=55;
		map64to61[59]=56;
		map64to61[60]=57;
		map64to61[61]=58;
		map64to61[62]=59;
		map64to61[63]=60;
		}
		
	int CodonLookup(int i){
		assert(i >= 0 && i < 64);
		return codonTable[i];
		}
	int Map64stateTo61state(int i){
		assert(i >= 0 && i < 64);
		return map64to61[i];
		}
	};

class CodonData : public HKYData {
	char **codonMatrix;
	int *codonCount;
	int *numCodonStates;
	int *codonNumber;
	int nCodonChar;
	int gapsIncludedNCodonChar;
	GeneticCode code;
	
	FLOAT_TYPE *empBaseFreqPos1;
	FLOAT_TYPE *empBaseFreqPos2;
	FLOAT_TYPE *empBaseFreqPos3;

public:
	CodonData() : HKYData(), codonMatrix(NULL), codonCount(NULL), numCodonStates(NULL), codonNumber(NULL),
		empBaseFreqPos1(NULL), empBaseFreqPos2(NULL), empBaseFreqPos3(NULL) {maxNumStates = 61;code.SetStandardCode();}
	~CodonData(){
		if( codonCount ) MEM_DELETE_ARRAY(codonCount); // count is of length nChar
		if( codonMatrix ) {
			int j;
			for( j = 0; j < nTax; j++ )
				MEM_DELETE_ARRAY(codonMatrix[j]); // matrix[j] is of length nChar
			MEM_DELETE_ARRAY(codonMatrix); // matrix is of length nTax
			}
		}

	void SetAminoAcid(){maxNumStates = 20;}
	void FillCodonMatrix(bool);
	void NewCodonMatrix( int taxa, int sites );
	int NChar() const {
		if(codonMatrix) return nCodonChar;
		else return nChar;
		}
	void Collapse();
	void Pack();
	void SetMatrix( int i, int j, unsigned char c ) {
		if(codonMatrix){
			assert(c > -1);
			if((i < nTax) && (j < nCodonChar) ) codonMatrix[i][j] = c; 
			}
		else if((i < nTax) && (j < nChar) ) matrix[i][j] = c; 
		}
	void SwapCharacters( int i, int j );
	int NNucChar() const{
		return nChar;
		}
	const int *GetCounts() const{
		if(codonMatrix) return codonCount;
		else return count;
		}
	int Count(int j) const{
		if(codonMatrix) return ( codonCount && (j < nCodonChar) ? codonCount[j] : 0 ); 
		else return ( count && (j < nChar) ? count[j] : 0 ); 
		}
	void CalcEmpiricalFreqs();
	void CalcF1x4Freqs();
	void CalcF3x4Freqs();
	int ComparePatterns( const int i, const int j ) const;
	int PatternType( int k , int *c, unsigned char *s) const;
	void GetEmpiricalFreqs(FLOAT_TYPE *f) const{
		for(int i=0;i<maxNumStates;i++) f[i]=empStateFreqs[i];
		}
	const char *GetCodonRow(int i)const{
		assert( codonMatrix );
		assert( i >= 0 );
		assert( i < NTax() );
		return codonMatrix[i];
		}
	unsigned char Matrix( int i, int j ) const {
		assert( i >= 0 );
		assert( i < nTax );
		assert( j >= 0 );
		if( codonMatrix ){
			assert( j < nCodonChar );
			return (unsigned char) codonMatrix[i][j];
			}
		else{
			assert( j < nChar );
			return (unsigned char) matrix[i][j];
			}
		}
	unsigned char NucMatrix( int i, int j ) const {
		assert( i >= 0 );
		assert( i < nTax );
		assert( j >= 0 );
		assert( j < nChar );
		return matrix[i][j];
		}
	void SetVertMitoCode() {code.SetVertMitoCode();}
};

#endif
