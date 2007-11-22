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

#ifndef _SEQUENCE_DATA_
#define _SEQUENCE_DATA_

#include <vector>
using namespace std;

#include "defs.h"
#include "datamatr.h"
//#include "model.h"

class SequenceData : public DataMatrix{
public:
	SequenceData() : DataMatrix()
		{ maxNumStates=4; strcpy( info, "DNA" ); empStateFreqs=NULL;}
	SequenceData( int ntax, int nchar ) : DataMatrix( ntax, nchar )
		{ maxNumStates=4; strcpy( info, "DNA" ); empStateFreqs=NULL;}
	virtual ~SequenceData() {
		if(empStateFreqs != NULL) delete []empStateFreqs;
		}

protected:
	FLOAT_TYPE *empStateFreqs;
	// overrides of base class's virtual fuctions
	virtual unsigned char CharToDatum( char ch ) = 0;
	virtual unsigned char CharToBitwiseRepresentation( char ch );
	virtual char	DatumToChar( unsigned char d );
	virtual unsigned char	FirstState() { return 0; }
	virtual unsigned char	LastState() { return 3; }
	virtual int	NumStates(int) const { return 4; }

public:
	virtual void CreateMatrixFromNCL(GarliReader &reader) = 0;
	virtual void CalcEmpiricalFreqs() = 0;
	virtual void GetEmpiricalFreqs(FLOAT_TYPE *f) const{
		for(int i=0;i<maxNumStates;i++) f[i]=empStateFreqs[i];
		}
	};

inline unsigned char SequenceData::CharToDatum( char ch )
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

inline unsigned char SequenceData::CharToBitwiseRepresentation( char ch )
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
       	default  : throw ErrorException("Unknown nucleotide %c!", ch);
		}
	return datum;
}

inline char SequenceData::DatumToChar( unsigned char d )
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

class NucleotideData : public SequenceData{

	vector<char*> ambigStrings;
#ifdef OPEN_MP
	vector<unsigned*> ambigToCharMap;
#endif

public:		
	NucleotideData() : SequenceData() {}
	NucleotideData( int ntax, int nchar ) : SequenceData( ntax, nchar ) {}
	~NucleotideData() {
		for(vector<char*>::iterator delit=ambigStrings.begin();delit!=ambigStrings.end();delit++)
			delete [](*delit);
#ifdef OPEN_MP
		for(vector<unsigned*>::iterator delit=ambigToCharMap.begin();delit!=ambigToCharMap.end();delit++)
			delete [](*delit);
#endif
		}

	unsigned char CharToDatum(char d);
	void CalcEmpiricalFreqs();
	void CreateMatrixFromNCL(GarliReader &reader);
	void MakeAmbigStrings();
	char *GetAmbigString(int i) const{
		return ambigStrings[i];
		}
#ifdef OPEN_MP
	unsigned *GetAmbigToCharMap(int i) const{
		return ambigToCharMap[i];
		}
#endif
	};

class GeneticCode{
	//mapping from codon number (ordered AAA, AAC, AAG, AAT, ACA, etc) to 
	//amino acid number (0-19). Stop codons are 20.
	int codonTable[64];
	int map64to61[64];

	public:
	GeneticCode(){
		SetStandardCode();
		}
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

class CodonData : public SequenceData {

	GeneticCode code;
	//these are for use in the F1x4 or F3x4 methods of calculating the state freqs
	FLOAT_TYPE empBaseFreqsPos1[4];
	FLOAT_TYPE empBaseFreqsPos2[4];
	FLOAT_TYPE empBaseFreqsPos3[4];

	FLOAT_TYPE empBaseFreqsAllPos[4];
	int empType; //codon table = 0
				//F1x4 = 1
				//F3x4 = 2

public:
	CodonData() : SequenceData(){
		maxNumStates = 61;
		code.SetStandardCode();
		empType = 0;
		}

	CodonData(const NucleotideData *dat, int genCode) : SequenceData(){
		assert(dat->Dense() == false);
		maxNumStates = 61;
		if(genCode == 0) code.SetStandardCode();
		else code.SetVertMitoCode();
		FillCodonMatrixFromDNA(dat);
		CopyNamesFromOtherMatrix(dat);
		empType = 0;
		}

	~CodonData(){}

	void FillCodonMatrixFromDNA(const NucleotideData *);
	unsigned char CharToDatum(char c) {
		//this shouldn't be getting called, as it makes no sense for codon data
		assert(0);
		return 0;
		}
	void CreateMatrixFromNCL(GarliReader &reader){

		//this also should not be getting called.  The codon matrix
		//is created from a DNA matrix that has been read in, possibly
		//by the NCL
		assert(0);
		}
	void SetEmpType(int t) {empType = t;}
	void CalcEmpiricalFreqs();
	void CalcF1x4Freqs();
	void CalcF3x4Freqs();
	//int ComparePatterns( const int i, const int j ) const;
	void SetVertMitoCode() {code.SetVertMitoCode();}
};


class AminoacidData : public SequenceData{

public:
	AminoacidData() : SequenceData(){
		maxNumStates = 20;
		}

	AminoacidData(const NucleotideData *dat, int genCode) : SequenceData(){
		maxNumStates = 20;
		GeneticCode c;
		if(genCode == 0) c.SetStandardCode();
		else c.SetVertMitoCode();
		FillAminoacidMatrixFromDNA(dat, &c);
		CopyNamesFromOtherMatrix(dat);
		}
	void FillAminoacidMatrixFromDNA(const NucleotideData *dat, GeneticCode *code);
	void CalcEmpiricalFreqs();
	unsigned char CharToDatum(char d);
	void CreateMatrixFromNCL(GarliReader &reader);

	};

#endif
