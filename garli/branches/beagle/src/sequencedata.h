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
	virtual unsigned char CharToBitwiseRepresentation( char ch ) const;
	virtual char	DatumToChar( unsigned char d ) const;
	virtual unsigned char	FirstState() const { return 0; }
	virtual unsigned char	LastState() const { return 3; }
	virtual int	NumStates(int) const { return 4; }

public:
	virtual void CreateMatrixFromNCL(const NxsCharactersBlock *) = 0;
	virtual void CreateMatrixFromNCL(const NxsCharactersBlock *, NxsUnsignedSet &charset) = 0;
	virtual void CalcEmpiricalFreqs() = 0;
	virtual void GetEmpiricalFreqs(FLOAT_TYPE *f) const{
		assert(empStateFreqs);
		for(int i=0;i<maxNumStates;i++) f[i]=empStateFreqs[i];
		}
	};

inline unsigned char SequenceData::CharToDatum( char ch ){
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

inline unsigned char SequenceData::CharToBitwiseRepresentation( char ch ) const{
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

inline char SequenceData::DatumToChar( unsigned char d ) const
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
	NucleotideData() : SequenceData() {fullyAmbigChar = 15;}
	NucleotideData( int ntax, int nchar ) : SequenceData( ntax, nchar ) {fullyAmbigChar = 15;}
	~NucleotideData() {
		for(vector<char*>::iterator delit=ambigStrings.begin();delit!=ambigStrings.end();delit++)
			delete [](*delit);
#ifdef OPEN_MP
		for(vector<unsigned*>::iterator delit=ambigToCharMap.begin();delit!=ambigToCharMap.end();delit++)
			delete [](*delit);
#endif
		}

	unsigned char CharToDatum(char d) ;
	void CalcEmpiricalFreqs();
	void CreateMatrixFromNCL(const NxsCharactersBlock *);
	void CreateMatrixFromNCL(const NxsCharactersBlock *charblock, NxsUnsignedSet &charset);
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
	int map64toNonStops[64];
	vector<int> stops;
	//this holds the correspondence between the state indeces and actual codons
	//for display purposes.  Stops are removed and thus any mapIndexToCodonDisplay[index]
	//gives the codon for that index
	vector<string> mapIndexToCodonDisplay;

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
		
		map64toNonStops[0]=0;
		map64toNonStops[1]=1;
		map64toNonStops[2]=2;
		map64toNonStops[3]=3;
		map64toNonStops[4]=4;
		map64toNonStops[5]=5;
		map64toNonStops[6]=6;
		map64toNonStops[7]=7;
		map64toNonStops[8]=8;
		map64toNonStops[9]=9;
		map64toNonStops[10]=10;
		map64toNonStops[11]=11;
		map64toNonStops[12]=12;
		map64toNonStops[13]=13;
		map64toNonStops[14]=14;
		map64toNonStops[15]=15;
		map64toNonStops[16]=16;
		map64toNonStops[17]=17;
		map64toNonStops[18]=18;
		map64toNonStops[19]=19;
		map64toNonStops[20]=20;
		map64toNonStops[21]=21;
		map64toNonStops[22]=22;
		map64toNonStops[23]=23;
		map64toNonStops[24]=24;
		map64toNonStops[25]=25;
		map64toNonStops[26]=26;
		map64toNonStops[27]=27;
		map64toNonStops[28]=28;
		map64toNonStops[29]=29;
		map64toNonStops[30]=30;
		map64toNonStops[31]=31;
		map64toNonStops[32]=32;
		map64toNonStops[33]=33;
		map64toNonStops[34]=34;
		map64toNonStops[35]=35;
		map64toNonStops[36]=36;
		map64toNonStops[37]=37;
		map64toNonStops[38]=38;
		map64toNonStops[39]=39;
		map64toNonStops[40]=40;
		map64toNonStops[41]=41;
		map64toNonStops[42]=42;
		map64toNonStops[43]=43;
		map64toNonStops[44]=44;
		map64toNonStops[45]=45;
		map64toNonStops[46]=46;
		map64toNonStops[47]=47;
		map64toNonStops[48]=-1;
		map64toNonStops[49]=48;
		map64toNonStops[50]=-1;
		map64toNonStops[51]=49;
		map64toNonStops[52]=50;
		map64toNonStops[53]=51;
		map64toNonStops[54]=52;
		map64toNonStops[55]=53;
		map64toNonStops[56]=-1;
		map64toNonStops[57]=54;
		map64toNonStops[58]=55;
		map64toNonStops[59]=56;
		map64toNonStops[60]=57;
		map64toNonStops[61]=58;
		map64toNonStops[62]=59;
		map64toNonStops[63]=60;

		stops.clear();
		stops.push_back(48);
		stops.push_back(50);
		stops.push_back(56);

		FillIndexToCodonDisplayMap();
		}
		
	void SetVertMitoCode(){
		SetStandardCode();
		codonTable[56] = 18; //TGA
		codonTable[8] = 20;  //AGA
		codonTable[10] = 20;  //AGG
		codonTable[12] = 10;  //ATA

		map64toNonStops[0]=0;
		map64toNonStops[1]=1;
		map64toNonStops[2]=2;
		map64toNonStops[3]=3;
		map64toNonStops[4]=4;
		map64toNonStops[5]=5;
		map64toNonStops[6]=6;
		map64toNonStops[7]=7;
		map64toNonStops[8]=-1;
		map64toNonStops[9]=8;
		map64toNonStops[10]=-1;
		map64toNonStops[11]=9;
		map64toNonStops[12]=10;
		map64toNonStops[13]=11;
		map64toNonStops[14]=12;
		map64toNonStops[15]=13;
		map64toNonStops[16]=14;
		map64toNonStops[17]=15;
		map64toNonStops[18]=16;
		map64toNonStops[19]=17;
		map64toNonStops[20]=18;
		map64toNonStops[21]=19;
		map64toNonStops[22]=20;
		map64toNonStops[23]=21;
		map64toNonStops[24]=22;
		map64toNonStops[25]=23;
		map64toNonStops[26]=24;
		map64toNonStops[27]=25;
		map64toNonStops[28]=26;
		map64toNonStops[29]=27;
		map64toNonStops[30]=28;
		map64toNonStops[31]=29;
		map64toNonStops[32]=30;
		map64toNonStops[33]=31;
		map64toNonStops[34]=32;
		map64toNonStops[35]=33;
		map64toNonStops[36]=34;
		map64toNonStops[37]=35;
		map64toNonStops[38]=36;
		map64toNonStops[39]=37;
		map64toNonStops[40]=38;
		map64toNonStops[41]=39;
		map64toNonStops[42]=40;
		map64toNonStops[43]=41;
		map64toNonStops[44]=42;
		map64toNonStops[45]=43;
		map64toNonStops[46]=44;
		map64toNonStops[47]=45;
		map64toNonStops[48]=-1;
		map64toNonStops[49]=46;
		map64toNonStops[50]=-1;
		map64toNonStops[51]=47;
		map64toNonStops[52]=48;
		map64toNonStops[53]=49;
		map64toNonStops[54]=50;
		map64toNonStops[55]=51;
		map64toNonStops[56]=52;
		map64toNonStops[57]=53;
		map64toNonStops[58]=54;
		map64toNonStops[59]=55;
		map64toNonStops[60]=56;
		map64toNonStops[61]=57;
		map64toNonStops[62]=58;
		map64toNonStops[63]=59;

		stops.clear();
		stops.push_back(8);
		stops.push_back(10);
		stops.push_back(48);
		stops.push_back(50);

		FillIndexToCodonDisplayMap();
		}
		
	void SetInvertMitoCode(){
		SetStandardCode();
		codonTable[56] = 18; //TGA
		codonTable[8] = 15;  //AGA
		codonTable[10] = 15;  //AGG
		codonTable[12] = 10;  //ATA

		map64toNonStops[0]=0;
		map64toNonStops[1]=1;
		map64toNonStops[2]=2;
		map64toNonStops[3]=3;
		map64toNonStops[4]=4;
		map64toNonStops[5]=5;
		map64toNonStops[6]=6;
		map64toNonStops[7]=7;
		map64toNonStops[8]=8;
		map64toNonStops[9]=9;
		map64toNonStops[10]=10;
		map64toNonStops[11]=11;
		map64toNonStops[12]=12;
		map64toNonStops[13]=13;
		map64toNonStops[14]=14;
		map64toNonStops[15]=15;
		map64toNonStops[16]=16;
		map64toNonStops[17]=17;
		map64toNonStops[18]=18;
		map64toNonStops[19]=19;
		map64toNonStops[20]=20;
		map64toNonStops[21]=21;
		map64toNonStops[22]=22;
		map64toNonStops[23]=23;
		map64toNonStops[24]=24;
		map64toNonStops[25]=25;
		map64toNonStops[26]=26;
		map64toNonStops[27]=27;
		map64toNonStops[28]=28;
		map64toNonStops[29]=29;
		map64toNonStops[30]=30;
		map64toNonStops[31]=31;
		map64toNonStops[32]=32;
		map64toNonStops[33]=33;
		map64toNonStops[34]=34;
		map64toNonStops[35]=35;
		map64toNonStops[36]=36;
		map64toNonStops[37]=37;
		map64toNonStops[38]=38;
		map64toNonStops[39]=39;
		map64toNonStops[40]=40;
		map64toNonStops[41]=41;
		map64toNonStops[42]=42;
		map64toNonStops[43]=43;
		map64toNonStops[44]=44;
		map64toNonStops[45]=45;
		map64toNonStops[46]=46;
		map64toNonStops[47]=47;
		map64toNonStops[48]=-1;
		map64toNonStops[49]=48;
		map64toNonStops[50]=-1;
		map64toNonStops[51]=49;
		map64toNonStops[52]=50;
		map64toNonStops[53]=51;
		map64toNonStops[54]=52;
		map64toNonStops[55]=53;
		map64toNonStops[56]=54;
		map64toNonStops[57]=55;
		map64toNonStops[58]=56;
		map64toNonStops[59]=57;
		map64toNonStops[60]=58;
		map64toNonStops[61]=59;
		map64toNonStops[62]=60;
		map64toNonStops[63]=61;

		stops.clear();
		stops.push_back(48);
		stops.push_back(50);
		
		FillIndexToCodonDisplayMap();
		}

	int CodonLookup(int i){
		assert(i >= 0 && i < 64);
		return codonTable[i];
		}
	int Map64stateToNonStops(int i){
		assert(i >= 0 && i < 64);
		assert(map64toNonStops[i] != -1);
		return map64toNonStops[i];
		}
	void FillIndexToCodonDisplayMap(){
		//this assumes that the correct genetic code has already been set
		mapIndexToCodonDisplay.clear();
		char nucs[4] = {'A', 'C', 'G', 'T'};
		//char cod[3];
		char *cod = new char[4];
		for(int f = 0;f < 4;f++){
			for(int s = 0;s < 4;s++){
				for(int t = 0;t < 4;t++){
					if(CodonLookup(f * 16 + s * 4 + t) != 20){ 
						sprintf(cod, "%c%c%c\0", nucs[f], nucs[s], nucs[t]);
						mapIndexToCodonDisplay.push_back(cod);
						}
					}
				}
			}
		delete []cod;
		}
	const string LookupCodonDisplayFromIndex(int index) const{
		return mapIndexToCodonDisplay[index];
		}

	unsigned NumStates() const {return (unsigned) mapIndexToCodonDisplay.size();}
	};

class CodonData : public SequenceData {

	GeneticCode code;
	//these are for use in the F1x4 or F3x4 methods of calculating the state freqs
	FLOAT_TYPE empBaseFreqsPos1[4];
	FLOAT_TYPE empBaseFreqsPos2[4];
	FLOAT_TYPE empBaseFreqsPos3[4];

	FLOAT_TYPE empBaseFreqsAllPos[4];
	enum{ NOT_EMPIRICAL	= 0,
		  CODON_TABLE	= 1,
		  F1X4			= 2,
		  F3X4			= 3
		}empType;
//	int empType; //codon table = 0
				//F1x4 = 1
				//F3x4 = 2

public:
	CodonData() : SequenceData(){
		maxNumStates = 61;
		code.SetStandardCode();
		empType = NOT_EMPIRICAL;
		fullyAmbigChar = maxNumStates;
		}

	CodonData(const NucleotideData *dat, int genCode) : SequenceData(){
		assert(dat->Dense() == false);
		if(genCode == 0){
			code.SetStandardCode();
			maxNumStates = 61;
			}
		else if(genCode == 1){
			code.SetVertMitoCode();
			maxNumStates = 60;
			}
		else{
			code.SetInvertMitoCode();
			maxNumStates = 62;
			}
		FillCodonMatrixFromDNA(dat);
		CopyNamesFromOtherMatrix(dat);
		empType = NOT_EMPIRICAL;
		fullyAmbigChar = maxNumStates;
		}

	~CodonData(){}

	void FillCodonMatrixFromDNA(const NucleotideData *);
	unsigned char CharToDatum(char c)  {
		//this shouldn't be getting called, as it makes no sense for codon data
		assert(0);
		return 0;
		}
	void CreateMatrixFromNCL(const NxsCharactersBlock *){
		//this also should not be getting called.  The codon matrix
		//is created from a DNA matrix that has been read in, possibly
		//by the NCL
		assert(0);
		}
	void CreateMatrixFromNCL(const NxsCharactersBlock *, NxsUnsignedSet &charset){
		//this also should not be getting called.  The codon matrix
		//is created from a DNA matrix that has been read in, possibly
		//by the NCL
		assert(0);
		}
	GeneticCode* GetCode() {return &code;}
	//void SetEmpType(int t) {empType = t;}
	void SetF1X4Freqs(){empType = F1X4;}
	void SetF3X4Freqs(){empType = F3X4;}
	void SetCodonTableFreqs(){empType = CODON_TABLE;}
	void CalcEmpiricalFreqs();
	void CalcF1x4Freqs();
	void CalcF3x4Freqs();
	void BaseFreqXPositionReport();
	//int ComparePatterns( const int i, const int j ) const;
	void SetVertMitoCode() {code.SetVertMitoCode();}
	void SetInvertMitoCode() {code.SetInvertMitoCode();}
};


class AminoacidData : public SequenceData{

public:
	AminoacidData() : SequenceData(){
		maxNumStates = 20;
		fullyAmbigChar = maxNumStates;
		}

	AminoacidData(const NucleotideData *dat, int genCode) : SequenceData(){
		maxNumStates = 20;
		GeneticCode c;
		if(genCode == 0) c.SetStandardCode();
		else if(genCode == 1) c.SetVertMitoCode();
		else c.SetInvertMitoCode();
		FillAminoacidMatrixFromDNA(dat, &c);
		CopyNamesFromOtherMatrix(dat);
		fullyAmbigChar = maxNumStates;
		}
	void FillAminoacidMatrixFromDNA(const NucleotideData *dat, GeneticCode *code);
	void CalcEmpiricalFreqs();
	unsigned char CharToDatum(char d) ;
	void CreateMatrixFromNCL(const NxsCharactersBlock *);
	void CreateMatrixFromNCL(const NxsCharactersBlock *, NxsUnsignedSet &charset);
	};

#endif


