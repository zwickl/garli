// GARLI version 2.0 source code
// Copyright 2005-2011 Derrick J. Zwickl
// email: garli.support@gmail.com
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
		{ maxNumStates=4; strcpy( info, "DNA" ); empStateFreqs=NULL; numConditioningPatterns = 0;}
	SequenceData( int ntax, int nchar ) : DataMatrix( ntax, nchar )
		{ maxNumStates=4; strcpy( info, "DNA" ); empStateFreqs=NULL; numConditioningPatterns = 0;}
	virtual ~SequenceData() {
		if(empStateFreqs != NULL) delete []empStateFreqs;
		}

protected:
	FLOAT_TYPE *empStateFreqs;
	// overrides of base class's virtual fuctions
	virtual unsigned char CharToDatum( char ch ) const = 0;
	virtual unsigned char CharToBitwiseRepresentation( char ch ) const;
	virtual char	DatumToChar( unsigned char d ) const;
	virtual unsigned char	FirstState() const { return 0; }
	virtual unsigned char	LastState() const { return 3; }
	virtual int	NumStates(int) const { return 4; }

public:
	virtual void CreateMatrixFromNCL(const NxsCharactersBlock *, NxsUnsignedSet &charset) = 0;
	virtual void CalcEmpiricalFreqs() = 0;
	virtual void GetEmpiricalFreqs(FLOAT_TYPE *f) const{
		assert(empStateFreqs);
		for(int i=0;i<maxNumStates;i++) f[i]=empStateFreqs[i];
		}
	virtual void AddDummyRootToExistingMatrix();

	virtual bool IsNucleotide() const {return false;}
	virtual bool IsCodon() const {return false;}
	virtual bool IsAminoAcid() const {return false;}
	virtual bool IsNState() const {return false;}
	virtual bool IsOrientedGap() const {return false;}
	};

inline unsigned char SequenceData::CharToDatum( char ch ) const{
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
	}
	else
		throw ErrorException("Unknown character \"%c\" in SequenceData::CharToDatum", ch); 

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
	bool IsNucleotide() const {return true;}

	unsigned char CharToDatum(char d) const;
	void CalcEmpiricalFreqs();
	void CreateMatrixFromNCL(const NxsCharactersBlock *charblock, NxsUnsignedSet &charset);
	void MakeAmbigStrings();
	void AddDummyRootToExistingMatrix();
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
	//except for two-serine models, when they are 21 and the serine with two codons is state 20
	int codonTable[64];
	int map64toNonStops[64];
	vector<int> stops;

	//this holds the correspondence between the state indeces and actual codons
	//for display purposes.  Stops are removed and thus any mapIndexToCodonDisplay[index]
	//gives the codon for that index
	vector<string> mapIndexToCodonDisplay;

	public:
	enum{
		STANDARD= 0,
		VERTMITO = 1,
		INVERTMITO = 2,
		STANDARDTWOSERINE = 3,
		VERTMITOTWOSERINE = 4,
		INVERTMITOTWOSERINE = 5
		}codeName;
	
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
		
	void SetStandardTwoSerineCode(){
		//because the stops don't change location, I don't think that anything else needs to be changed here

		//the two lone serines become the 20th state
		codonTable[ 9 ]= 20;  //AGC
		codonTable[ 11 ]= 20; //AGT

		//the three stop codons become the 21st state
		codonTable[ 48 ]= 21;
		codonTable[ 50 ]= 21;
		codonTable[ 56 ]= 21;
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

	//this should be called AFTER SetVertMitoCode()
	void SetVertMitoTwoSerineCode(){
		//because the stops don't change location, I don't think that anything else needs to be changed here

		//the two lone serines become the 20th state
		codonTable[ 9 ]= 20;  //AGC
		codonTable[ 11 ]= 20; //AGT

		//the four stop codons become the 21st state
		codonTable[8] = 21;  //AGA
		codonTable[10] = 21;  //AGG
		codonTable[ 48 ]= 21;
		codonTable[ 50 ]= 21;
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

	//this should be called AFTER SetInvertMitoCode()
	void SetInvertMitoTwoSerineCode(){
		//because the stops don't change location, I don't think that anything else needs to be changed here

		//the two lone serines become the 20th state
		codonTable[ 9 ]= 20;  //AGC
		codonTable[ 11 ]= 20; //AGT

		//the two stop codons become the 21st state
		codonTable[ 48 ]= 21;
		codonTable[ 50 ]= 21;
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
						sprintf(cod, "%c%c%c", nucs[f], nucs[s], nucs[t]);
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

	int NumStates() const {return mapIndexToCodonDisplay.size();}
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

	CodonData(const NucleotideData *dat, int genCode, bool ignoreStops=false) : SequenceData(){
		assert(dat->Dense() == false);
		if(genCode == GeneticCode::STANDARD){
			code.SetStandardCode();
			maxNumStates = 61;
			}
		else if(genCode == GeneticCode::VERTMITO){
			code.SetVertMitoCode();
			maxNumStates = 60;
			}
		else if(genCode == GeneticCode::INVERTMITO){
			code.SetInvertMitoCode();
			maxNumStates = 62;
			}
		else{
			throw ErrorException("Sorry, only the standard, vert mito and invert mito codes can be used with codon models");
			}
		usePatternManager = dat->GetUsePatternManager();
		FillCodonMatrixFromDNA(dat, ignoreStops);
		CopyNamesFromOtherMatrix(dat);
		empType = NOT_EMPIRICAL;
		fullyAmbigChar = maxNumStates;
		}

	~CodonData(){}
	bool IsCodon() const {return true;}

	void FillCodonMatrixFromDNA(const NucleotideData *, bool ignoreStops);
	unsigned char CharToDatum(char c) const{
		//this shouldn't be getting called, as it makes no sense for codon data
		assert(0);
		return 0;
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

	AminoacidData(const NucleotideData *dat, int genCode, bool ignoreStops=false) : SequenceData(){
		maxNumStates = 20;
		GeneticCode c;
		if(genCode == GeneticCode::STANDARD) c.SetStandardCode();
		else if(genCode == GeneticCode::VERTMITO) c.SetVertMitoCode();
		else if(genCode == GeneticCode::INVERTMITO) c.SetInvertMitoCode();
		else{
			if(genCode == GeneticCode::STANDARDTWOSERINE){
				c.SetStandardTwoSerineCode();
				}
			else if(genCode == GeneticCode::VERTMITOTWOSERINE){
				c.SetVertMitoCode();
				c.SetVertMitoTwoSerineCode();
				}
			else if(genCode == GeneticCode::INVERTMITOTWOSERINE){
				c.SetInvertMitoCode();
				c.SetInvertMitoTwoSerineCode();
				}
			else assert(0);
			maxNumStates = 21;	
			}
		usePatternManager = dat->GetUsePatternManager();
		FillAminoacidMatrixFromDNA(dat, &c, ignoreStops);
		CopyNamesFromOtherMatrix(dat);
		fullyAmbigChar = maxNumStates;
		}
	bool IsAminoAcid() const {return true;}
	void FillAminoacidMatrixFromDNA(const NucleotideData *dat, GeneticCode *code, bool ignoreStops);
	void CalcEmpiricalFreqs();
	unsigned char CharToDatum(char d) const;
	void CreateMatrixFromNCL(const NxsCharactersBlock *, NxsUnsignedSet &charset);
	};

class DataPartition {
private:
	vector<SequenceData *> dataSubsets;
	int nTax;
public:
	void AddSubset(SequenceData* sub){
		dataSubsets.push_back(sub);
		nTax = sub->NTax();
		}
	SequenceData *GetSubset(int num) const{
		if(num < 0 || (num < dataSubsets.size()) == false) 
			throw ErrorException("Tried to access invalid subset number");
		return dataSubsets[num];
		}
	void Delete(){
		for(vector<SequenceData *>::iterator it = dataSubsets.begin();it != dataSubsets.end(); it++)
			delete *it;
		dataSubsets.clear();
		}
	int NTax() const {return nTax;}
	int NumSubsets() const {return dataSubsets.size();}
	void BeginNexusTreesBlock(string &trans) const {dataSubsets[0]->BeginNexusTreesBlock(trans);}
	void BeginNexusTreesBlock(ofstream &out) const {dataSubsets[0]->BeginNexusTreesBlock(out);}
	NxsString TaxonLabel(int t) const {return dataSubsets[0]->TaxonLabel(t);}
	int TaxonNameToNumber(NxsString name) const{return dataSubsets[0]->TaxonNameToNumber(name);}

	void AddDummyRoots(){
		nTax++;
		for(int p = 0;p < NumSubsets();p++){
			dataSubsets[p]->AddDummyRootToExistingMatrix();
			assert(nTax == dataSubsets[p]->NTax());
			}
		}
	int BootstrapReweight(int seedToUse, FLOAT_TYPE resampleProportion){
		int nextSeed = seedToUse;
		for(int p = 0;p < NumSubsets();p++){
			outman.UserMessage("\tSubset %d: Random seed for bootstrap reweighting: %d", p + 1, nextSeed);
			SequenceData *curData = GetSubset(p);
			nextSeed = curData->BootstrapReweight(nextSeed, resampleProportion);
			}
		return nextSeed;
		}
	};

class DataSubsetInfo{
public:
	int garliSubsetNum;
	int charblockNum;
	string charblockName;
	int partitionSubsetNum;
	string partitionSubsetName;
	enum type{
		NUCLEOTIDE = 0,
		AMINOACID = 1,
		CODON = 2, 
		NSTATE= 3,
		NSTATEV = 4,
		ORDNSTATE = 5,
		ORDNSTATEV = 6,
		ORIENTEDGAP = 7,
		BINARY = 8,
		BINARY_NOT_ALL_ZEROS = 9
		}readAs, usedAs;
	int totalCharacters;
	int uniqueCharacters;
	string outputNames[10];//{"Nucleotide data", "Amino acid data", "Codon data"};
	DataSubsetInfo(int gssNum, int cbNum, string cbName, int psNum, string psName, type rAs, type uAs) :
		garliSubsetNum(gssNum), charblockNum(cbNum), charblockName(cbName), partitionSubsetNum(psNum), partitionSubsetName(psName), readAs(rAs), usedAs(uAs){
			outputNames[NUCLEOTIDE]="Nucleotide data";
			outputNames[AMINOACID]="Amino acid data";
			outputNames[CODON]="Codon data";
			outputNames[NSTATE]="Standard k-state data";
			outputNames[NSTATEV]="Standard k-state data, variable only";
			outputNames[ORDNSTATE]="Standard ordered k-state data";
			outputNames[ORDNSTATEV]="Standard ordered k-state data, variable only";
			outputNames[ORIENTEDGAP]="Gap-coded data, oriented with respect to time";
			outputNames[BINARY]="Binary data";
			outputNames[BINARY_NOT_ALL_ZEROS]="Binary data, no constant state 0 chars";
			}
	void Report(){
		outman.UserMessage("GARLI data subset %d", garliSubsetNum+1);
		outman.UserMessage("\tCHARACTERS block #%d (\"%s\")", charblockNum+1, charblockName.c_str());
		if(partitionSubsetNum >= 0) outman.UserMessage("\tCHARPARTITION subset #%d (\"%s\")", partitionSubsetNum+1, partitionSubsetName.c_str());
		outman.UserMessage("\tData read as %s,\n\tmodeled as %s", outputNames[readAs].c_str(), outputNames[usedAs].c_str());
		}
	};
/*
//
// Mk type model, with binary data
class BinaryData : public SequenceData{
	public:
		BinaryData() : SequenceData(){
			maxNumStates = 2;
			}

		unsigned char CharToDatum(char d);
		char DatumToChar( unsigned char d );
		void CreateMatrixFromNCL(const NxsCharactersBlock *, NxsUnsignedSet &charset);
		void CalcEmpiricalFreqs(){
			//BINARY - this might actually make sense for gap encoding
			}
		//this is just a virtual overload that avoids doing anything if determine const is called with inappropriate data
		void DetermineConstantSites(){};
	};

inline unsigned char BinaryData::CharToDatum( char ch ){
	unsigned char datum;

	if( ch == '0' || ch == '-' )
		datum = 0;
	else if( ch == '1' || ch == '+' )
		datum = 1;
	else if( ch == '?' )
		datum = 2;
	else
      THROW_BADSTATE(ch);

	return datum;
	}

inline char BinaryData::DatumToChar( unsigned char d ){
	char ch = 'X';	// ambiguous

	if( d == 2 )
		ch = '?';
	else if( d == 0 )
		ch = '0';
	else if( d == 1 )
		ch = '1';

	return ch;
	}
*/
//
// Mk or Mkv type model, with n-state data
class NStateData : public SequenceData{
	public:
		enum{
			ALL = 0,
			ONLY_VARIABLE = 1,
			ONLY_INFORM = 2,
			BINARY = 3,
			BINARY_NOT_ALL_ZEROS = 4
			}datatype;
		enum{
			UNORDERED = 0,
			ORDERED = 1
			}modeltype;
		NStateData() : SequenceData(){
			maxNumStates = 99;
			}
		NStateData(int ns) : SequenceData(){
			maxNumStates = ns;
			}
		//NStateData(int ns, bool isMkv, bool isOrdered) : SequenceData(){'
		NStateData(int ns, bool isOrdered, bool isBinary, bool isConditioned) : SequenceData(){
			if(isBinary){
				if(isConditioned)
					datatype = BINARY_NOT_ALL_ZEROS;
				else
					datatype = BINARY;
				}
			else if(isConditioned)
				datatype = ONLY_VARIABLE;
			else 
				datatype = ALL;
			if(isOrdered)
				modeltype = ORDERED;
			else
				modeltype = UNORDERED;
			maxNumStates = ns;
			}
		bool IsNState() const {return true;}
		void SetNumStates(int ns){maxNumStates = ns;}

		virtual unsigned char CharToDatum(char d) const;
		virtual char DatumToChar( unsigned char d ) const;
		virtual void CreateMatrixFromNCL(const NxsCharactersBlock *, NxsUnsignedSet &charset);
		void CalcEmpiricalFreqs(){
			//BINARY - this might actually make sense for gap encoding
			}
		//this is a virtual overload for NState because it might have to deal with the conditioning chars, which shouldn't be included in the resampling
		int BootstrapReweight(int restartSeed, FLOAT_TYPE resampleProportion);
		//this is just a virtual overload that avoids doing anything if determine const is called with inappropriate data
		void DetermineConstantSites(){};
	};

inline unsigned char NStateData::CharToDatum( char ch ) const{
	unsigned char datum;

	if( ch == '0')
		datum = 0;
	else if( ch == '1')
		datum = 1;
	else if( ch == '2')
		 datum = 2;
	else if( ch == '3')
		 datum = 3;
	else if( ch == '4')
		 datum = 4;
	else if( ch == '5')
		 datum = 5;
	else if( ch == '6')
		 datum = 6;
	else if( ch == '7')
		 datum = 7;
	else if( ch == '8')
		 datum = 8;
	else if( ch == '9')
		 datum = 9;
	else if( ch == '?')
		datum = 99;
	else
      throw ErrorException("Unknown character \"%c\" in NStateData::CharToDatum", ch); 

	return datum;
	}

inline char NStateData::DatumToChar( unsigned char d ) const{
	//NSTATE - not sure how this should work, but it isn't that important anyway
	
	char ch = 'X';	// ambiguous
/*
	if( d == 2 )
		ch = '?';
	else if( d == 0 )
		ch = '0';
	else if( d == 1 )
		ch = '1';
*/
	return ch;
	}

class OrientedGapData : public NStateData{
	public:
		OrientedGapData() : NStateData(){
			maxNumStates = 2;
			}
		OrientedGapData(int ns) : NStateData(){
			assert(0);
			}
		OrientedGapData(int ns, bool isMkv) : NStateData(){
			assert(0);
			}
		bool IsOrientedGap() const {return true;}
		void SetNumStates(int ns){maxNumStates = ns;}

		virtual void CreateMatrixFromNCL(const NxsCharactersBlock *, NxsUnsignedSet &charset);
		void CalcEmpiricalFreqs(){
			//BINARY - this might actually make sense for gap encoding
			}
		//this is just a virtual overload that avoids doing anything if determine const is called with inappropriate data
		void DetermineConstantSites(){};
	};


#endif


