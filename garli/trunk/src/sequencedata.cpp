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

#include "defs.h"
#include "sequencedata.h"
#include "garlireader.h"

#undef DEBUG_CALCFREQ
#undef DEBUG_CALCPRMATRIX
#undef DEBUGGING_PRMATRICES
#undef DEBUGGING_PATTERN_PROBS

#if defined( DEBUGGING_PATTERN_PROBS )
#	include <io.h>
#endif

void NucleotideData::CalcEmpiricalFreqs(){

	empStateFreqs=new FLOAT_TYPE[4];//this is a member of the class, and where the final freqs will be stored
	empStateFreqs[0]=empStateFreqs[1]=empStateFreqs[2]=empStateFreqs[3]=0.0;
	
	//these are all temporary and local
	FLOAT_TYPE freqSumNoAmbig[4] = {0.0, 0.0, 0.0, 0.0};
	FLOAT_TYPE freqSumAmbig[4]  = {0.0, 0.0, 0.0, 0.0};
	FLOAT_TYPE nonAmbigTotal = 0.0;
	FLOAT_TYPE ambigTotal = 0.0;

	vector<char> ambigStates;
	vector<int> ambigCounts;
	for( int i = 0; i < NTax(); i++ ) {
		for( int j = 0; j < NChar(); j++ ) {
			char thischar=(char) Matrix( i, j );
			int nstates=0;
			//first figure out how many states we've got
			if(thischar & 1) nstates++;
			if(thischar & 2) nstates++;
			if(thischar & 4) nstates++;
			if(thischar & 8) nstates++;			

			if(nstates==1){
				if(thischar & 1)
					freqSumNoAmbig[0] += (FLOAT_TYPE) Count(j);
				if(thischar & 2)
					freqSumNoAmbig[1] += (FLOAT_TYPE) Count(j);
				if(thischar & 4)
					freqSumNoAmbig[2] += (FLOAT_TYPE) Count(j);
				if(thischar & 8) 
					freqSumNoAmbig[3] += (FLOAT_TYPE) Count(j);	
				nonAmbigTotal += Count(j);
				}
			else if(nstates < 4){
				//now divide the states up to the bases
				//division will be equal for this pass, and refined below
				if(thischar & 1)
					freqSumAmbig[0] += (FLOAT_TYPE) Count(j)/nstates;
				if(thischar & 2)
					freqSumAmbig[1] += (FLOAT_TYPE) Count(j)/nstates;
				if(thischar & 4)
					freqSumAmbig[2] += (FLOAT_TYPE) Count(j)/nstates;
				if(thischar & 8) 
					freqSumAmbig[3] += (FLOAT_TYPE) Count(j)/nstates;
				ambigTotal += Count(j);

				//these will store a list of the ambiguous characters so that iterations
				//below don't require going through the whole dataset again
				ambigStates.push_back(thischar);
				ambigCounts.push_back(Count(j));					
				}
			}
		}
	
	for(int j=0;j<4;j++)
		empStateFreqs[j] = (freqSumNoAmbig[j] + freqSumAmbig[j]) / (nonAmbigTotal + ambigTotal);

	//now iterate to refine the emp freqs to account for partial ambiguity
	if(ambigStates.size() > 0){
		bool continueIterations;
		do{
			continueIterations = false;
			freqSumAmbig[0]=freqSumAmbig[1]=freqSumAmbig[2]=freqSumAmbig[3]=0.0;
			for(unsigned i=0;i<ambigStates.size();i++){
				FLOAT_TYPE fracSum = 0.0;
				int nstates = 0;
				char thischar = ambigStates[i];
				
				if(thischar & 1)
					fracSum += empStateFreqs[0];
				if(thischar & 2)
					fracSum += empStateFreqs[1];
				if(thischar & 4)
					fracSum += empStateFreqs[2];
				if(thischar & 8)
					fracSum += empStateFreqs[3];
				
				//this time they are allocated to the bases in proportion to the total
				//frequencies from the last iteration
				if(thischar & 1)
					freqSumAmbig[0] += (FLOAT_TYPE) ambigCounts[i] * (empStateFreqs[0]/fracSum);
				if(thischar & 2)
					freqSumAmbig[1] += (FLOAT_TYPE) ambigCounts[i] * (empStateFreqs[1]/fracSum);
				if(thischar & 4)
					freqSumAmbig[2] += (FLOAT_TYPE) ambigCounts[i] * (empStateFreqs[2]/fracSum);
				if(thischar & 8) 
					freqSumAmbig[3] += (FLOAT_TYPE) ambigCounts[i] * (empStateFreqs[3]/fracSum);
				}
			FLOAT_TYPE tempFreqs[4] = {0.0, 0.0, 0.0, 0.0};
			for(int j=0;j<4;j++){
				tempFreqs[j] = (freqSumNoAmbig[j] + freqSumAmbig[j]) / (nonAmbigTotal + ambigTotal);
				if(fabs(tempFreqs[j] - empStateFreqs[j]) >  max(1.0e-8, GARLI_FP_EPS * 2.0)) continueIterations = true;
				empStateFreqs[j] = tempFreqs[j];
				}
			}while(continueIterations);
		}	
	

#if defined( DEBUG_CALCFREQ )
	cerr << endl << "Frequency of A: " << p[0] << endl;
	cerr << "Frequency of C: " << p[1] << endl;
	cerr << "Frequency of G: " << p[2] << endl;
	cerr << "Frequency of T: " << p[3] << endl;
	cerr << "Total         : " << ( p[0] + p[1] + p[2] + p[3] ) << endl;

	cerr << endl << "Program stopped after calculating base frequencies because" << endl;
	cerr << "DEBUG_CALCFREQS was #define'd in source code file \"mlhky.cpp\" " << endl;

	cerr << endl << "Press Enter key to continue..." << endl;
	char ch = '\0';
	cin.get(ch);
	exit(0);
#endif
}

void AminoacidData::CalcEmpiricalFreqs(){
	empStateFreqs=new FLOAT_TYPE[maxNumStates];//this is a member of the class, and where the final freqs will be stored
	
	for(int i=0;i<maxNumStates;i++) empStateFreqs[i] = 0.0;
	FLOAT_TYPE total = 0.0;
	//for codons and aminoacids this will assume no ambiguity
	
	for( int i = 0; i < NTax(); i++ ) {
		for( int j = 0; j < nChar; j++ ) {
			char thischar= matrix[i][j];

			if(thischar != maxNumStates){
				assert(thischar > -1);
				empStateFreqs[thischar] += (FLOAT_TYPE) Count(j);
				total += (FLOAT_TYPE) Count(j);
				}
			}
		}
	//check whether this might be nucleotide data in disguise
	if((empStateFreqs[0]+empStateFreqs[1]+empStateFreqs[5]+empStateFreqs[16])/total > 0.90) throw ErrorException("Model specified as aminoacid, but nucleotide data found!");
		
	FLOAT_TYPE freqTot = 0.0;
	bool allPresent = true;
	for(int j=0;j<maxNumStates;j++) if(empStateFreqs[j] == ZERO_POINT_ZERO) allPresent = false;
	if(!allPresent){
		outman.UserMessage("WARNING: Not all amino acids were observed in this dataset.\n\tOne pseudo-count will be added to each amino acid for calculation of the\n\tempirical frequencies. You should probably use\n\ta statefrequencies setting other than emprical.\n");
		for(int j=0;j<maxNumStates;j++) empStateFreqs[j] += ONE_POINT_ZERO;
		total += (FLOAT_TYPE) maxNumStates;
		}
	for(int j=0;j<maxNumStates;j++){
		empStateFreqs[j] /= total;
		freqTot += empStateFreqs[j];
		}
	assert(fabs(freqTot - 1.0) < 1e-5);
	}

void CodonData::CalcEmpiricalFreqs(){
	if(empType == NOT_EMPIRICAL) return;

	empStateFreqs=new FLOAT_TYPE[maxNumStates];//this is a member of the class, and where the final freqs will be stored

	if(empType == F1X4){
		CalcF1x4Freqs();
		BaseFreqXPositionReport();
		return;
		}
	if(empType == F3X4){
		CalcF3x4Freqs();
		BaseFreqXPositionReport();
		return;
		}

	//we must be using the actual observed frequencies
	assert(empType == CODON_TABLE);
	for(int i=0;i<maxNumStates;i++) empStateFreqs[i] = 0.0;
	FLOAT_TYPE total = 0.0;
	//for codons and aminoacids this will assume no ambiguity
	
	for( int i = 0; i < NTax(); i++ ) {
		for( int j = 0; j < nChar; j++ ) {
			char thischar= matrix[i][j];

			if(thischar != maxNumStates){
				assert(thischar > -1);
				empStateFreqs[thischar] += (FLOAT_TYPE) Count(j);
				total += (FLOAT_TYPE) Count(j);
				}
			}
		}
	FLOAT_TYPE freqTot = 0.0;

	bool allPresent = true;
	for(int j=0;j<maxNumStates;j++) if(empStateFreqs[j] == ZERO_POINT_ZERO) allPresent = false;
	if(!allPresent){
		outman.UserMessage("WARNING: Not all allowable codons were observed in this dataset.\n\tOne pseudo-count will be added to each codon for caluclation of the\n\tempirical frequencies. You should probably use\n\tstatefrequencies = f1x4 or f3x4 instead of empirical.\n");
		for(int j=0;j<maxNumStates;j++) empStateFreqs[j] += ONE_POINT_ZERO;
		total += (FLOAT_TYPE) maxNumStates;
		}

	for(int j=0;j<maxNumStates;j++){
		empStateFreqs[j] /= total;
		freqTot += empStateFreqs[j];
		}
	assert(fabs(freqTot - 1.0) < 1e-5);
	}

void NucleotideData::MakeAmbigStrings(){
	//this will populate the ambigStrings vector with the data in the typical ambiguity format
	
	ambigStrings.reserve(NTax());

	for(int i=0;i<NTax();i++){
		unsigned char* thisdata=GetRow(i);
		
		//run through all the characters once just to see how many states we have so we can allocate
		//the correct length of array to hold the string
		int totalStates=0;
		for(int j=0;j<NChar();j++){
			char thisbase=thisdata[j];
			int numstates=0;
			if(thisbase&1){
				numstates++;
				}
			if(thisbase&2){
				numstates++;
				}
			if(thisbase&4){
				numstates++;
				}
			if(thisbase&8){
				numstates++;
				}
			if(numstates!=4) totalStates += numstates;
			//remember that if we have ambiguity we need an extra character to hold the number of states
			//and if we have total ambiguity (numstates=0) we also need a character to hold that
			if(numstates>1 || numstates==0 || numstates==4) totalStates++;
			}
	
		char *thisString=new char[totalStates];

#ifdef OPEN_MP
		unsigned *thisMap=new unsigned[NChar()];
#endif

		//now do it for real
		int index=0;
		for(int j=0;j<NChar();j++){

#ifdef OPEN_MP
			thisMap[j]=index;
#endif
			char thisbase=thisdata[j];
			int numstates=0;
			char thiscode;
			if(thisbase&1){
				numstates++;
				thiscode=0;
				}
			if(thisbase&2){
				numstates++;
				thiscode=1;
				}
			if(thisbase&4){
				numstates++;
				thiscode=2;
				}
			if(thisbase&8){
				numstates++;
				thiscode=3;
				}
			
			if(numstates==1){
				thisString[index++]=thiscode;
				}
			else if(numstates==4||numstates==0){
				thisString[index++] = -4;
				}
			else{
				thisString[index++] = -numstates;
				if(thisbase&1) thisString[index++] = 0;
				if(thisbase&2) thisString[index++] = 1;
				if(thisbase&4) thisString[index++] = 2;			
				if(thisbase&8) thisString[index++] = 3;		
				}
			}
		ambigStrings.push_back(thisString);
#ifdef OPEN_MP
		ambigToCharMap.push_back(thisMap);
#endif
		}
	}

void CodonData::FillCodonMatrixFromDNA(const NucleotideData *dnaData){
	//first we need to convert the nucleotide data to codons numbered 0-60 or 61 and assign them back to the terminals
	//codons are ordered AAA, AAC, AAG, AAT, ACA, ... TTT
	short pos1, pos2, pos3;

	nonZeroCharCount = nChar = dnaData->NChar()/3;
	nTax = dnaData->NTax();
	if(dnaData->NChar() % 3 != 0) throw ErrorException("Codon datatype specified, but number of nucleotides not divisible by 3!");  
	NewMatrix(nTax, nChar);

	//this will just map from the bitwise format to the index format (A, C, G, T = 0, 1, 2, 3)
	//partial ambiguity is mapped to total ambiguity currently
	short bitwiseToIndexFormat[16] = {15,0,1,15,2,15,15,15,3,15,15,15,15,15,15,15};

	//keep track of the empirical base freqs at the codon positions, for possible use 
	//in the F1x4 or F3x4 methods of calculating the equilibrium codon freqs
	empBaseFreqsPos1[0]=empBaseFreqsPos1[1]=empBaseFreqsPos1[2]=empBaseFreqsPos1[3]=ZERO_POINT_ZERO;
	empBaseFreqsPos2[0]=empBaseFreqsPos2[1]=empBaseFreqsPos2[2]=empBaseFreqsPos2[3]=ZERO_POINT_ZERO;
	empBaseFreqsPos3[0]=empBaseFreqsPos3[1]=empBaseFreqsPos3[2]=empBaseFreqsPos3[3]=ZERO_POINT_ZERO;

	FLOAT_TYPE total = ZERO_POINT_ZERO;

	int tax=0, thisCodonNum;
	for(int tax=0;tax<NTax();tax++){
		bool firstAmbig = true;
		for(int cod=0;cod<nChar;cod++){
			short p1 = dnaData->Matrix(tax, cod*3);
			short p2 = dnaData->Matrix(tax, cod*3+1);
			short p3 = dnaData->Matrix(tax, cod*3+2);

			pos1 = bitwiseToIndexFormat[p1];
			pos2 = bitwiseToIndexFormat[p2];
			pos3 = bitwiseToIndexFormat[p3];
			
			thisCodonNum=(pos1)*16 + (pos2)*4 + pos3;
			
			if(pos1==15||pos2==15||pos3==15){//check for gaps or ambiguity
				if(pos1+pos2+pos3 != 45){
					//warn about gaps or ambiguity in codons
					if(firstAmbig){
						outman.UserMessageNoCR("Gaps or ambiguity codes found within codon for taxon %s.\n\tCodons coded as missing for that taxon: ", dnaData->TaxonLabel(tax));
						firstAmbig = false;
						}
					outman.UserMessageNoCR("%d ", cod+1);
					}
				thisCodonNum=64;
				}
			else{
				empBaseFreqsPos1[pos1] += ONE_POINT_ZERO;
				empBaseFreqsPos2[pos2] += ONE_POINT_ZERO;
				empBaseFreqsPos3[pos3] += ONE_POINT_ZERO;
				total += ONE_POINT_ZERO;
				}

			char prot;
			//note that a return code of 20 from the codon lookup indicates a stop codon, but a protein code of 20 generally means total ambiguity
			if(thisCodonNum != 64){
				prot = code.CodonLookup(thisCodonNum);
				if(prot == 20){
					string c;
					char b[4]={'A','C','G','T'};
					c += b[pos1];
					c += b[pos2];
					c += b[pos3];
					throw ErrorException("stop codon %s found at codon site %d (nuc site %d) in taxon %s.  Bailing out.", c.c_str(), cod+1, cod*3+1,  dnaData->TaxonLabel(tax));
					}
				}

			if(thisCodonNum == 64)//missing or ambiguous 
				matrix[tax][cod] = maxNumStates;
			else 
				matrix[tax][cod] = code.Map64stateToNonStops(thisCodonNum);
			}
		if(firstAmbig == false) outman.UserMessage("");
		}
	for(int b=0;b<4;b++){
		empBaseFreqsAllPos[b] = (empBaseFreqsPos1[b] + empBaseFreqsPos2[b] + empBaseFreqsPos3[b]) / (3.0 * total);

		empBaseFreqsPos1[b] /= total;
		empBaseFreqsPos2[b] /= total;
		empBaseFreqsPos3[b] /= total;
		}
	}

void AminoacidData::FillAminoacidMatrixFromDNA(const NucleotideData *dnaData, GeneticCode *code){
	//first we need to convert the nucleotide data to codons, and then translate the codons to AA's
	//codons are ordered AAA, AAC, AAG, AAT, ACA, ... TTT
	short pos1, pos2, pos3;

	nonZeroCharCount = nChar = dnaData->NChar()/3;
	nTax = dnaData->NTax();
	if(dnaData->NChar() % 3 != 0) throw ErrorException("Codon to Aminoacid translation specified, but number of nucleotides not divisible by 3!");  
	NewMatrix(nTax, nChar);

	//this will just map from the bitwise format to the index format (A, C, G, T = 0, 1, 2, 3)
	//partial ambiguity is mapped to total ambiguity currently
	short bitwiseToIndexFormat[16] = {15,0,1,15,2,15,15,15,3,15,15,15,15,15,15,15};
	
	int tax=0, thisCodonNum;
	for(int tax=0;tax<NTax();tax++){
		bool firstAmbig = true;
		for(int cod=0;cod<nChar;cod++){
			short p1 = dnaData->Matrix(tax, cod*3);
			short p2 = dnaData->Matrix(tax, cod*3+1);
			short p3 = dnaData->Matrix(tax, cod*3+2);

			pos1 = bitwiseToIndexFormat[p1];
			pos2 = bitwiseToIndexFormat[p2];
			pos3 = bitwiseToIndexFormat[p3];
			
			thisCodonNum=(pos1)*16 + (pos2)*4 + pos3;
			
			if(pos1==15||pos2==15||pos3==15){//check for gaps
				if(pos1+pos2+pos3 != 45){
					//warn about gaps or ambiguity in codons
					if(firstAmbig){
						outman.UserMessageNoCR("Gaps or ambiguity codes found within codon for taxon %s.\n\tAminoacids coded as missing for that taxon: ", dnaData->TaxonLabel(tax));
						firstAmbig = false;
						}
					outman.UserMessageNoCR("%d ", cod+1);
					}
				thisCodonNum=64;
				}

			char prot;
			//note that a return code of 20 from the codon lookup indicates a stop codon, but a protein code of 20 generally means total ambiguity
			if(thisCodonNum != 64){
				prot = code->CodonLookup(thisCodonNum);
				if(prot == 20){
					string c;
					char b[4]={'A','C','G','T'};
					c += b[pos1];
					c += b[pos2];
					c += b[pos3];
					throw ErrorException("stop codon %s found at codon site %d (nuc site %d) in taxon %s.  Bailing out.", c.c_str(), cod+1, cod*3+1,  dnaData->TaxonLabel(tax));
					}
				}
			else prot = 20;

			matrix[tax][cod] = prot;
			}
		if(firstAmbig == false) outman.UserMessage("");
		}	
	}

void CodonData::CalcF1x4Freqs(){
	//this assumes that the empirical base freqs have already been calculated in FillCodonMatrixFromDNA
	assert(fabs(empBaseFreqsAllPos[0] + empBaseFreqsAllPos[1] + empBaseFreqsAllPos[2] + empBaseFreqsAllPos[3] - 1.0) < 1.0e-4);

	FLOAT_TYPE total = ZERO_POINT_ZERO;

	int stops=0;
	for(int base1=0;base1<4;base1++){
		for(int base2=0;base2<4;base2++){
			for(int base3=0;base3<4;base3++){
				if(code.CodonLookup(base1*16+base2*4+base3) != 20){
					empStateFreqs[base1*16+base2*4+base3 - stops] = empBaseFreqsAllPos[base1] * empBaseFreqsAllPos[base2] * empBaseFreqsAllPos[base3];
					total += empStateFreqs[base1*16+base2*4+base3 - stops];
					}
				else stops++;
				}
			}
		}
	//now normalize, because the stop codons will make the total of the 60 or 61 allowed codons < 1.0
	for(int s=0;s<maxNumStates;s++) empStateFreqs[s] /= total;
	}

void CodonData::CalcF3x4Freqs(){
	//this assumes that the empirical base freqs have already been calculated in FillCodonMatrixFromDNA
	assert(fabs(empBaseFreqsPos1[0] + empBaseFreqsPos1[1] + empBaseFreqsPos1[2] + empBaseFreqsPos1[3] - 1.0) < 1.0e-4);

	if((empBaseFreqsPos1[0] == ZERO_POINT_ZERO) ||
	   (empBaseFreqsPos1[1] == ZERO_POINT_ZERO) ||
	   (empBaseFreqsPos1[2] == ZERO_POINT_ZERO) ||
	   (empBaseFreqsPos1[3] == ZERO_POINT_ZERO) ||
	   (empBaseFreqsPos2[0] == ZERO_POINT_ZERO) ||
	   (empBaseFreqsPos2[1] == ZERO_POINT_ZERO) ||
	   (empBaseFreqsPos2[2] == ZERO_POINT_ZERO) ||
	   (empBaseFreqsPos2[3] == ZERO_POINT_ZERO) ||
	   (empBaseFreqsPos3[0] == ZERO_POINT_ZERO) ||
	   (empBaseFreqsPos3[1] == ZERO_POINT_ZERO) ||
	   (empBaseFreqsPos3[2] == ZERO_POINT_ZERO) ||
	   (empBaseFreqsPos3[3] == ZERO_POINT_ZERO)) throw ErrorException("All bases were not observed at all codon positions!\n\tYou probably shouldn't be using F3x4 to estimate equilibrium codon frequencies.\n\tUse F1x4 or equal state frequencies.");

	FLOAT_TYPE total = ZERO_POINT_ZERO;
	
	int stops=0;
	for(int base1=0;base1<4;base1++){
		for(int base2=0;base2<4;base2++){
			for(int base3=0;base3<4;base3++){
				if(code.CodonLookup(base1*16+base2*4+base3) != 20){
					empStateFreqs[base1*16+base2*4+base3 - stops] = empBaseFreqsPos1[base1] * empBaseFreqsPos2[base2] * empBaseFreqsPos3[base3];
					total += empStateFreqs[base1*16+base2*4+base3 - stops];
					}
				else stops++;
				}
			}
		}
	//now normalize, because the stop codons will make the total of the 60 or 61 allowed codons < 1.0
	for(int s=0;s<maxNumStates;s++) empStateFreqs[s] /= total;
	}

void CodonData::BaseFreqXPositionReport(){
	//positional frequency report
	outman.UserMessage("Base usage at codon positions:");
	outman.UserMessage("       %10c%10c%10c%10c", 'A', 'C', 'G', 'T');
	outman.UserMessage(" pos 1   %10.5f%10.5f%10.5f%10.5f", empBaseFreqsPos1[0], empBaseFreqsPos1[1], empBaseFreqsPos1[2], empBaseFreqsPos1[3]);
	outman.UserMessage(" pos 2   %10.5f%10.5f%10.5f%10.5f", empBaseFreqsPos2[0], empBaseFreqsPos2[1], empBaseFreqsPos2[2], empBaseFreqsPos2[3]);
	outman.UserMessage(" pos 3   %10.5f%10.5f%10.5f%10.5f", empBaseFreqsPos3[0], empBaseFreqsPos3[1], empBaseFreqsPos3[2], empBaseFreqsPos3[3]);
	outman.UserMessage(" all pos %10.5f%10.5f%10.5f%10.5f\n", empBaseFreqsAllPos[0], empBaseFreqsAllPos[1], empBaseFreqsAllPos[2], empBaseFreqsAllPos[3]);
	}

//
// ComparePatterns returns:
//	 0		complete identity
//	-1		if i less than j
//	 1		if i greater than j
//
/*
int CodonData::ComparePatterns( const int i, const int j ) const{
	int cmp = 0;

	for( int k = 0; k < nTax; k++ ) {
		int same = ( Matrix( k, i ) == Matrix( k, j ) );
		if( !same )	{
			FLOAT_TYPE diff = ( (FLOAT_TYPE)Matrix( k, i ) - (FLOAT_TYPE)Matrix( k, j ) );
			cmp = ( diff < 0.0 ? -1 : 1 );
			break;
		}
	}
	return cmp;
}
*/

// PatternType determines whether pattern k is constant, informative, or autoapomorphic
//
/*
int CodonData::PatternType( int k , int *c, unsigned char *s) const
{
	if( k >= NChar() )
		return 0;
	int i, j, retval;

	for( i = 0; i < nTax; i++ )
		c[i] = 1;

	for( i = 0; i < nTax; i++ )
		s[i] = Matrix( i, k );

	// sort elements of s
	for( i = 0; i < nTax-1; i++ ) {
		for( j = i+1; j < nTax; j++ ) {
			if( s[i] > s[j] ) {
				unsigned char tmp = s[i];
				s[i] = s[j];
				s[j] = tmp;
			}
		}
	}
	
	// add counts of duplicate elements of s to first instance
	int nStates = 1; 	// treats ? as a new state
	bool ambig = false;	// will be true if any ? found
	bool allMissing = true;
	i = 0;
	for( j = 1; j < nTax; j++ ) {
		if(s[j]!=maxNumStates)  allMissing=false;
		if( s[j] == s[i] ) {
			c[i]++;
			c[j]--;
			}
		else {
			i = j;
			nStates++;
			}
		}

	//DJZ 10/28/03 changing this to allow for invariant sites.  Sites which contain 
	//some missing data but are otherwise constant must be marked as such because they 
	//will be considered constant for the purposes of invariant sites calcs.
	//also marking sites that are all missing

//	if( nStates == 1 )
	if( nStates == 1 /*|| (nStates==2 && missing)*/
/*
)
		retval = PT_CONSTANT;
	else if( nStates == 2 && ( c[0] == 1 || c[0] == nTax-1 ) )
		retval = PT_AUTAPOMORPHIC | PT_VARIABLE;
	else if( nStates < nTax )
		retval = PT_INFORMATIVE | PT_VARIABLE;
	else
		retval = PT_VARIABLE;

//	MEM_DELETE_ARRAY(s); // s is of length nTax
//	MEM_DELETE_ARRAY(c); // c is of length nTax

//	numStates[k] = ( missing ? nStates-1 : nStates );
	if(allMissing) nStates=0;
	numCodonStates[k] = nStates;

	return retval;
}
*/
/*
void CodonData::NewCodonMatrix( int taxa, int sites ){

	// delete data matrix and count and number arrays
	if( codonMatrix ) {
		int j;
		for( j = 0; j < nTax; j++ )
			MEM_DELETE_ARRAY(codonMatrix[j]);
		MEM_DELETE_ARRAY(codonMatrix); 
		}

	if( codonCount ) {
		MEM_DELETE_ARRAY(codonCount); 
		}
	if( numCodonStates ) {
		MEM_DELETE_ARRAY(numCodonStates); 
		}
	if( codonNumber ) {
        MEM_DELETE_ARRAY(codonNumber); 
        }

	// create new data matrix, and new count and number arrays
	// all counts are initially 1, and characters are numbered
	// sequentially from 0 to nChar-1
	if( taxa > 0 && sites > 0 ) {
		MEM_NEW_ARRAY(codonMatrix,char*,taxa);
		MEM_NEW_ARRAY(codonCount,int,sites);
		MEM_NEW_ARRAY(numCodonStates,int,sites);
		MEM_NEW_ARRAY(codonNumber,int,sites);

		for( int j = 0; j < sites; j++ ) {
			codonCount[j] = 1;
			numCodonStates[j] = 1;
			codonNumber[j] = j;
		}
		for( int i = 0; i < taxa; i++ ) {
			codonMatrix[i]=new char[sites];
			//MEM_NEW_ARRAY(matrix[i],unsigned char,sites);
			//memset( matrix[i], 0xff, taxa*sizeof(unsigned char) );
			memset( codonMatrix[i], 0xff, sites*sizeof(unsigned char) );
		}
	}

	// set dimension variables to new values
	nTax = taxa;
	gapsIncludedNChar = sites;
	}
*/

/*
void CodonData::Pack()
{
	int i, j, newNChar = 0;

	// determine dimensions of new matrix
	for( j = 0; j < nCodonChar; j++ ) {
		if( codonCount[j] )
			newNChar++;
	}

	// create new matrix and count and number arrays and fill
	char** newMatrix;
        MEM_NEW_ARRAY(newMatrix,char*,nTax);
	int* newCount;
        MEM_NEW_ARRAY(newCount,int,newNChar);
	int* newNumStates;
        MEM_NEW_ARRAY(newNumStates,int,newNChar);

	for( i = 0; i < nTax; i++ )
		 MEM_NEW_ARRAY(newMatrix[i],char,newNChar);


	i = 0;
	for( j = 0; j < nCodonChar; j++ ) {
		if( codonCount[j] ) {
			for( int k = 0; k < nTax; k++ )
				newMatrix[k][i] = codonMatrix[k][j];
			newCount[i] = codonCount[j];
			newNumStates[i] = numCodonStates[j];
			//newNumber[i] = number[j];
			i++;
			}
		else{//as we remove columns, shift all the greater numbers over
			for(int c=0;c<gapsIncludedNCodonChar;c++){
				if(codonNumber[c] >= i) codonNumber[c]--;
				}
			}
		}

	// delete old matrix and count and number arrays
	if( codonCount ) MEM_DELETE_ARRAY(codonCount); // count has length nChar
	if( numCodonStates ) MEM_DELETE_ARRAY(numCodonStates); // numStates has length nChar
//	if( number ) MEM_DELETE_ARRAY(number); // number has length nChar
	if( codonMatrix ) {
		for( i = 0; i < nTax; i++ )
			MEM_DELETE_ARRAY(codonMatrix[i]); // matrix[i] has length nChar
		MEM_DELETE_ARRAY(codonMatrix); // matrix has length nTax
        }

	// set count, number and matrix to their new counterparts
	codonCount = newCount;
	numCodonStates = newNumStates;
	codonMatrix = newMatrix;
	nCodonChar = newNChar;	
}
*/

unsigned char NucleotideData::CharToDatum( char ch ){
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

unsigned char AminoacidData::CharToDatum(char ch){
	char datum = 20;
	switch(ch){
		case'A'	: datum=	0	;break;
		case'B' :{
			outman.UserMessage("NOTE: unsupported amino acid abiguity code 'B' changed to full ambiguity");
			datum=	20;
			break;
			}
		case'C' : datum=	1	;break;
		case'D' : datum=	2	;break;
		case'E' : datum=	3	;break;
		case'F' : datum=	4	;break;
		case'G' : datum=	5	;break;
		case'H' : datum=	6	;break;
		case'I' : datum=	7	;break;
		case'K' : datum=	8	;break;
		case'L' : datum=	9	;break;
		case'M' : datum=	10	;break;
		case'N' : datum=	11	;break;
		case'P' : datum=	12	;break;
		case'Q' : datum=	13	;break;
		case'R' : datum=	14	;break;
		case'S' : datum=	15	;break;
		case'T' : datum=	16	;break;
		case'V'	: datum=	17	;break;
		case'W'	: datum=	18	;break;
		case'X' : datum=	20	;break;
		case'Y' : datum=	19	;break;
		case'Z' :{
			outman.UserMessage("NOTE: unsupported amino acid abiguity code 'Z' changed to full ambiguity");
			datum=	20;
			break;
			}
		case'-'	: datum=   20  ;break;
		case'?'	: datum=   20  ;break;
		case'a'	: datum=	0	;break;
		case'b' :{
			outman.UserMessage("NOTE: unsupported amino acid abiguity code 'b' changed to full ambiguity");
			datum=	20;
			break;
			}
		case'c' : datum=	1	;break;
		case'd' : datum=	2	;break;
		case'e' : datum=	3	;break;
		case'f' : datum=	4	;break;
		case'g' : datum=	5	;break;
		case'h' : datum=	6	;break;
		case'i' : datum=	7	;break;
		case'k' : datum=	8	;break;
		case'l' : datum=	9	;break;
		case'm' : datum=	10	;break;
		case'n' : datum=	11	;break;
		case'p' : datum=	12	;break;
		case'q' : datum=	13	;break;
		case'r' : datum=	14	;break;
		case's' : datum=	15	;break;
		case't' : datum=	16	;break;
		case'v'	: datum=	17	;break;
		case'w'	: datum=	18	;break;
		case'x' : datum=	20	;break;
		case'y' : datum=	19	;break;
		case'z' :{
			outman.UserMessage("NOTE: unsupported amino acid abiguity code 'z' changed to full ambiguity");
			datum=	20;
			break;
			}
		default : throw ErrorException("Unknown Amino Acid %c!", ch);
		}
	return datum;
	}

void NucleotideData::CreateMatrixFromNCL(const NxsCharactersBlock *charblock){
	//deprecated - use the 2 param version
	assert(0);

	if(charblock->GetDataType() != NxsCharactersBlock::dna 
		&& charblock->GetDataType() != NxsCharactersBlock::rna 
		&& charblock->GetDataType() != NxsCharactersBlock::nucleotide )
		throw ErrorException("Tried to create nucleotide matrix from non-nucleotide data.\n\t(Check your datatype setting.)");

	if(charblock->GetNumActiveChar() < charblock->GetNChar()){
		NxsUnsignedSet excluded = charblock->GetExcludedIndexSet();
		string exset = NxsSetReader::GetSetAsNexusString(excluded);
		outman.UserMessage("Excluded characters: %s\n\t", exset.c_str());
		}

	int numOrigTaxa = charblock->GetNTax();
	int numActiveTaxa = charblock->GetNumActiveTaxa();
	int numOrigChar = charblock->GetNChar();
	int numActiveChar = charblock->GetNumActiveChar();

	NewMatrix( numActiveTaxa, numActiveChar );

	// read in the data, including taxon names
	int i=0;
	for( int origTaxIndex = 0; origTaxIndex < numOrigTaxa; origTaxIndex++ ) {
		if(charblock->IsActiveTaxon(origTaxIndex)){
			//Now storing names as escaped Nexus values - this means:
			//if they have underscores - store with underscores
			//if they have spaces within single quotes - store with underscores
			//if they have punctuation within single parens (including spaces) - store with single quotes maintained
			NxsString tlabel = charblock->GetTaxonLabel(origTaxIndex);
			SetTaxonLabel(i, NxsString::GetEscaped(tlabel).c_str());
			
			int j = 0;
			for( int origIndex = 0; origIndex < numOrigChar; origIndex++ ) {
				if(charblock->IsActiveChar(origIndex)){	
					unsigned char datum = '\0';
					if(charblock->IsGapState(origTaxIndex, origIndex) == true) datum = 15;
					else if(charblock->IsMissingState(origTaxIndex, origIndex) == true) datum = 15;
					else{
						int nstates = charblock->GetNumStates(origTaxIndex, origIndex);
						for(int s=0;s<nstates;s++){
							datum += CharToBitwiseRepresentation(charblock->GetState(origTaxIndex, origIndex, s));
							}
						}
					SetMatrix( i, j++, datum );
					}
				}
			i++;
			}
		}
	}

void NucleotideData::CreateMatrixFromNCL(const NxsCharactersBlock *charblock, NxsUnsignedSet &charset){
	
	if(charblock->GetDataType() != NxsCharactersBlock::dna 
		&& charblock->GetDataType() != NxsCharactersBlock::rna 
		&& charblock->GetDataType() != NxsCharactersBlock::nucleotide )
		throw ErrorException("Tried to create nucleotide matrix from non-nucleotide data.\n\t(Check your datatype setting.)");

	int numOrigTaxa = charblock->GetNTax();
	int numActiveTaxa = charblock->GetNumActiveTaxa();

	if(charset.empty()){
		//the charset was empty, implying that all characters in this block will go into a single matrix
		for(int i = 0;i < charblock->GetNumChar();i++)
			charset.insert(i);
		}

	//deal with any exclusions
	NxsUnsignedSet excluded = charblock->GetExcludedIndexSet();
	const NxsUnsignedSet *realCharSet = & charset;
	NxsUnsignedSet charsetMinusExcluded;
	if (!excluded.empty()) {
		string exsetName = NxsSetReader::GetSetAsNexusString(excluded);
		outman.UserMessage("Excluded characters: %s\n\t", exsetName.c_str());
		set_difference(charset.begin(), charset.end(), excluded.begin(), excluded.end(), inserter(charsetMinusExcluded, charsetMinusExcluded.begin()));
		realCharSet = &charsetMinusExcluded;
		}

	int numOrigChar = charset.size();
	int numActiveChar = realCharSet->size();

	if(numActiveChar == 0){
		throw ErrorException("Sorry, fully excluded characters blocks or partition subsets are not currently supported.");
		}

	NewMatrix( numActiveTaxa, numActiveChar );

	// read in the data, including taxon names
	int i=0;
	for( int origTaxIndex = 0; origTaxIndex < numOrigTaxa; origTaxIndex++ ) {
		if(charblock->IsActiveTaxon(origTaxIndex)){
			//Now storing names as escaped Nexus values - this means:
			//if they have underscores - store with underscores
			//if they have spaces within single quotes - store with underscores
			//if they have punctuation within single parens (including spaces) - store with single quotes maintained
			NxsString tlabel = charblock->GetTaxonLabel(origTaxIndex);
			SetTaxonLabel(i, NxsString::GetEscaped(tlabel).c_str());
			int j = 0;

			for(NxsUnsignedSet::const_iterator cit = realCharSet->begin(); cit != realCharSet->end();cit++){	
				unsigned char datum = '\0';
				if(charblock->IsGapState(origTaxIndex, *cit) == true) datum = 15;
				else if(charblock->IsMissingState(origTaxIndex, *cit) == true) datum = 15;
				else{
					int nstates = charblock->GetNumStates(origTaxIndex, *cit);
					for(int s=0;s<nstates;s++){
						datum += CharToBitwiseRepresentation(charblock->GetState(origTaxIndex, *cit, s));
						}
					}
				SetMatrix( i, j++, datum );
				}
			i++;
			}
		}
	}

void AminoacidData::CreateMatrixFromNCL(const NxsCharactersBlock *charblock){
	//deprecated - use the 2 param version
	assert(0);

	if(charblock->GetDataType() != NxsCharactersBlock::protein)
		throw ErrorException("Tried to create amino acid matrix from non-amino acid data.\n\t(Did you mean to use datatype = codon-aminoacid?)");

	if(charblock->GetNumActiveChar() < charblock->GetNChar()){
		outman.UserMessageNoCR("Excluded characters:\n\t");
		for(int c=0;c<charblock->GetNCharTotal();c++)
			if(charblock->IsExcluded(c)) outman.UserMessageNoCR("%d ", c+1);
		outman.UserMessage("");
		}
	
	int numOrigTaxa = charblock->GetNTax();
	int numActiveTaxa = charblock->GetNumActiveTaxa();
	int numOrigChar = charblock->GetNChar();
	int numActiveChar = charblock->GetNumActiveChar();

	NewMatrix( numActiveTaxa, numActiveChar );

	// read in the data, including taxon names
	int i=0;
	for( int origTaxIndex = 0; origTaxIndex < numOrigTaxa; origTaxIndex++ ) {
		if(charblock->IsActiveTaxon(origTaxIndex)){
			//Now storing names as escaped Nexus values - this means:
			//if they have underscores - store with underscores
			//if they have spaces within single quotes - store with underscores
			//if they have punctuation within single parens (including spaces) - store with single quotes maintained
			NxsString tlabel = charblock->GetTaxonLabel(origTaxIndex);
			SetTaxonLabel(i, NxsString::GetEscaped(tlabel).c_str());
			
			int j = 0;
			bool firstAmbig = true;
			for( int origIndex = 0; origIndex < numOrigChar; origIndex++ ) {
				if(charblock->IsActiveChar(origIndex)){	
					unsigned char datum = '\0';
					if(charblock->IsGapState(origTaxIndex, origIndex) == true) datum = 20;
					else if(charblock->IsMissingState(origTaxIndex, origIndex) == true) datum = 20;
					else{
						int nstates = charblock->GetNumStates(origTaxIndex, origIndex);
						//assert(nstates == 1);
						//need to deal with the possibility of multiple states represented in matrix
						//just convert to full ambiguity
						if(nstates == 1)
							datum = CharToDatum(charblock->GetState(origTaxIndex, origIndex, 0));
						else{
							if(firstAmbig){
								outman.UserMessageNoCR("Partially ambiguous characters of taxon %s converted to full ambiguity:\n\t", TaxonLabel(origTaxIndex));
								firstAmbig = false;
								}
							outman.UserMessageNoCR("%d ", origIndex+1);
							datum = CharToDatum('?');
							}
						}
					SetMatrix( i, j++, datum );
					}
				}
			if(firstAmbig == false) outman.UserMessage("");
			i++;
			}
		}
	}

void AminoacidData::CreateMatrixFromNCL(const NxsCharactersBlock *charblock, NxsUnsignedSet &charset){

	if(charblock->GetDataType() != NxsCharactersBlock::protein)
		throw ErrorException("Tried to create amino acid matrix from non-amino acid data.\n\t(Did you mean to use datatype = codon-aminoacid?)");

	int numOrigTaxa = charblock->GetNTax();
	int numActiveTaxa = charblock->GetNumActiveTaxa();

	if(charset.empty()){
		//the charset was empty, implying that all characters in this block will go into a single matrix (actually, for nstate
		//might be split anyway).  Create an effective charset that contains all of the characters, which will be filtered
		//for exclusions and for the right number of max states
		for(int i = 0;i < charblock->GetNumChar();i++)
			charset.insert(i);
		}

	//deal with any exclusions
	NxsUnsignedSet excluded = charblock->GetExcludedIndexSet();
	const NxsUnsignedSet *realCharSet = & charset;
	NxsUnsignedSet charsetMinusExcluded;
	if (!excluded.empty()) {
		string exsetName = NxsSetReader::GetSetAsNexusString(excluded);
		outman.UserMessage("Excluded characters: %s\n\t", exsetName.c_str());
		set_difference(charset.begin(), charset.end(), excluded.begin(), excluded.end(), inserter(charsetMinusExcluded, charsetMinusExcluded.begin()));
		realCharSet = &charsetMinusExcluded;
		}

	int numOrigChar = charset.size();
	int numActiveChar = realCharSet->size();

	if(numActiveChar == 0){
		throw ErrorException("Sorry, fully excluded characters blocks or partition subsets are not currently supported.");
		}

	NewMatrix( numActiveTaxa, numActiveChar );

	// read in the data, including taxon names
	int i=0;
	for( int origTaxIndex = 0; origTaxIndex < numOrigTaxa; origTaxIndex++ ) {
		if(charblock->IsActiveTaxon(origTaxIndex)){
			//Now storing names as escaped Nexus values - this means:
			//if they have underscores - store with underscores
			//if they have spaces within single quotes - store with underscores
			//if they have punctuation within single parens (including spaces) - store with single quotes maintained
			NxsString tlabel = charblock->GetTaxonLabel(origTaxIndex);
			SetTaxonLabel(i, NxsString::GetEscaped(tlabel).c_str());
			
			int j = 0;
			bool firstAmbig = true;
//			for( int origIndex = 0; origIndex < numOrigChar; origIndex++ ) {
			for(NxsUnsignedSet::const_iterator cit = realCharSet->begin(); cit != realCharSet->end();cit++){	
				unsigned char datum = '\0';
				if(charblock->IsGapState(origTaxIndex, *cit) == true) datum = 20;
				else if(charblock->IsMissingState(origTaxIndex, *cit) == true) datum = 20;
				else{
					int nstates = charblock->GetNumStates(origTaxIndex, *cit);
					//assert(nstates == 1);
					//need to deal with the possibility of multiple states represented in matrix
					//just convert to full ambiguity
					if(nstates == 1)
						datum = CharToDatum(charblock->GetState(origTaxIndex, *cit, 0));
					else{
						if(firstAmbig){
							outman.UserMessageNoCR("Partially ambiguous characters of taxon %s converted to full ambiguity:\n\t", TaxonLabel(origTaxIndex));
							firstAmbig = false;
							}
						outman.UserMessageNoCR("%d ", *cit+1);
						datum = CharToDatum('?');
						}
					}
				SetMatrix( i, j++, datum );
				}
			if(firstAmbig == false) outman.UserMessage("");
			i++;
			}
		}
	}

