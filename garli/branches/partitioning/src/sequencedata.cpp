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
#include "rng.h"
#include <iterator>

extern rng rnd;
extern OutputManager outman;
extern bool FloatingPointEquals(const FLOAT_TYPE first, const FLOAT_TYPE sec, const FLOAT_TYPE epsilon);

#undef DEBUG_CALCFREQ
#undef DEBUG_CALCPRMATRIX
#undef DEBUGGING_PRMATRICES
#undef DEBUGGING_PATTERN_PROBS

#if defined( DEBUGGING_PATTERN_PROBS )
#	include <io.h>
#endif

//this depends on the fact that a spare taxon was allocated
void SequenceData::AddDummyRootToExistingMatrix(){
	assert(nTaxAllocated > nTax);

	nTax++;
	SetTaxonLabel( nTax - 1, "ROOT");
	for(int c = 0;c < nChar;c++){	
		SetMatrix( nTax - 1, c, maxNumStates);
		}
	}

//this depends on the fact that a spare taxon was allocated
void NucleotideData::AddDummyRootToExistingMatrix(){
	assert(nTaxAllocated > nTax);

	nTax++;
	SetTaxonLabel( nTax - 1, "ROOT");
	for(int c = 0;c < nChar;c++){	
		SetMatrix( nTax - 1, c, 15);
		}
	}

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
		outman.UserMessage("WARNING: Not all amino acids were observed in this dataset.\n\tOne pseudo-count will be added to each amino acid\n\tfor calculation of the empirical frequencies. You\n\tshould probably use a statefrequencies setting other \n\tthan emprical.\n");
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

	string wtsetName = dnaData->WeightsetName();
	if(wtsetName.length() > 0)
		throw ErrorException("Cannot use wtsets with DNA to Codon translation! Remove wtset \"%s\".", wtsetName.c_str());  

	nonZeroCharCount = nChar = dnaData->NChar()/3;
	nTax = dnaData->NTax();
	if(dnaData->NChar() % 3 != 0) throw ErrorException("Codon datatype specified, but number of nucleotides not divisible by 3!");  
	NewMatrix(nTax, nChar);
	patman.Initialize(nTax, maxNumStates);

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
	
	//copy matrix into alternative PatternManager for pattern sorting
	if(usePatternManager){
		Pattern thisPat;
		for(int cod=0;cod<nChar;cod++){
			for(int tax=0;tax<NTax();tax++){
				thisPat.AddChar(matrix[tax][cod]);
				}
			//numbers directly copy over, and should actually be in 0->nchar order now anyway
			thisPat.siteNumbers.push_back(number[cod]);
			thisPat.SetCount(1);
			patman.AddPattern(thisPat);
			thisPat.Reset();
			}
		}
	}

void AminoacidData::FillAminoacidMatrixFromDNA(const NucleotideData *dnaData, GeneticCode *code){
	//first we need to convert the nucleotide data to codons, and then translate the codons to AA's
	//codons are ordered AAA, AAC, AAG, AAT, ACA, ... TTT
	short pos1, pos2, pos3;

	string wtsetName = dnaData->WeightsetName();
	if(wtsetName.length() > 0)
		throw ErrorException("Cannot use wtsets with DNA to Aminoacid translation! Remove wtset \"%s\" or do the translation yourself.", wtsetName.c_str());  

	nonZeroCharCount = nChar = dnaData->NChar()/3;
	nTax = dnaData->NTax();
	if(dnaData->NChar() % 3 != 0) throw ErrorException("Codon to Aminoacid translation specified, but number of nucleotides not divisible by 3!");  
	NewMatrix(nTax, nChar);
	patman.Initialize(nTax, maxNumStates);

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
			//note that a return code of 20 (or 21 for the two serine model) from the codon lookup indicates a stop codon, but a protein code of 20 generally means total ambiguity
			if(thisCodonNum != 64){
				prot = code->CodonLookup(thisCodonNum);
				if(prot == maxNumStates){
					string c;
					char b[4]={'A','C','G','T'};
					c += b[pos1];
					c += b[pos2];
					c += b[pos3];
					throw ErrorException("stop codon %s found at codon site %d (nuc site %d) in taxon %s.  Bailing out.", c.c_str(), cod+1, cod*3+1,  dnaData->TaxonLabel(tax));
					}
				}
			else prot = maxNumStates;

			matrix[tax][cod] = prot;
			}
		if(firstAmbig == false) outman.UserMessage("");
		}	

	//copy matrix into alternative PatternManager for pattern sorting
	if(usePatternManager){
		Pattern thisPat;
		for(int cod=0;cod<nChar;cod++){
			for(int tax=0;tax<NTax();tax++){
				thisPat.AddChar(matrix[tax][cod]);
				}
			//numbers directly copy over, and should actually be in 0->nchar order now anyway
			thisPat.siteNumbers.push_back(number[cod]);
			thisPat.SetCount(1);
			patman.AddPattern(thisPat);
			thisPat.Reset();
			}
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

void NucleotideData::CreateMatrixFromNCL(const NxsCharactersBlock *charblock, NxsUnsignedSet &charset){
	
	if(charblock->GetDataType() != NxsCharactersBlock::dna 
		&& charblock->GetDataType() != NxsCharactersBlock::rna 
		&& charblock->GetDataType() != NxsCharactersBlock::nucleotide )
		throw ErrorException("Tried to create nucleotide matrix from non-nucleotide data.\n\tCheck the datatype settings in your datafile in the characters\n\tor data block and the datatype setting in your Garli config file.");

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
	patman.Initialize(numActiveTaxa, maxNumStates);

	//get weightset if one was specified
	vector<int> charWeights;
	if(useDefaultWeightsets){
		wtsetName = GarliReader::GetDefaultIntWeightSet(charblock, charWeights);
		if(charWeights.size() > 0){
			assert(charWeights.size() == charblock->GetNumChar());
			outman.UserMessage("\tFound wtset \"%s\" with data, applying...", wtsetName.c_str());
			}
		}

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
				if(i == 0)
					SetOriginalDataNumber(j, *cit);
				unsigned char datum = '\0';
				if(charblock->IsGapState(origTaxIndex, *cit) == true) datum = 15;
				else if(charblock->IsMissingState(origTaxIndex, *cit) == true) datum = 15;
				else{
					int nstates = charblock->GetNumStates(origTaxIndex, *cit);
					for(int s=0;s<nstates;s++){
						datum += CharToBitwiseRepresentation(charblock->GetState(origTaxIndex, *cit, s));
						}
					}
				if(i == 0 && charWeights.size() > 0)
					SetCount(j, charWeights[*cit]);
				SetMatrix( i, j, datum );
				j++;
				}
			i++;
			}
		}

	if(usePatternManager){
		//optionally also read into the alternative pattern manager, this is taxa loop inside char loop
		bool haveWeights = !charWeights.empty();
		Pattern thisPat;
		int charNum = 0;
		for(NxsUnsignedSet::const_iterator cit = realCharSet->begin(); cit != realCharSet->end();cit++){	
			int tax = 0;
			for( int origTaxIndex = 0; origTaxIndex < numOrigTaxa; origTaxIndex++ ) {
				if(charblock->IsActiveTaxon(origTaxIndex)){
					unsigned char datum = '\0';
					if(charblock->IsGapState(origTaxIndex, *cit) == true) datum = 15;
					else if(charblock->IsMissingState(origTaxIndex, *cit) == true) datum = 15;
					else{
						int nstates = charblock->GetNumStates(origTaxIndex, *cit);
						for(int s=0;s<nstates;s++){
							datum += CharToBitwiseRepresentation(charblock->GetState(origTaxIndex, *cit, s));
							}
						}
					thisPat.AddChar(datum);
					}
				}
			thisPat.siteNumbers.push_back(charNum++);
			thisPat.SetCount((haveWeights ? charWeights[*cit] : 1));
			patman.AddPattern(thisPat);
			thisPat.Reset();
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
	patman.Initialize(numActiveTaxa, maxNumStates);

	//get weightset if one was specified
	vector<int> charWeights;
	if(useDefaultWeightsets){
		wtsetName = GarliReader::GetDefaultIntWeightSet(charblock, charWeights);
		if(charWeights.size() > 0){
			assert(charWeights.size() == charblock->GetNumChar());
			outman.UserMessage("\tFound wtset \"%s\" with data, applying...", wtsetName.c_str());
			}
		}

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
			for(NxsUnsignedSet::const_iterator cit = realCharSet->begin(); cit != realCharSet->end();cit++){	
				if(i == 0)
					SetOriginalDataNumber(j, *cit);
				unsigned char datum = '\0';
				if(charblock->IsGapState(origTaxIndex, *cit) == true) datum = maxNumStates;
				else if(charblock->IsMissingState(origTaxIndex, *cit) == true) datum = maxNumStates;
				else{
					int nstates = charblock->GetNumStates(origTaxIndex, *cit);
					//need to deal with the possibility of multiple states represented in matrix
					//just convert to full ambiguity
					if(nstates == 1)
						datum = CharToDatum(charblock->GetState(origTaxIndex, *cit, 0));
					else{
						if(firstAmbig){
							outman.UserMessageNoCR("\tPart ambig. char's of taxon %s converted to full ambiguity:\n\t  char ", TaxonLabel(origTaxIndex));
							firstAmbig = false;
							}
						outman.UserMessageNoCR(" %d ", *cit+1);
						datum = CharToDatum('?');
						}
					}
				if(i == 0 && charWeights.size() > 0)
					SetCount(j, charWeights[*cit]);
				SetMatrix( i, j, datum );
				j++;
				}
			if(firstAmbig == false) outman.UserMessage("");
			i++;
			}
		}
	//read the same data into the alternate pattern sorting machinery, which only makes sense looping over tax within char
	if(usePatternManager){
		Pattern thisPat;
		bool haveWeights = !charWeights.empty();
		int charNum = 0;
		for(NxsUnsignedSet::const_iterator cit = realCharSet->begin(); cit != realCharSet->end();cit++){	
			int tax = 0;
			for( int origTaxIndex = 0; origTaxIndex < numOrigTaxa; origTaxIndex++ ) {
				if(charblock->IsActiveTaxon(origTaxIndex)){
					unsigned char datum = '\0';
					if(charblock->IsGapState(origTaxIndex, *cit) == true) datum = maxNumStates;
					else if(charblock->IsMissingState(origTaxIndex, *cit) == true) datum = maxNumStates;
					else{
						int nstates = charblock->GetNumStates(origTaxIndex, *cit);

						//need to deal with the possibility of multiple states represented in matrix
						//just convert to full ambiguity
						if(nstates == 1)
							datum = CharToDatum(charblock->GetState(origTaxIndex, *cit, 0));
						else{
							datum = CharToDatum('?');
							}
						}
					thisPat.AddChar(datum);
					}
				}
			thisPat.siteNumbers.push_back(charNum++);
			thisPat.SetCount((haveWeights ? charWeights[*cit] : 1));
			patman.AddPattern(thisPat);
			thisPat.Reset();
			}
		}
	}
/*
void BinaryData::CreateMatrixFromNCL(const NxsCharactersBlock *charblock, NxsUnsignedSet &origCharset){
	if(charblock->GetDataType() != NxsCharactersBlock::standard)
		throw ErrorException("Tried to create binary matrix from non-standard data.\n\t(Did you mean to use datatype = binary?)");

	//this creates a copy of the charset that we can screw with here without hosing the one that was passed in,
	//which might be needed elsewhere
	NxsUnsignedSet charset = origCharset;

	int numOrigTaxa = charblock->GetNTax();
	int numActiveTaxa = charblock->GetNumActiveTaxa();

	if(charset.empty()){
		//the charset was empty, implying that all characters in this block will go into a single matrix (actually, for nstate
		//might be split anyway).  Create an effective charset that contains all of the characters, which will be filtered
		//for exclusions and for the right number of max states
		for(int i = 0;i < charblock->GetNumIncludedChars();i++)
			charset.insert(i);
		}

	NxsUnsignedSet excluded = charblock->GetExcludedIndexSet();
	const NxsUnsignedSet *realCharSet = & charset;
	NxsUnsignedSet charsetMinusExcluded;
	if (!excluded.empty()) {
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
			//store the taxon names based on NCL's "escaped" version, which will properly deal
			//with whether quotes are necessary, etc.  No conversion needed at output.
			NxsString tlabel = charblock->GetTaxonLabel(origTaxIndex);
			SetTaxonLabel(i, NxsString::GetEscaped(tlabel).c_str());
			
			int j = 0;
			bool firstAmbig = true;
			for(NxsUnsignedSet::const_iterator cit = realCharSet->begin(); cit != realCharSet->end();cit++){
				if(i == 0)
					SetOriginalDataNumber(j, *cit);
				unsigned char datum = '\0';
				if(charblock->IsGapState(origTaxIndex, *cit) == true) datum = 2;
				else if(charblock->IsMissingState(origTaxIndex, *cit) == true) datum = 2;
				else{
					int nstates = charblock->GetNumStates(origTaxIndex, *cit);
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
*/
void NStateData::CreateMatrixFromNCL(const NxsCharactersBlock *charblock, NxsUnsignedSet &origCharset){
	if(charblock->GetDataType() != NxsCharactersBlock::standard)
		throw ErrorException("Tried to create n-state matrix from non-standard data.\n\t(Did you mean to use datatype = standard?)");

	//this creates a copy of the charset that we can screw with here without hosing the one that was passed in,
	//which might be needed elsewhere
	NxsUnsignedSet charset = origCharset;

	int numOrigTaxa = charblock->GetNTax();
	int numActiveTaxa = charblock->GetNumActiveTaxa();

	if(charset.empty()){
		//the charset was empty, implying that all characters in this block will go into a single matrix (actually, for nstate
		//might be split anyway).  Create an effective charset that contains all of the characters, which will be filtered
		//for exclusions and for the right number of max states
		for(int i = 0;i < charblock->GetNumChar();i++)
			charset.insert(i);
		}

	NxsUnsignedSet excluded = charblock->GetExcludedIndexSet();
	NxsUnsignedSet *realCharSet = & charset;
	NxsUnsignedSet charsetMinusExcluded;
	if (!excluded.empty()) {
		set_difference(charset.begin(), charset.end(), excluded.begin(), excluded.end(), inserter(charsetMinusExcluded, charsetMinusExcluded.begin()));
		realCharSet = &charsetMinusExcluded;
	}	

	int numOrigChar = charset.size();
	int numActiveChar = realCharSet->size();

	if(numActiveChar == 0){
		throw ErrorException("Sorry, fully excluded characters blocks or partition subsets are not currently supported.");
		}

	//first count the number of characters with the number of observed states that was specified for
	//this matrix, create a matrix with those dimensions  and grab them from the charblock and make a matrix.
	//If not, just return and the function that called this should be able to check if any characters were actually read, and act accordingly
	//remove_if(realCharSet->begin(), realCharSet->end(), charblock->GetObsNumStates);

	NxsUnsignedSet consts;
	NxsUnsignedSet missing;
	for(NxsUnsignedSet::iterator cit = realCharSet->begin(); cit != realCharSet->end();){
		unsigned num = *cit;
		cit++;
		int ns = charblock->GetObsNumStates(num, false);
		if(ns == 1)
			consts.insert(num);
		//the maxNumStates == 2 part here is so that the message is only output when reading the first standard data matrix
		else if(ns == 0 && maxNumStates == 2)
			missing.insert(num);
		if(datatype == BINARY || datatype == BINARY_NOT_ALL_ZEROS){
			if(ns > 2){
				throw ErrorException("More than two character states found in binary data (character %d)!", num + 1);  
				}
			}
		//In all cases zero-state characters are tossed.  If not binary data, toss any char with #states not equal to the 
		//prespecified maxNumStates for this matrix.  For binary data accept all chars with > 0 states (because any with > 2
		//should already have caused an error above)
		if(ns == 0 || (ns != maxNumStates && !(datatype == BINARY || datatype == BINARY_NOT_ALL_ZEROS))){
			realCharSet->erase(num);
			}
		}
	if(missing.size() > 0){
		string str = NxsSetReader::GetSetAsNexusString(missing);
		outman.UserMessage("\tNOTE: entirely missing characters removed from matrix: %s", str.c_str());
		}

	//verify that we're not breaking the assumptions of these datatypes
	if(consts.size() > 0 && (datatype == ONLY_VARIABLE || datatype == BINARY_NOT_ALL_ZEROS)){
		string c = NxsSetReader::GetSetAsNexusString(consts);
		if(datatype == BINARY_NOT_ALL_ZEROS){
			for(NxsUnsignedSet::iterator cit = consts.begin(); cit != consts.end();){
				int num = *cit;
				cit++;
				std::set<NxsDiscreteStateCell> states = charblock->GetNamedStateSetOfColumn(num);
				assert(states.size() == 1);
				if(states.find(0) == states.end())
					consts.erase(num);
				}
			if(consts.size() > 0){
				string c = NxsSetReader::GetSetAsNexusString(consts);
				throw ErrorException("Constant characters of state 0 are not allowed when using the binarynotallzeros datatype (as opposed to plain binary).\nChange to datatype = binary\n\tor exclude them by adding this to your nexus datafile:\nbegin assumptions;\nexset * const = %s;\nend;", c.c_str());  
				}
			}
		else{
			string c = NxsSetReader::GetSetAsNexusString(consts);
			throw ErrorException("Constant characters are not allowed when using the Mkv\n\tmodel (as opposed to Mk), because it assumes that all\n\tcharacters are variable.  Change to datatype = standard\n\tor exclude them by adding this to your nexus datafile:\nbegin assumptions;\nexset * const = %s;\nend;", c.c_str());
			}
		}
	//maxNumStates = 2 here is only so that the message is output when creating the first standard matrix
	else if(consts.size() > 0 && !(datatype == BINARY) && maxNumStates == 2){
		string c = NxsSetReader::GetSetAsNexusString(consts);
		outman.UserMessage("\t****\n\tWARNING - Constant characters found in standard data matrix (sites %s)", c.c_str());
		outman.UserMessage("\tCurrently these will be ignored because including them in the likelihood");
		outman.UserMessage("\tcalculations would require knowledge of how many states were possible for");
		outman.UserMessage("\tthose columns (i.e., 1 state was observed, but was that out of 2 possible,");
		outman.UserMessage("\tor 3 or 4, etc)\n\t****");
		}

	if(realCharSet->size() == 0)
		return;

	//Make room for dummy conditioning (generally constant) character(s) here.
	//For anything besides BINARY_NOT_ALL_ZEROS the # will be equal to maxNumStates
	//although for symetrical Mkv that many are not needed because they are all equal
	//it defaults to zero in the constructor
	if(datatype == ONLY_VARIABLE || datatype == BINARY_NOT_ALL_ZEROS){
		if(datatype == BINARY_NOT_ALL_ZEROS)
			numConditioningPatterns = 1;
		else
			numConditioningPatterns = maxNumStates;
		}

	NewMatrix( numActiveTaxa, realCharSet->size() + numConditioningPatterns);

	map<int, int> nclStateIndexToGarliState;
	vector< map<int, int> > stateMaps;

	bool recodeSkipped = false;
	if(modeltype == UNORDERED && !(datatype == BINARY || datatype == BINARY_NOT_ALL_ZEROS))
		recodeSkipped = true;

	if(recodeSkipped){
		//Recode characters that skip states (assuming numerical order of states) to not skip any.  i.e., recode a
		//char with states 0 1 5 7 to 0 1 2 3 and assume that it has 4 states
		//With assumptions block "options gapmode=newstate" things get even more confusing.  GetNamedStateSetOfColumn
		//returns the gap as a code of -2, in which case the mapping would be -2 0 1 5 7 -> 0 1 2 3 4 5 
		if(datatype == ONLY_VARIABLE)
			//add in the conditioning patterns such that the character numbers match up later
			for(int i = 0;i < numConditioningPatterns;i++)
				stateMaps.push_back(nclStateIndexToGarliState);

		for(NxsUnsignedSet::const_iterator cit = realCharSet->begin(); cit != realCharSet->end();cit++){
			set<int> stateSet = charblock->GetNamedStateSetOfColumn(*cit);
			int myIndex = 0;
			for(set<int>::iterator sit = stateSet.begin();sit != stateSet.end();sit++){
				nclStateIndexToGarliState.insert(pair<int, int>(*sit, myIndex++));
				}
			stateMaps.push_back(nclStateIndexToGarliState);
			nclStateIndexToGarliState.clear();
			}
		}
	else{//for ordered data we don't want to remove unobserved states
		if(charblock->GetGapModeSetting() == CharactersBlock::GAP_MODE_NEWSTATE){
			throw ErrorException("Cannot use ordered Mk/Mkv data with gapmode=newstate. Recode the state or choose unordered.");
			}
		}

	// read in the data, including taxon names
	int effectiveTax=0;
	for( int origTaxIndex = 0; origTaxIndex < numOrigTaxa; origTaxIndex++ ) {
		if(charblock->IsActiveTaxon(origTaxIndex)){
			//store the taxon names based on NCL's "escaped" version, which will properly deal
			//with whether quotes are necessary, etc.  No conversion needed at output.
			NxsString tlabel = charblock->GetTaxonLabel(origTaxIndex);
			SetTaxonLabel(effectiveTax, NxsString::GetEscaped(tlabel).c_str());
			
			int effectiveChar = 0;
			//add the dummy constant character(s).  This will be one of each possible constant
			//state, except for BINARY_NOT_ALL_ZEROS, where it will be only state 0
			if(numConditioningPatterns > 0){
				for(int s = 0; s < (datatype == BINARY_NOT_ALL_ZEROS ? 1 : maxNumStates); s ++){
					if(effectiveTax == 0)
						SetOriginalDataNumber(s, -1);
					SetMatrix( effectiveTax, effectiveChar++, s);
					}
				}

			bool firstAmbig = true;
			for(NxsUnsignedSet::const_iterator cit = realCharSet->begin(); cit != realCharSet->end();cit++){
				if(effectiveTax == 0)
					SetOriginalDataNumber(effectiveChar, *cit);
				unsigned char datum = '\0';
				if(charblock->IsGapState(origTaxIndex, *cit) == true){
					if(datatype == BINARY || datatype == BINARY_NOT_ALL_ZEROS)
						throw ErrorException("Cannot use gap characters with binary datatype.  Recode to 0 and 1");
					//if gapmode=newstate is on (default is gapmode=missing) then need handle the gap properly
					//changes in NCL should now have it correctly reporting the number of states with gaps {in, ex}cluded
					if(charblock->GetGapModeSetting() == CharactersBlock::GAP_MODE_NEWSTATE){
						if(recodeSkipped){
							datum = stateMaps[effectiveChar][NXS_GAP_STATE_CODE];
							}
						else{
							assert(0);
							datum = maxNumStates - 1;
							}
						}
					else{
						datum = maxNumStates;
						}
					}
				else if(charblock->IsMissingState(origTaxIndex, *cit) == true){
					datum = maxNumStates;
					}
				else{
					int nstates = charblock->GetNumStates(origTaxIndex, *cit);
					//need to deal with the possibility of multiple states represented in matrix
					//just convert to full ambiguity
					if(nstates == 1){
						int nclIndex = charblock->GetStateIndex(origTaxIndex, *cit, 0);
						if(recodeSkipped)
							datum = stateMaps[effectiveChar][nclIndex];
						else 
							datum = nclIndex;
						}
					else{
						if(firstAmbig){
							outman.UserMessageNoCR("\tPart ambig. char's of taxon %s converted to full ambiguity:\n\t  char ", TaxonLabel(origTaxIndex));
							firstAmbig = false;
							}
						outman.UserMessageNoCR(" %d ", *cit+1);
						datum = maxNumStates;
						}
					}
				SetMatrix( effectiveTax, effectiveChar++, datum);
				}
			if(firstAmbig == false) outman.UserMessage("");
			effectiveTax++;
			}
		}
	//verify that every allowed state was observed for each character
#ifndef NDEBUG
	bool found;
	if(recodeSkipped){
		for(int c = numConditioningPatterns;c < nChar;c++){
			for(int s = 0;s < maxNumStates;s++){
				found = false;
				for(int t = 0;t < nTax;t++){
					if(Matrix(t, c) == s){
						found = true;
						break;
						}
					}
				if(!found){
					outman.UserMessage("\nWARNING - some state in a %d-state character appeared only as part\n\tof an ambiguity code, e.g., a column with states 0, 1 and (12).", maxNumStates);
					outman.UserMessage("\tThe ambiguity code will be treated as missing data,\n\tbut the character will still be considered to have %d states.\n", maxNumStates);
					}
				}
			}
		}
#endif
	}

//this is a virtual overload for NState because it might have to deal with the dummy char, which shouldn't be included in the resampling
int NStateData::BootstrapReweight(int seedToUse, FLOAT_TYPE resampleProportion){
	//a seed is passed in and used for the reweighting - Either for restarting or not
	//Either way we'll return the seed at the end of the reweighting, to be stored as the Population::nextBootstrapSeed

	//which allows exactly the same bootstraped datasets to be used in multiple runs, but with different
	//settings for the actual search
	if(resampleProportion >= 5.0) outman.UserMessage("WARNING: The resampleproportion setting is the proportion to resample,\nNOT the percentage (1.0 = 100%%).\nThe value you specified (%.2f) is a very large proportion.", resampleProportion);

	int originalSeed = rnd.seed();
	rnd.set_seed(seedToUse);

	//for nstate data this will include the conditioning chars, but they will
	//have a resample prob of zero
	FLOAT_TYPE *cumProbs = new FLOAT_TYPE[nChar];
	
	assert(origCounts[0] > 0 && origCounts[1] > 0);

	for(int i = 0;i < numConditioningPatterns;i++)
		cumProbs[i] = ZERO_POINT_ZERO;
	cumProbs[numConditioningPatterns]=(FLOAT_TYPE) origCounts[numConditioningPatterns] / ((FLOAT_TYPE) totalNChar - numConditioningPatterns);
	for(int i=numConditioningPatterns + 1;i<nChar;i++){
		cumProbs[i] = cumProbs[i-1] + (FLOAT_TYPE) origCounts[i] / ((FLOAT_TYPE) totalNChar - numConditioningPatterns);
		assert(origCounts[i] > 0);
		}
	for(int q=numConditioningPatterns;q<nChar;q++) 
		count[q]=0;
	assert(FloatingPointEquals(cumProbs[nChar-1], ONE_POINT_ZERO, 1e-6));
	cumProbs[nChar-1] = ONE_POINT_ZERO;
		
	//ofstream deb("counts.log", ios::app);

	//round to nearest int
	int numToSample = (int) (((FLOAT_TYPE)totalNChar * resampleProportion) + 0.5);
	if(numToSample != totalNChar) outman.UserMessage("Resampling %d characters (%.2f%%).\n", numToSample, resampleProportion*100);

	for(int c=0;c<numToSample;c++){
		FLOAT_TYPE p=rnd.uniform();
		int pat=0; 
		while(p > cumProbs[pat]) pat++;
		count[pat]++;
		}
/*
	for(int i = 0;i < nChar;i++)
		deb << i << "\t" << cumProbs[i] << "\t" << origCounts[i] << "\t" << count[i] <<  endl;
*/
	//take a count of the number of chars that were actually resampled
	nonZeroCharCount = 0;
	int numZero = 0;
	int totCounts = 0;
	for(int d = numConditioningPatterns;d < nChar;d++){
		if(count[d] > 0) {
			nonZeroCharCount++;
			totCounts += count[d];
			}
		else 
			numZero++;
		}
	if(datatype == ONLY_VARIABLE) 
		assert(count[0] == 1);
	assert(totCounts == totalNChar);
	assert(nonZeroCharCount + numZero == nChar - numConditioningPatterns);

	delete []cumProbs;
	int nextSeed = rnd.seed();
	rnd.set_seed(originalSeed);
	return nextSeed;
	}

void OrientedGapData::CreateMatrixFromNCL(const NxsCharactersBlock *charblock, NxsUnsignedSet &origCharset){
	if(charblock->GetDataType() != NxsCharactersBlock::standard)
		throw ErrorException("Tried to create n-state matrix from non-standard data.\n\t(Did you mean to use datatype = nstate?)");

	//this creates a copy of the charset that we can screw with here without hosing the one that was passed in,
	//which might be needed elsewhere
	NxsUnsignedSet charset = origCharset;

	int numOrigTaxa = charblock->GetNTax();
	int numActiveTaxa = charblock->GetNumActiveTaxa();

	if(charset.empty()){
		//the charset was empty, implying that all characters in this block will go into a single matrix (actually, for nstate
		//might be split anyway).  Create an effective charset that contains all of the characters, which will be filtered
		//for exclusions and for the right number of max states
		for(int i = 0;i < charblock->GetNumChar();i++)
			charset.insert(i);
		}

	NxsUnsignedSet excluded = charblock->GetExcludedIndexSet();
	NxsUnsignedSet *realCharSet = & charset;
	NxsUnsignedSet charsetMinusExcluded;
	if (!excluded.empty()) {
		set_difference(charset.begin(), charset.end(), excluded.begin(), excluded.end(), inserter(charsetMinusExcluded, charsetMinusExcluded.begin()));
		realCharSet = &charsetMinusExcluded;
	}	

	int numOrigChar = charset.size();
	int numActiveChar = realCharSet->size();

	if(numActiveChar == 0){
		throw ErrorException("Sorry, fully excluded characters blocks or partition subsets are not currently supported.");
		}

	if(realCharSet->size() == 0)
		return;

//	the dummy root is now taken care of outside of here in a non-datatype specific way
//	int myEffectiveTaxa = numActiveTaxa + 1;

	bool allGapChar = true;

	//Make room for dummy conditioning (all zero) character here.
	//it defaults to zero in the constructor
	if(datatype == ONLY_VARIABLE || allGapChar){
		numConditioningPatterns = 1;
		}

	//make room for a dummy constant character here
	NewMatrix( numActiveTaxa, realCharSet->size() + numConditioningPatterns);

	// read in the data, including taxon names
	int effectiveTax=0;
	for( int origTaxIndex = 0; origTaxIndex < numOrigTaxa; origTaxIndex++ ) {
		if(charblock->IsActiveTaxon(origTaxIndex)){
			//store the taxon names based on NCL's "escaped" version, which will properly deal
			//with whether quotes are necessary, etc.  No conversion needed at output.
			NxsString tlabel = charblock->GetTaxonLabel(origTaxIndex);
			SetTaxonLabel(effectiveTax, NxsString::GetEscaped(tlabel).c_str());
			
			int effectiveChar = 0;
			//add the dummy character
			if(numConditioningPatterns > 0){
				if(effectiveTax == 0)
					SetOriginalDataNumber(0, -1);
				if(tlabel != "ROOT")
					SetMatrix( effectiveTax, effectiveChar++, 0);
				else
					SetMatrix( effectiveTax, effectiveChar++, maxNumStates);
				}

			bool firstAmbig = true;
			for(NxsUnsignedSet::const_iterator cit = realCharSet->begin(); cit != realCharSet->end();cit++){	
				if(effectiveTax == 0)
					SetOriginalDataNumber(effectiveChar, *cit);
				unsigned char datum = '\0';
				if(charblock->IsGapState(origTaxIndex, *cit) == true)
					datum = 0;
				else if(charblock->IsMissingState(origTaxIndex, *cit) == true){
					datum = maxNumStates;
					}
				else{
					int nstates = charblock->GetNumStates(origTaxIndex, *cit);
					if(nstates == 1){
						int nclIndex = charblock->GetStateIndex(origTaxIndex, *cit, 0);
						datum = nclIndex;
						}
					else{
						if(firstAmbig){
							outman.UserMessageNoCR("\tPart ambig. char's of taxon %s converted to full ambiguity:\n\t  char ", TaxonLabel(origTaxIndex));
							firstAmbig = false;
							}
						outman.UserMessageNoCR(" %d ", *cit+1);
						datum = maxNumStates;
						}
					}
				SetMatrix( effectiveTax, effectiveChar++, datum);
				}
			if(firstAmbig == false) outman.UserMessage("");
			effectiveTax++;
			}
		}
	}
