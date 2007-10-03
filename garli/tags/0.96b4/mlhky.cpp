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

#include "defs.h"
#include "mlhky.h"

#undef DEBUG_CALCFREQ
#undef DEBUG_CALCPRMATRIX
#undef DEBUGGING_PRMATRICES
#undef DEBUGGING_PATTERN_PROBS

#if defined( DEBUGGING_PATTERN_PROBS )
#	include <io.h>
#endif

void HKYData::CalcEmpiricalFreqs(){
	if(maxNumStates == 20){
		CalcEmpiricalAAFreqs();
		return;
		}
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
				if(fabs(tempFreqs[j] - empStateFreqs[j]) > 1.0e-8) continueIterations = true;
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

void HKYData::CalcEmpiricalAAFreqs(){
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
		for(int j=0;j<maxNumStates;j++) empStateFreqs[j] += ONE_POINT_ZERO;
		total += (FLOAT_TYPE) maxNumStates;
		}
	for(int j=0;j<maxNumStates;j++){
		empStateFreqs[j] /= total;
		freqTot += empStateFreqs[j];
		}
	assert(fabs(freqTot - 1.0) < 1e-5);
	}


void HKYData::MakeAmbigStrings(){
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

void CodonData::FillCodonMatrix(bool translateAminoAcids){
	//first we need to convert the nucleotide data to codons numbered 0-61 and assign them back to the terminals
	//codons are ordered AAA, AAC, AAG, AAT, ACA, ... TTT
	short pos1, pos2, pos3;

	nCodonChar = NNucChar()/3;
	if(NNucChar() % 3 != 0) throw ErrorException("Codon datatype specified, but number of nucleotides not divisible by 3!");  
	NewCodonMatrix(NTax(), nCodonChar);

	//this will just map from the bitwise format to the index format (A, C, G, T = 0, 1, 2, 3)
	//partial ambiguity is mapped to total ambiguity currently
	short bitwiseToIndexFormat[16] = {15,0,1,15,2,15,15,15,3,15,15,15,15,15,15,15};
	
	int tax=0, thisCodonNum;
	//for(termit=termTax->begin();termit!=termTax->end();termit++){
	for(int tax=0;tax<NTax();tax++){
		//nucseq=(*termit)->characters->rawData;
		for(int cod=0;cod<nCodonChar;cod++){
			short p1 = NucMatrix(tax, cod*3);
			short p2 = NucMatrix(tax, cod*3+1);
			short p3 = NucMatrix(tax, cod*3+2);

			pos1 = bitwiseToIndexFormat[p1];
			pos2 = bitwiseToIndexFormat[p2];
			pos3 = bitwiseToIndexFormat[p3];
			
			thisCodonNum=(pos1)*16 + (pos2)*4 + pos3;
			
/*			if(thisCodonNum==48||thisCodonNum==50||thisCodonNum==56){//check for stop codons
				throw ErrorException("stop codon found at codon site %d in taxon %s.  Bailing out.", cod,  TaxonLabel(tax));
				}
*/
			if(pos1==15||pos2==15||pos3==15){//check for gaps
				if(pos1+pos2+pos3 != 45){
					//warn about gaps or ambiguity in codons
					outman.UserMessage("Warning!: gaps or ambiguity codes found in codon site %d for taxon %s", cod, TaxonLabel(tax));
					outman.UserMessage("\tThis codon coded as missing (ambiguity not currently supported).");
					}
				thisCodonNum=64;
				}

			char prot;
			//note that a return code of 20 from the codon lookup indicates a stop codon, but a protein code of 20 generally means total ambiguity
			if(thisCodonNum != 64){
				prot = code.CodonLookup(thisCodonNum);
				if(prot == 20)
					throw ErrorException("stop codon found at codon site %d in taxon %s.  Bailing out.", cod,  TaxonLabel(tax));
				}
			else prot = 20;

			if(translateAminoAcids){
				codonMatrix[tax][cod] = prot;
				}
			else{
				if(thisCodonNum == 64)//missing or ambiguous 
					codonMatrix[tax][cod] = 61;
				else 
					codonMatrix[tax][cod] = code.Map64stateTo61state(thisCodonNum);
				}
			}
		}	
	}

void CodonData::CalcEmpiricalFreqs(){
	empStateFreqs=new FLOAT_TYPE[maxNumStates];//this is a member of the class, and where the final freqs will be stored
	
	for(int i=0;i<maxNumStates;i++) empStateFreqs[i] = 0.0;
	FLOAT_TYPE total = 0.0;
	//for codons and aminoacids this will assume no ambiguity
	
	for( int i = 0; i < NTax(); i++ ) {
		for( int j = 0; j < nCodonChar; j++ ) {
			char thischar= codonMatrix[i][j];

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
		for(int j=0;j<maxNumStates;j++) empStateFreqs[j] += ONE_POINT_ZERO;
		total += (FLOAT_TYPE) maxNumStates;
		}
	for(int j=0;j<maxNumStates;j++){
		empStateFreqs[j] /= total;
		freqTot += empStateFreqs[j];
		}
	assert(fabs(freqTot - 1.0) < 1e-5);
	}

void CodonData::CalcF1x4Freqs(){
	
	empBaseFreqPos1=new FLOAT_TYPE[4];
	empBaseFreqPos1[0]=empBaseFreqPos1[1]=empBaseFreqPos1[2]=empBaseFreqPos1[3]=0.0;
	
	FLOAT_TYPE total = ZERO_POINT_ZERO;
	for( int i = 0; i < NTax(); i++ ) {
		for( int j = 0; j < NNucChar(); j++ ) {
			char thischar=(char) NucMatrix( i, j );

			if(thischar & 1)
				empBaseFreqPos1[0] += (FLOAT_TYPE) Count(j);
			if(thischar & 2)
				empBaseFreqPos1[1] += (FLOAT_TYPE) Count(j);
			if(thischar & 4)
				empBaseFreqPos1[2] += (FLOAT_TYPE) Count(j);
			if(thischar & 8) 
				empBaseFreqPos1[3] += (FLOAT_TYPE) Count(j);	
			total += Count(j);
			}
		}
	
	for(int j=0;j<4;j++)
		empBaseFreqPos1[j] /= (total);
		
	int stops=0;
	for(int base1=0;base1<4;base1++){
		for(int base2=0;base2<4;base2++){
			for(int base3=0;base3<4;base3++){
				if(code.CodonLookup(base1*16+base2*4+base3) != 20) empStateFreqs[base1*16+base2*4+base3 - stops] = empBaseFreqPos1[base1] * empBaseFreqPos1[base2] * empBaseFreqPos1[base3];
				else stops++;
				}
			}
		}
	}

void CodonData::CalcF3x4Freqs(){
	
	empBaseFreqPos1=new FLOAT_TYPE[4];
	empBaseFreqPos1[0]=empBaseFreqPos1[1]=empBaseFreqPos1[2]=empBaseFreqPos1[3]=0.0;
	empBaseFreqPos2=new FLOAT_TYPE[4];
	empBaseFreqPos2[0]=empBaseFreqPos2[1]=empBaseFreqPos2[2]=empBaseFreqPos2[3]=0.0;
	empBaseFreqPos3=new FLOAT_TYPE[4];
	empBaseFreqPos3[0]=empBaseFreqPos3[1]=empBaseFreqPos3[2]=empBaseFreqPos3[3]=0.0;

	FLOAT_TYPE total1 = ZERO_POINT_ZERO;
	FLOAT_TYPE total2 = ZERO_POINT_ZERO;
	FLOAT_TYPE total3 = ZERO_POINT_ZERO;
	
	for( int i = 0; i < NTax(); i++ ) {
		for( int j = 0; j < NNucChar(); j++ ) {
			char thischar=(char) NucMatrix( i, j );
			if(j%3 == 1){
				if(thischar & 1)
					empBaseFreqPos1[0] += (FLOAT_TYPE) Count(j);
				if(thischar & 2)
					empBaseFreqPos1[1] += (FLOAT_TYPE) Count(j);
				if(thischar & 4)
					empBaseFreqPos1[2] += (FLOAT_TYPE) Count(j);
				if(thischar & 8) 
					empBaseFreqPos1[3] += (FLOAT_TYPE) Count(j);
				total1 += Count(j);
				}
			if(j%3 == 2){
				if(thischar & 1)
					empBaseFreqPos2[0] += (FLOAT_TYPE) Count(j);
				if(thischar & 2)
					empBaseFreqPos2[1] += (FLOAT_TYPE) Count(j);
				if(thischar & 4)
					empBaseFreqPos2[2] += (FLOAT_TYPE) Count(j);
				if(thischar & 8) 
					empBaseFreqPos2[3] += (FLOAT_TYPE) Count(j);
				total2 += Count(j);
				}					
			if(j%3 == 3){
				if(thischar & 1)
					empBaseFreqPos3[0] += (FLOAT_TYPE) Count(j);
				if(thischar & 2)
					empBaseFreqPos3[1] += (FLOAT_TYPE) Count(j);
				if(thischar & 4)
					empBaseFreqPos3[2] += (FLOAT_TYPE) Count(j);
				if(thischar & 8) 
					empBaseFreqPos3[3] += (FLOAT_TYPE) Count(j);
				total3 += Count(j);
				}					
			}
		}
	
	for(int j=0;j<4;j++){
		empBaseFreqPos1[j] /= (total1);
		empBaseFreqPos2[j] /= (total2);
		empBaseFreqPos3[j] /= (total3);
		}
		
	int stops=0;
	for(int base1=0;base1<4;base1++){
		for(int base2=0;base2<4;base2++){
			for(int base3=0;base3<4;base3++){
				if(code.CodonLookup(base1*16+base2*4+base3) != 20) empStateFreqs[base1*16+base2*4+base3 - stops] = empBaseFreqPos1[base1] * empBaseFreqPos2[base2] * empBaseFreqPos3[base3];
				else stops++;
				}
			}
		}
	}

//
// ComparePatterns returns:
//	 0		complete identity
//	-1		if i less than j
//	 1		if i greater than j
//
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

// PatternType determines whether pattern k is constant, informative, or autoapomorphic
//
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
	if( nStates == 1 /*|| (nStates==2 && missing)*/)
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
		int max = maxNumStates;
		for( int k = 0; k <= max; k++ )
			stateDistr[k] = 0.0;
		}

	// set dimension variables to new values
	nTax = taxa;
	gapsIncludedNChar = sites;
	}

//
// Collapse merges like patterns
//
void CodonData::Collapse(){
	int i = 0, j = 1;

	Sort();

	while( i < NChar() ) {
		while( j < NChar() && ComparePatterns( i, j ) == 0 ) {
			// pattern j same as pattern i
			codonCount[i] += codonCount[j];
			codonCount[j] = 0;
			j++;
			}
		i = j++;
		}
		
	//DJZ 10/28/03 get rid of all missing patterns	
	int q=nCodonChar-1;
	gapsIncludedNCodonChar = nCodonChar;
	while(numCodonStates[q]==0){
		for(i=0;i<gapsIncludedNCodonChar;i++) if(codonNumber[i]==q) codonNumber[i]=-1;
		codonCount[q--]=0;
		//when all missing columns are deleted, remove them from the total number of characters
		nCodonChar--;
		}
	
	Pack();
	}

void CodonData::SwapCharacters( int i, int j )
{
	unsigned char tmp;
	for( int k = 0; k < nTax; k++ ) {
		tmp = Matrix( k, i );
		SetMatrix( k, i, Matrix( k, j ) );
		SetMatrix( k, j, tmp );
		
	}
	//DJZ also swap the nStates array
	int s=numCodonStates[i];
	numCodonStates[i]=numCodonStates[j];
	numCodonStates[j]=s;

	//DJZ 2/14/06 and the number array
	for(int c=0;c<gapsIncludedNChar;c++){
		if(codonNumber[c] == i) codonNumber[c]=j;
		else if(codonNumber[c] == j) codonNumber[c]=i;
		}
}


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
