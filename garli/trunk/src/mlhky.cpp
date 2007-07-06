// GARLI version 0.952b2 source code
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






