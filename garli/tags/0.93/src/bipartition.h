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

#ifndef BIPARTITION
#define BIPARTITION

#include <cmath>
#include <string>

using namespace std;

class Bipartition{
	public:	
	unsigned int *rep;
	static int nBlocks;
	static int blockBits;
	static int ntax;
	static unsigned int largestBlockDigit;
	static unsigned int allBitsOn;
	static char* str;
	static unsigned int partialBlockMask;//this can be used to mask out the bits that
										// aren't used in the last block.  This becomes
										//important if we start doing complements

	Bipartition(){
		rep=new unsigned int[nBlocks];
		}
	~Bipartition();

	void SetPartialBlockMask(){
		partialBlockMask=0;
		int bit=largestBlockDigit;
		for(int b=0;b<ntax%blockBits;b++){
			partialBlockMask += bit;
			bit = bit >> 1;
			}
		if(ntax%blockBits == 0) partialBlockMask=allBitsOn;
		}
	
	void operator+=(const Bipartition *rhs){
		for(int i=0;i<nBlocks;i++){
			rep[i]|=rhs->rep[i];
			}
		}
	void operator=(const Bipartition *rhs){
		memcpy(rep, rhs->rep, nBlocks*sizeof(int));
		}
	bool EqualsEquals(const Bipartition *rhs){
		for(int i=0;i<nBlocks;i++){
			if(rep[i]!=rhs->rep[i]) return false;
			}
		return true;
		}

	bool EqualsEqualsComplement(const Bipartition *rhs){
		int i;
		for(i=0;i<nBlocks-1;i++){
			if(rep[i]!=~rhs->rep[i]) return false;
			}
			
		if(rep[i]!=(~rhs->rep[i])&partialBlockMask) return false;
		
		return true;
		}

	int FirstPresentTaxon(){
		int blk=0;
		unsigned int tmp=rep[blk];
		while(tmp==0) tmp=rep[++blk];
		
		int t=blockBits*blk+1;
		while( ! (tmp & largestBlockDigit)){
			tmp = tmp << 1;
			t++;
			}
		return t;
		}

	int FirstNonPresentTaxon(){
		int blk=0;
		unsigned int tmp=rep[blk];
		while(tmp==allBitsOn) tmp=rep[++blk];
		
		int t=blockBits*blk+1;
		while((tmp & largestBlockDigit)){
			tmp = tmp << 1;
			t++;
			}
		return t;
		}

	bool ContainsTaxon(int t){
		unsigned int tmp=rep[ t / blockBits];
		if(tmp & largestBlockDigit>>((t-1)%blockBits)) return true;
		return false;
		}
	
	bool IsIncompatible(const Bipartition *target){
		for(int i=0;i<nBlocks;i++){
			if((target->rep[i] | rep[i]) != target->rep[i]) return true;
			}
		return false;
		}

	bool IsIncompatibleComplement(const Bipartition *target){
		int i;
		for(i=0;i<nBlocks-1;i++){
			if((~target->rep[i] | rep[i]) != ~target->rep[i]) return true;
			}
		if(((~target->rep[i])&partialBlockMask | rep[i]) != (~target->rep[i])&partialBlockMask) return true;
		return false;
		}

	Bipartition *TerminalBipart(int taxNum){
		for(int i=0;i<nBlocks;i++) rep[i]=0;
		//8-19-05 fixed this bit of stupidity
		//rep[(taxNum)/blockBits]|=(largestBlockDigit>>((taxNum-1)%blockBits));
		rep[(taxNum-1)/blockBits]|=(largestBlockDigit>>((taxNum-1)%blockBits));
		return this;
		}
	char * Output(){
		for(int i=0;i<nBlocks;i++){
			unsigned int t=rep[i];
			for(int j=0;j<blockBits;j++){
				if(t&largestBlockDigit) str[i*blockBits+j]='*';
				else str[i*blockBits+j]='-';
				t=t<<1;
				}
			}
		return str;
		}
	};

class BipartitionSet{
	Bipartition *allParts;
	};

#endif

