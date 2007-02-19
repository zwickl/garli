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

#include "defs.h"
#include "bipartition.h"
#include "reconnode.h"

int Bipartition::nBlocks;
int Bipartition:: blockBits;
int Bipartition::ntax;
unsigned int Bipartition::largestBlockDigit;
unsigned int Bipartition::allBitsOn;
char * Bipartition::str;
unsigned int Bipartition::partialBlockMask;

Bipartition::~Bipartition(){
	if(rep!=NULL) delete []rep;
	rep=NULL;
	}	

//note that this is "less than" for sorting purposes, not in a subset sense
bool BipartitionLessThan(const Bipartition &lhs, const Bipartition &rhs){
	int i;
	for(i=0;i<Bipartition::nBlocks-1;i++){
		if(lhs.rep[i] > rhs.rep[i]) return false;
		else if(lhs.rep[i] < rhs.rep[i]) return true;
		}
		
	if(((lhs.rep[i]) & lhs.partialBlockMask) > ((rhs.rep[i]) & lhs.partialBlockMask)) return false;
	else return true;
	}

void Bipartition::SetBipartitionStatics(int nt){
	Bipartition::blockBits=sizeof(int)*8;
	Bipartition::ntax=nt;
	Bipartition::nBlocks=(int)ceil((FLOAT_TYPE)nt/(FLOAT_TYPE)Bipartition::blockBits);
	Bipartition::largestBlockDigit=1<<(Bipartition::blockBits-1);
	Bipartition::allBitsOn=(unsigned int)(pow(2.0, Bipartition::blockBits)-1);
	Bipartition::str=new char[nt+1];
	Bipartition::str[nt] = '\0';
	Bipartition::SetPartialBlockMask();	
	}

void Bipartition::SetPartialBlockMask(){
		partialBlockMask=0;
		unsigned int bit=largestBlockDigit;
		for(int b=0;b<ntax%blockBits;b++){
			partialBlockMask += bit;
			bit = bit >> 1;
			}
		if(ntax%blockBits == 0) partialBlockMask=allBitsOn;
		}